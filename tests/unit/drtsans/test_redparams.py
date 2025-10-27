# standard library imports
from copy import deepcopy
import json
import os
from os.path import join as path_join
from pathlib import Path
import shutil
import tempfile

# third-party imports
from mantid.kernel import amend_config
import pytest

# drtsans imports
from drtsans import configdir
from drtsans.instruments import instrument_standard_names, instrument_standard_name
from drtsans.redparams import (
    DefaultJson,
    generate_json_files,
    load_schema,
    ReductionParameterError,
    ReductionParameters,
    ReferenceResolver,
    reduction_parameters,
    resolver_common,
    Suggest,
    update_reduction_parameters,
    validate_reduction_parameters,
)


@pytest.fixture(scope="module")
def redparams_data():
    schema_common = json.loads(
        r"""
    {
      "$schema": "http://json-schema.org/draft-07/schema#",

      "definitions": {
        "runNumberTypes": {
          "anyOf": [{"type": "string", "minLength": 1, "pattern": "^[1-9][0-9]*$"},
                    {"type": "integer", "minimum": 1},
                    {"type": "array", "items": {"type": ["integer", "string"]}}],
          "preferredType": "int"
          },
        "transmissionValueTypes": {
          "anyOf": [{"type": "string", "pattern": "^$|^0?.[0-9]*$"},
                    {"type": "number", "exclusiveMinimum": 0, "maximum": 1}],
                    "preferredType": "float"
        },
        "safeStringPositiveFloat": {
          "anyOf": [{"type": "string", "pattern": "^$|^[0-9]*.[0-9]*$"},
                    {"type": "number", "exclusiveMinimum": 0}],
          "preferredType": "float"
        }
      },

      "instrumentName": {
        "type": "string",
        "description": "The name of the instrument. Valid values are BIOSANS, EQSANS, and GPSANS",
        "enum": ["BIOSANS", "EQSANS", "GPSANS"],
        "examples": ["BIOSANS", "EQSANS", "GPSANS"]
      },
      "iptsNumber": {
        "anyOf": [{"type": "string", "minLength": 1, "pattern": "^[1-9][0-9]*$"},
                  {"type": "integer", "minimum": 1}],
        "preferredType": "int",
        "description": "The IPTS number for the data files",
        "examples": ["24769"]
      },

      "sample": {
        "type": "object",
        "properties": {
          "runNumber": {
            "$ref": "common.json#/definitions/runNumberTypes",
            "description": "The run number(s) for the sample"
          },
          "transmission": {
            "type": "object",
            "properties": {
              "runNumber": {
                "$ref": "common.json#/definitions/runNumberTypes",
                "description": "The run number(s) for the transmission sample."
              },
              "value": {"$ref": "common.json#/definitions/transmissionValueTypes"}
            },
            "maxProperties": 2,
            "required": ["runNumber", "value"],
            "description": "The transmission for the sample"
          }
        },
        "maxProperties": 2,
        "required": ["transmission", "runNumber"]
      },

      "configuration": {
        "outputDir": {"type": "string", "description": "Output folder"},
        "useTimeSlice": {
          "type": "boolean",
          "useentry": "timeSliceInterval",
          "description": "Indicate whether the data should be processed as time slices"
        },
        "timeSliceInterval": {
          "$ref": "common.json#/definitions/safeStringPositiveFloat",
          "description": "Interval for time slicing"
        }
      }
    }"""
    )

    schema_common_dir = tempfile.mkdtemp(dir="/tmp")
    schema_common_file = os.path.join(schema_common_dir, "common.json")
    with open(schema_common_file, "w") as file_handle:
        json.dump(schema_common, file_handle)

    schema_instrument = json.loads(
        r"""
    {
      "$schema": "http://json-schema.org/draft-07/schema#",
      "type": "object",
      "properties": {
        "instrumentName": {"$ref": "common.json#/instrumentName", "default": "BIOSANS"},
        "iptsNumber": {"$ref": "common.json#/iptsNumber", "default": 42},
        "sample":     {"$ref": "common.json#/sample"},

        "configuration": {
          "type": "object",
          "properties": {
            "outputDir": {"$ref":  "common.json#/configuration/outputDir"},
            "timeSliceInterval": {"$ref":  "common.json#/configuration/timeSliceInterval"}
          },
          "maxProperties": 2,
          "required": ["outputDir", "timeSliceInterval"]
        }

      },
      "maxProperties": 4,
      "required": ["instrumentName", "iptsNumber", "sample", "configuration"]
    }"""
    )

    reduction_parameters = json.loads(
        r"""
    {"instrumentName": "BIOSANS",
     "iptsNumber": 42,
      "sample": {"runNumber": 24, "transmission": {"value": 0.95}},
      "configuration": {"outputDir": "/tmp", "timeSliceInterval": 20}
    }"""
    )

    yield {
        "schema_common": schema_common,
        "schema_common_file": schema_common_file,
        "schema_instrument": schema_instrument,
        "reduction_parameters": reduction_parameters,
    }
    shutil.rmtree(schema_common_dir)  # clean up before leaving module's scope


@pytest.mark.parametrize("instrument_name", instrument_standard_names())
def test_load_schema(instrument_name):
    r"""correct loading of the schema for each instrument"""
    schema = load_schema(instrument_name)
    assert instrument_standard_name(schema["properties"]["instrumentName"]["default"]) == instrument_name


class TestReferenceResolver:
    def test_derefence(self, redparams_data):
        resolver = ReferenceResolver(redparams_data["schema_common_file"])
        unresolved = {
            "background": {
                "iptsNumber": {"$ref": "common.json#/iptsNumber"},
                "default": 12345,
            }
        }
        resolved = resolver.dereference(unresolved)
        compared = {
            "background": {
                "iptsNumber": {
                    "anyOf": [
                        {"type": "string", "minLength": 1, "pattern": "^[1-9][0-9]*$"},
                        {"type": "integer", "minimum": 1},
                    ],
                    "preferredType": "int",
                    "description": "The IPTS number for the data files",
                    "examples": ["24769"],
                },
                "default": 12345,
            }
        }
        assert resolved == compared
        unresolved = {
            "configuration": {
                "outputDir": {"type": "string", "description": "Output folder"},
                "timeSliceInterval": {
                    "$ref": "common.json#/definitions/safeStringPositiveFloat",
                    "description": "Interval for time slicing",
                },
            }
        }
        resolved = resolver.dereference(unresolved)
        compared = {
            "configuration": {
                "outputDir": {"type": "string", "description": "Output folder"},
                "timeSliceInterval": {
                    "description": "Interval for time slicing",
                    "preferredType": "float",
                    "anyOf": [
                        {"type": "string", "pattern": "^$|^[0-9]*.[0-9]*$"},
                        {"type": "number", "exclusiveMinimum": 0},
                    ],
                },
            }
        }
        assert resolved == compared
        unresolved = {
            "$schema": "http://json-schema.org/draft-07/schema#",
            "type": "object",
            "configuration": {
                "type": "object",
                "properties": {"outputDir": {"$ref": "common.json#/configuration/outputDir"}},
                "required": ["outputDir"],
            },
        }
        resolved = resolver.dereference(unresolved)
        compared = {
            "$schema": "http://json-schema.org/draft-07/schema#",
            "type": "object",
            "configuration": {
                "type": "object",
                "properties": {"outputDir": {"type": "string", "description": "Output folder"}},
                "required": ["outputDir"],
            },
        }
        assert resolved == compared
        # reference common.json#/configuration/timeSliceInterval has inside another reference!
        unresolved = {
            "$schema": "http://json-schema.org/draft-07/schema#",
            "type": "object",
            "configuration": {
                "type": "object",
                "properties": {"timeSliceInterval": {"$ref": "common.json#/configuration/timeSliceInterval"}},
            },
        }
        resolved = resolver.dereference(unresolved)
        compared = {
            "$schema": "http://json-schema.org/draft-07/schema#",
            "type": "object",
            "configuration": {
                "type": "object",
                "properties": {
                    "timeSliceInterval": {
                        "description": "Interval for time slicing",
                        "preferredType": "float",
                        "anyOf": [
                            {"type": "string", "pattern": "^$|^[0-9]*.[0-9]*$"},
                            {"type": "number", "exclusiveMinimum": 0},
                        ],
                    }
                },
            },
        }
        assert resolved == compared
        unresolved = {
            "$schema": "http://json-schema.org/draft-07/schema#",
            "type": "object",
            "properties": {"sample": {"$ref": "common.json#/sample"}},
        }
        resolved = resolver.dereference(unresolved)
        compared = {
            "$schema": "http://json-schema.org/draft-07/schema#",
            "type": "object",
            "properties": {
                "sample": {
                    "type": "object",
                    "properties": {
                        "runNumber": {
                            "description": "The run number(s) for the sample",
                            "anyOf": [
                                {
                                    "type": "string",
                                    "minLength": 1,
                                    "pattern": "^[1-9][0-9]*$",
                                },
                                {"type": "integer", "minimum": 1},
                                {
                                    "type": "array",
                                    "items": {"type": ["integer", "string"]},
                                },
                            ],
                            "preferredType": "int",
                        },
                        "transmission": {
                            "type": "object",
                            "properties": {
                                "runNumber": {
                                    "description": "The run number(s) for the transmission sample.",
                                    "anyOf": [
                                        {
                                            "type": "string",
                                            "minLength": 1,
                                            "pattern": "^[1-9][0-9]*$",
                                        },
                                        {"type": "integer", "minimum": 1},
                                        {
                                            "type": "array",
                                            "items": {"type": ["integer", "string"]},
                                        },
                                    ],
                                    "preferredType": "int",
                                },
                                "value": {
                                    "anyOf": [
                                        {"type": "string", "pattern": "^$|^0?.[0-9]*$"},
                                        {
                                            "type": "number",
                                            "exclusiveMinimum": 0,
                                            "maximum": 1,
                                        },
                                    ],
                                    "preferredType": "float",
                                },
                            },
                            "maxProperties": 2,
                            "required": ["runNumber", "value"],
                            "description": "The transmission for the sample",
                        },
                    },
                    "maxProperties": 2,
                    "required": ["transmission", "runNumber"],
                }
            },
        }
        assert resolved == compared
        unresolved = {
            "$schema": "http://json-schema.org/draft-07/schema#",
            "type": "object",
            "properties": {
                "sample": {"$ref": "common.json#/sample"},
                "configuration": {
                    "type": "object",
                    "properties": {"timeSliceInterval": {"$ref": "common.json#/configuration/timeSliceInterval"}},
                },
            },
            "maxProperties": 2,
            "required": ["sample", "configuration"],
        }
        resolved = resolver.dereference(unresolved)
        compared = {
            "$schema": "http://json-schema.org/draft-07/schema#",
            "type": "object",
            "properties": {
                "sample": {
                    "type": "object",
                    "properties": {
                        "runNumber": {
                            "description": "The run number(s) for the sample",
                            "anyOf": [
                                {
                                    "type": "string",
                                    "minLength": 1,
                                    "pattern": "^[1-9][0-9]*$",
                                },
                                {"type": "integer", "minimum": 1},
                                {
                                    "type": "array",
                                    "items": {"type": ["integer", "string"]},
                                },
                            ],
                            "preferredType": "int",
                        },
                        "transmission": {
                            "type": "object",
                            "properties": {
                                "runNumber": {
                                    "description": "The run number(s) for the transmission sample.",
                                    "anyOf": [
                                        {
                                            "type": "string",
                                            "minLength": 1,
                                            "pattern": "^[1-9][0-9]*$",
                                        },
                                        {"type": "integer", "minimum": 1},
                                        {
                                            "type": "array",
                                            "items": {"type": ["integer", "string"]},
                                        },
                                    ],
                                    "preferredType": "int",
                                },
                                "value": {
                                    "anyOf": [
                                        {"type": "string", "pattern": "^$|^0?.[0-9]*$"},
                                        {
                                            "type": "number",
                                            "exclusiveMinimum": 0,
                                            "maximum": 1,
                                        },
                                    ],
                                    "preferredType": "float",
                                },
                            },
                            "maxProperties": 2,
                            "required": ["runNumber", "value"],
                            "description": "The transmission for the sample",
                        },
                    },
                    "maxProperties": 2,
                    "required": ["transmission", "runNumber"],
                },
                "configuration": {
                    "type": "object",
                    "properties": {
                        "timeSliceInterval": {
                            "description": "Interval for time slicing",
                            "anyOf": [
                                {"type": "string", "pattern": "^$|^[0-9]*.[0-9]*$"},
                                {"type": "number", "exclusiveMinimum": 0},
                            ],
                            "preferredType": "float",
                        }
                    },
                },
            },
            "maxProperties": 2,
            "required": ["sample", "configuration"],
        }
        assert resolved == compared
        # array of references
        unresolved = {
            "$schema": "http://json-schema.org/draft-07/schema#",
            "type": "object",
            "properties": {
                "example-multi-type": {
                    "oneOf": [
                        {"$ref": "common.json#/configuration/timeSliceInterval"},
                        {"$ref": "common.json#/configuration/outputDir"},
                    ]
                }
            },
        }
        resolved = resolver.dereference(unresolved)
        compared = {
            "$schema": "http://json-schema.org/draft-07/schema#",
            "type": "object",
            "properties": {
                "example-multi-type": {
                    "oneOf": [
                        {
                            "description": "Interval for time slicing",
                            "anyOf": [
                                {"type": "string", "pattern": "^$|^[0-9]*.[0-9]*$"},
                                {"type": "number", "exclusiveMinimum": 0},
                            ],
                            "preferredType": "float",
                        },
                        {"type": "string", "description": "Output folder"},
                    ]
                }
            },
        }
        assert resolved == compared

    # common.json can resolve itself, thus we include it in the list below
    @pytest.mark.parametrize("file_name", ["common"] + instrument_standard_names())
    def test_derefence_schemas(self, file_name):
        r"""Check all schemas can be resolved using common.json"""
        schema_dir = os.path.join(configdir, "schema")
        json_file = os.path.join(schema_dir, f"{file_name}.json")
        with open(json_file, "r") as file_handle:
            to_resolve = json.load(file_handle)
            resolver_common.dereference(to_resolve)


class TestSuggest:
    property_names = [
        "instrumentName",
        "iptsNumber",
        "sample",
        "thickness",
        "transmission",
        "runNumber",
        "runNumber",
    ]

    @pytest.mark.parametrize(
        "query, score",
        [("home", 0), ("Home", 1), ("hom", 1), ("ho_me", 1), ("Ho_E", 3)],
    )
    def test_levenshtein(self, query, score):
        assert Suggest.levenshtein(query, "home") == score

    def test_init(self):
        entries = Suggest(self.property_names)
        assert len(entries) == len(self.property_names) - 1

    @pytest.mark.parametrize("query, top_match", [("runumber", "runNumber"), ("Sample", "sample")])
    def test_top_match(self, query, top_match):
        entries = Suggest(self.property_names)
        assert entries.top_match(query) == top_match


@pytest.fixture(scope="module")
def default_json(redparams_data):
    resolver = ReferenceResolver(redparams_data["schema_common_file"])
    schema_resolved = resolver.dereference(redparams_data["schema_instrument"])
    return DefaultJson(schema_resolved)


class TestDefaultJson:
    def test_trim_schema(self, default_json):
        compared = {
            "instrumentName": "BIOSANS",
            "iptsNumber": 42,
            "sample": {
                "runNumber": None,
                "transmission": {"runNumber": None, "value": None},
            },
            "configuration": {"outputDir": None, "timeSliceInterval": None},
        }
        assert compared == default_json._json

    def test_str(self, default_json):
        compared = r"""#
# property-name (default value)
#
instrumentName = BIOSANS
iptsNumber = 42
sample:
    runNumber
    transmission:
        runNumber
        value
configuration:
    outputDir
    timeSliceInterval
"""
        assert compared == str(default_json)

    def test_dumps(self, default_json):
        compared = r"""{
  "instrumentName": "BIOSANS",
  "iptsNumber": 42,
  "sample": {
    "runNumber": null,
    "transmission": {
      "runNumber": null,
      "value": null
    }
  },
  "configuration": {
    "outputDir": null,
    "timeSliceInterval": null
  }
}"""
        assert default_json.dumps() == compared

    def test_dump(self, default_json):
        # open temporary file in 'read and write' mode
        _, json_file_path = tempfile.mkstemp(suffix=".json", dir="/tmp")
        default_json.dump(json_file_path)
        d = json.load(open(json_file_path, "r"))
        os.remove(json_file_path)
        compared = {
            "instrumentName": "BIOSANS",
            "iptsNumber": 42,
            "sample": {
                "runNumber": None,
                "transmission": {"runNumber": None, "value": None},
            },
            "configuration": {"outputDir": None, "timeSliceInterval": None},
        }
        assert d == compared

    def test_to_rest(self, default_json):
        d = default_json.to_rest()
        compared = r"""BIOSANS
=======

.. code-block:: python

   {
     "instrumentName": "BIOSANS",
     "iptsNumber": 42,
     "sample": {
       "runNumber": None,
       "transmission": {
         "runNumber": None,
         "value": None
       }
     },
     "configuration": {
       "outputDir": None,
       "timeSliceInterval": None
     }
   }

"""
        assert d == compared

    def test_property_names(self, default_json):
        compared = {
            "instrumentName",
            "iptsNumber",
            "sample",
            "configuration",
            "runNumber",
            "transmission",
            "runNumber",
            "value",
            "outputDir",
            "timeSliceInterval",
        }
        assert compared == default_json.property_names


class TestReductionParameters:
    def test_init(self, redparams_data):
        ReductionParameters(redparams_data["reduction_parameters"], redparams_data["schema_instrument"])

    def test_permissible(self, redparams_data):
        reduction_parameters_new = deepcopy(redparams_data["reduction_parameters"])
        reduction_parameters_new["configuration"]["entry_not_in_the_schema"] = None
        with pytest.raises(KeyError) as error_info:
            ReductionParameters(reduction_parameters_new, redparams_data["schema_instrument"])
        assert "entry_not_in_the_schema" in str(error_info.value)
        ReductionParameters(reduction_parameters_new, redparams_data["schema_instrument"], permissible=True)


class TestReductionParametersGPSANS:
    parameters_common = {
        "instrumentName": "GPSANS",
        "iptsNumber": 21981,
        "sample": {"runNumber": 9165, "thickness": 1.0},
        "outputFileName": "test_validator_datasource",
        "configuration": {"outputDir": "/tmp", "QbinType": "linear", "numQBins": 100},
    }
    parameters_all = reduction_parameters(parameters_common, validate=False)

    @pytest.mark.mount_eqsans
    @pytest.mark.parametrize(
        "validator_name, parameter_changes",
        [
            ("dataSource", {"sample": {"runNumber": 666999666}}),
            ("evaluateCondition", {"configuration": {"numQBins": None}}),
        ],
    )
    def test_validators(self, validator_name, parameter_changes, reference_dir, has_sns_mount):
        if not has_sns_mount:
            pytest.skip("No SNS mount")

        parameter_changes["dataDirectories"] = str(Path(reference_dir.gpsans))
        with pytest.raises(ReductionParameterError) as error_info:
            update_reduction_parameters(self.parameters_all, parameter_changes)
            assert validator_name in str(error_info.schema.keys())

    @pytest.mark.datarepo
    def test_scale_components(self, datarepo_dir):
        r"""Test the validation of the scaleComponents parameter"""

        parameters = deepcopy(self.parameters_all)
        # default value
        # (we need to tell pytest where to find file CG2_9165.nxs.h5)
        with amend_config(data_dir=datarepo_dir.gpsans):
            validate_reduction_parameters(parameters)
        assert parameters["configuration"]["scaleComponents"] == {
            "detector1": None,
            "midrange_detector": None,
            "wing_detector": None,
        }

        # no error cases
        parameters["configuration"]["scaleComponents"] = {"detector1": [1.0, 1.0, 1.0]}
        with amend_config(data_dir=datarepo_dir.gpsans):
            validate_reduction_parameters(parameters)

        # non-existent detector name
        with pytest.raises(KeyError) as exc_info:
            parameters["configuration"]["scaleComponents"] = {"impossible_detector_name": [1.0, 1.0, 1.0]}
            with amend_config(data_dir=datarepo_dir.gpsans):
                validate_reduction_parameters(parameters)
        assert "Parameter impossible_detector_name not found in the schema" in str(exc_info.value)

        # detectors of BIOSANS, not appropriate for GPSANS
        parameters["configuration"]["scaleComponents"] = {"midrange_detector": "", "wing_detector": ""}
        # no scaling factors don't raise an error
        with amend_config(data_dir=datarepo_dir.gpsans):
            validate_reduction_parameters(parameters)
        # attempting to scale a detector of BIOSANS when the instrument is GPSANS should raise an error
        for component in ["midrange_detector", "wing_detector"]:
            with pytest.raises(ReductionParameterError) as exc_info:
                parameters["configuration"]["scaleComponents"] = {component: [1.0, 1.0, 1.0]}
                with amend_config(data_dir=datarepo_dir.gpsans):
                    validate_reduction_parameters(parameters)
            assert f"{component} is not scalable" in str(exc_info.value)

        # wrong type for scaling value
        with pytest.raises(ReductionParameterError):
            parameters["configuration"]["scaleComponents"] = {"detector1": "one point three"}
            with amend_config(data_dir=datarepo_dir.gpsans):
                validate_reduction_parameters(parameters)

        # wrong scaling value
        with pytest.raises(ReductionParameterError):
            parameters["configuration"]["scaleComponents"] = {"detector1": [-1.0, 1.0, 1.0]}
            with amend_config(data_dir=datarepo_dir.gpsans):
                validate_reduction_parameters(parameters)

        # wrong number of scaling values
        with pytest.raises(ReductionParameterError):
            parameters["configuration"]["scaleComponents"] = {"detector1": [-1.0, 1.0]}
            with amend_config(data_dir=datarepo_dir.gpsans):
                validate_reduction_parameters(parameters)

    @pytest.mark.datarepo
    def test_permissible(self, datarepo_dir):
        parameters_new = update_reduction_parameters(
            self.parameters_common,
            {"dataDirectories": datarepo_dir.gpsans, "configuration": {"entry_not_in_the_schema": None}},
            validate=False,
            permissible=True,
        )
        reduction_parameters(parameters_new, permissible=True)


class TestReductionParametersBIOSANS:
    parameters_common = {
        "instrumentName": "BIOSANS",
        "iptsNumber": "23782",
        "sample": {"runNumber": "960", "transmission": {"runNumber": ""}},
        "outputFileName": "test_validator_biosans",
        "configuration": {
            "useDefaultMask": False,
            "outputDir": "/tmp",
            "QbinType": "linear",
            "numMainQBins": 100,
            "numWingQBins": 100,
            "numMidrangeQBins": 100,
        },
    }
    parameters_all = reduction_parameters(parameters_common, validate=False)

    @pytest.mark.mount_eqsans
    def test_validators_midrange_parameters_required(self, reference_dir, has_sns_mount):
        if not has_sns_mount:
            pytest.skip("No SNS mount")

        parameters = deepcopy(self.parameters_all)
        parameters["dataDirectories"] = str(Path(reference_dir.biosans))
        # remove all parameters related to the midrange detector
        config_no_midrange = {k: v for k, v in parameters["configuration"].items() if "Midrange" not in k}
        parameters["configuration"] = config_no_midrange
        with pytest.raises(ReductionParameterError) as error_info:
            validate_reduction_parameters(parameters)
        assert "'darkMidrangeFileName' is a required property" in str(error_info.value)

    @pytest.mark.mount_eqsans
    def test_validators_midrange_qmin_qmax(self, reference_dir, has_sns_mount):
        if not has_sns_mount:
            pytest.skip("No SNS mount")

        parameters = deepcopy(self.parameters_all)
        parameters["dataDirectories"] = str(Path(reference_dir.biosans))
        parameters["configuration"]["QminMidrange"] = 0.07
        parameters["configuration"]["QmaxMidrange"] = 0.05
        with pytest.raises(ReductionParameterError) as error_info:
            validate_reduction_parameters(parameters)
        assert "0.07 is not smaller than #configuration/QmaxMidrange" in str(error_info.value)

    @pytest.mark.mount_eqsans
    @pytest.mark.parametrize(
        "qmin_name, qmax_name",
        [
            ("overlapStitchQmin", "overlapStitchQmax"),
            ("wedge1overlapStitchQmin", "wedge1overlapStitchQmax"),
            ("wedge2overlapStitchQmin", "wedge2overlapStitchQmax"),
        ],
    )
    @pytest.mark.parametrize(
        "throws_error, qmin_value, qmax_value",
        [
            (True, None, [0.015]),  # different length
            (True, [0.015], None),  # different length
            (True, [0.01, 0.02], [0.015]),  # different length
            (True, [0.02], [0.01]),  # min > max
            (True, [0.01, 0.02, 0.03], [0.015, 0.025, 0.035]),  # lists too long
            # valid inputs:
            (False, None, None),
            (False, None, None),
            (False, [], []),
            (False, [], []),
            (False, 0.01, 0.015),
            (False, [0.01, 0.02], [0.015, 0.025]),
        ],
    )
    def test_overlap_stitch(
        self, reference_dir, qmin_name, qmax_name, throws_error, qmin_value, qmax_value, has_sns_mount
    ):
        if not has_sns_mount:
            pytest.skip("No SNS mount")

        parameters = deepcopy(self.parameters_all)
        parameters["dataDirectories"] = str(Path(reference_dir.biosans))
        parameters["configuration"][qmin_name] = qmin_value
        parameters["configuration"][qmax_name] = qmax_value
        if throws_error:
            with pytest.raises(ReductionParameterError) as error_info:
                validate_reduction_parameters(parameters)
            assert qmin_name in str(error_info.value)
        else:
            validate_reduction_parameters(parameters)

    @pytest.mark.datarepo
    def test_scale_components(self, datarepo_dir):
        r"""Test the validation of the scaleComponents parameter"""

        parameters = deepcopy(self.parameters_all)

        # check default value
        # (we need to tell pytest where to find file CG3_960.nxs.h5)
        with amend_config(data_dir=path_join(datarepo_dir.biosans, "pixel_calibration", "test_loader_algorithm")):
            validate_reduction_parameters(parameters)
        assert parameters["configuration"]["scaleComponents"] == {
            "detector1": None,
            "midrange_detector": None,
            "wing_detector": None,
        }

        # no error cases
        parameters["configuration"]["scaleComponents"] = {
            "detector1": [1.0, 1.0, 1.0],
            "wing_detector": "",
            "midrange_detector": None,
        }
        with amend_config(data_dir=path_join(datarepo_dir.biosans, "pixel_calibration", "test_loader_algorithm")):
            validate_reduction_parameters(parameters)

        # wrong detector name
        with pytest.raises(KeyError) as exc_info:
            parameters["configuration"]["scaleComponents"] = {"impossible_detector_name": [1.0, 1.0, 1.0]}
            with amend_config(data_dir=path_join(datarepo_dir.biosans, "pixel_calibration", "test_loader_algorithm")):
                validate_reduction_parameters(parameters)
        assert "Parameter impossible_detector_name not found in the schema" in str(exc_info.value)

        # wrong type for scaling value
        with pytest.raises(ReductionParameterError):
            parameters["configuration"]["scaleComponents"] = {"detector1": "one point three"}
            with amend_config(data_dir=path_join(datarepo_dir.biosans, "pixel_calibration", "test_loader_algorithm")):
                validate_reduction_parameters(parameters)

        # wrong scaling value
        with pytest.raises(ReductionParameterError):
            parameters["configuration"]["scaleComponents"] = {"detector1": [-1.0, 1.0, 1.0]}
            with amend_config(data_dir=path_join(datarepo_dir.biosans, "pixel_calibration", "test_loader_algorithm")):
                validate_reduction_parameters(parameters)

        # wrong number of scaling values
        with pytest.raises(ReductionParameterError):
            parameters["configuration"]["scaleComponents"] = {"detector1": [-1.0, 1.0]}
            with amend_config(data_dir=path_join(datarepo_dir.biosans, "pixel_calibration", "test_loader_algorithm")):
                validate_reduction_parameters(parameters)

    @pytest.mark.datarepo
    def test_permissible(self, datarepo_dir):
        parameters_new = update_reduction_parameters(
            self.parameters_common,
            {
                "dataDirectories": os.path.join(datarepo_dir.biosans, "pixel_calibration", "test_loader_algorithm"),
                "configuration": {"entry_not_in_the_schema": None},
            },
            validate=False,
            permissible=True,
        )
        reduction_parameters(parameters_new, permissible=True)


class TestReductionParametersEQSANS:
    parameters_common = {
        "instrumentName": "EQSANS",
        "iptsNumber": 20196,
        "sample": {"runNumber": 89157, "thickness": 1.0},
        "outputFileName": "test_validator_datasource",
        "configuration": {"outputDir": "/tmp", "QbinType": "linear", "numQBins": 100},
    }
    parameters_all = reduction_parameters(parameters_common, validate=False)

    @pytest.mark.datarepo
    def test_incohfit_parameters(self, datarepo_dir):
        parameters = deepcopy(self.parameters_all)
        del parameters["configuration"]["blockedBeamRunNumber"]
        with amend_config(data_dir=datarepo_dir.eqsans):
            validate_reduction_parameters(parameters)

        # assert incohfit_qmin/qmax/factor are null by default
        assert parameters["configuration"]["incohfit_qmin"] is None
        assert parameters["configuration"]["incohfit_qmax"] is None
        assert parameters["configuration"]["incohfit_factor"] is None

        # assert incohfit_qmin/qmax/factor can be set to a float
        parameters["configuration"]["incohfit_qmin"] = 0.01
        parameters["configuration"]["incohfit_qmax"] = 0.1
        parameters["configuration"]["incohfit_factor"] = 10
        with amend_config(data_dir=datarepo_dir.eqsans):
            validate_reduction_parameters(parameters)

        # assert incohfit_qmin/qmax/factor can be set to a list of floats
        parameters["configuration"]["incohfit_qmin"] = [0.01, 0.02]
        parameters["configuration"]["incohfit_qmax"] = [0.1, 0.2]
        parameters["configuration"]["incohfit_factor"] = [10, 20]
        with amend_config(data_dir=datarepo_dir.eqsans):
            validate_reduction_parameters(parameters)

        # assert incohfit_qmin/qmax/factor must be lists length 1 or 2
        with pytest.raises(ReductionParameterError):
            parameters["configuration"]["incohfit_qmin"] = [0.01, 0.02, 0.03]
            parameters["configuration"]["incohfit_qmax"] = []
            parameters["configuration"]["incohfit_factor"] = [10, 20, 30, 40]
            with amend_config(data_dir=datarepo_dir.eqsans):
                validate_reduction_parameters(parameters)

    @pytest.mark.datarepo
    def test_scale_components(self, datarepo_dir):
        r"""Test the validation of the scaleComponents parameter"""

        parameters = deepcopy(self.parameters_all)
        del parameters["configuration"]["blockedBeamRunNumber"]
        # default value
        # (we need to tell pytest where to find file EQSANS_89157.nxs.h5)
        with amend_config(data_dir=datarepo_dir.eqsans):
            validate_reduction_parameters(parameters)
        assert parameters["configuration"]["scaleComponents"] == {
            "detector1": None,
            "midrange_detector": None,
            "wing_detector": None,
        }

        # no error cases
        parameters["configuration"]["scaleComponents"] = {"detector1": [1.0, 1.0, 1.0]}
        with amend_config(data_dir=datarepo_dir.eqsans):
            validate_reduction_parameters(parameters)

        # non-existent detector name
        with pytest.raises(KeyError) as exc_info:
            parameters["configuration"]["scaleComponents"] = {"impossible_detector_name": [1.0, 1.0, 1.0]}
            with amend_config(data_dir=datarepo_dir.eqsans):
                validate_reduction_parameters(parameters)
        assert "Parameter impossible_detector_name not found in the schema" in str(exc_info.value)

        # detectors of BIOSANS, not appropriate for EQSANS
        parameters["configuration"]["scaleComponents"] = {"midrange_detector": "", "wing_detector": ""}
        # no scaling factors don't raise an error
        with amend_config(data_dir=datarepo_dir.eqsans):
            validate_reduction_parameters(parameters)
        # attempting to scale a detector of BIOSANS when the instrument is EQSANS should raise an error
        for component in ["midrange_detector", "wing_detector"]:
            with pytest.raises(ReductionParameterError) as exc_info:
                parameters["configuration"]["scaleComponents"] = {component: [1.0, 1.0, 1.0]}
                with amend_config(data_dir=datarepo_dir.eqsans):
                    validate_reduction_parameters(parameters)
            assert f"{component} is not scalable" in str(exc_info.value)

        # wrong type for scaling value
        with pytest.raises(ReductionParameterError):
            parameters["configuration"]["scaleComponents"] = {"detector1": "one point three"}
            with amend_config(data_dir=datarepo_dir.eqsans):
                validate_reduction_parameters(parameters)

        # wrong scaling value
        with pytest.raises(ReductionParameterError):
            parameters["configuration"]["scaleComponents"] = {"detector1": [-1.0, 1.0, 1.0]}
            with amend_config(data_dir=datarepo_dir.eqsans):
                validate_reduction_parameters(parameters)

        # wrong number of scaling values
        with pytest.raises(ReductionParameterError):
            parameters["configuration"]["scaleComponents"] = {"detector1": [-1.0, 1.0]}
            with amend_config(data_dir=datarepo_dir.eqsans):
                validate_reduction_parameters(parameters)

    @pytest.mark.parametrize(
        "timeSliceInterval, timeSlicePeriod, throws_error",
        [
            (10, 100, False),  # exact integer multiple
            (0.001, 1.00000000001, False),  # nearly exact integer ratio
            (-45, -100, True),  # negative values
            (1, 100.1, False),  # not an integer multiple
        ],
        ids=["exact integer multiple", "nearly exact integer ratio", "negative values", "not an integer multiple"],
    )
    def test_timeslice_parameters(self, datarepo_dir, timeSliceInterval, timeSlicePeriod, throws_error):
        parameters = deepcopy(self.parameters_all)
        del parameters["configuration"]["blockedBeamRunNumber"]
        # default values
        with amend_config(data_dir=datarepo_dir.eqsans):
            validate_reduction_parameters(parameters)
        assert parameters["configuration"]["useTimeSlice"] is False
        assert parameters["configuration"]["timeSliceInterval"] == 300
        assert parameters["configuration"]["timeSlicePeriod"] is None

        # set timeSliceInterval and timeSlicePeriod
        parameters["configuration"]["useTimeSlice"] = True
        parameters["configuration"]["timeSliceInterval"] = timeSliceInterval
        parameters["configuration"]["timeSlicePeriod"] = timeSlicePeriod

        with amend_config(data_dir=datarepo_dir.eqsans):
            if throws_error:
                with pytest.raises(ReductionParameterError):
                    validate_reduction_parameters(parameters)
            else:
                validate_reduction_parameters(parameters)

    @pytest.mark.datarepo
    def test_permissible(self, datarepo_dir):
        parameters_new = update_reduction_parameters(
            self.parameters_common,
            {"dataDirectories": datarepo_dir.eqsans, "configuration": {"entry_not_in_the_schema": None}},
            validate=False,
            permissible=True,
        )
        del parameters_new["configuration"]["blockedBeamRunNumber"]
        reduction_parameters(parameters_new, permissible=True)

    @pytest.mark.datarepo
    def test_loadOptions_additionalProperties_true(self, datarepo_dir):
        parameters = deepcopy(self.parameters_all)
        del parameters["configuration"]["blockedBeamRunNumber"]
        with amend_config(data_dir=datarepo_dir.eqsans):
            parameters["sample"]["loadOptions"]["additionalProperties"] = True
            validate_reduction_parameters(parameters)


def test_generate_json_files(tmpdir, cleanfile):
    directory = tmpdir.mkdir("generate_json_files")
    cleanfile(directory)
    generate_json_files(directory)


class TestReductionParameterError:
    @pytest.fixture
    def mock_runNumber_exists(self, monkeypatch):
        """Trick the validation to think that the runNumber exists"""
        fake_path = "fake/path"

        def mock_exists(path):
            return True

        def mock_abspath(path):
            return fake_path

        monkeypatch.setattr(os.path, "exists", mock_exists)
        monkeypatch.setattr(os.path, "abspath", mock_abspath)

    parameters_EQSANS = {
        "instrumentName": "EQSANS",
        "iptsNumber": 20196,
        "dataDirectories": None,
        "sample": {
            "runNumber": 89157,
            "thickness": 1.0,
            "transmission": {"runNumber": 59157, "value": 1.0},
            "loadOptions": {"additionalProperties": True},
        },
        "background": {
            "runNumber": 89158,
            "transmission": {"runNumber": 59158, "value": 1.0},
        },
        "emptyTransmission": {"runNumber": 89159, "value": 1.0},
        "beamCenter": {"runNumber": 89160},
        "outputFileName": "tmp",
        "configuration": {
            "outputDir": "/tmp",
            "instrumentConfigurationDir": "/SNS/EQSANS/shared/instrument_configuration",
            "useTimeSlice": False,
            "timeSliceInterval": 300,
            "incohfit_qmin": None,
            "incohfit_qmax": None,
            "QbinType": "log",
            "numQBins": 100,
            "darkFileName": None,
        },
    }

    def test_reduction_parameters_nominal(self, mock_runNumber_exists):
        reduction_parameters(self.parameters_EQSANS, validate=True)

    def test_reduction_parameters_missing_runNumber(self, mock_runNumber_exists):
        parameters = deepcopy(self.parameters_EQSANS)
        parameters["sample"]["runNumber"] = None
        with pytest.raises(ReductionParameterError) as error_info:
            reduction_parameters(parameters, validate=True)
        assert "runNumber" in str(error_info.value)

    def test_runNumber_wrong_type(self, mock_runNumber_exists):
        parameters = deepcopy(self.parameters_EQSANS)
        parameters["iptsNumber"] = 8.75
        with pytest.raises(ReductionParameterError) as error_info:
            reduction_parameters(parameters, validate=True)
            assert "Description: " in str(error_info.message)
        assert "iptsNumber" in str(error_info.value)

    def test_max_less_than_min(self, mock_runNumber_exists):
        parameters = deepcopy(self.parameters_EQSANS)
        parameters["configuration"]["incohfit_qmin"] = 2000
        parameters["configuration"]["incohfit_qmax"] = 500
        schema = load_schema("EQSANS")
        with pytest.raises(ReductionParameterError) as error_info:
            reduction_parameters(parameters, validate=True)
            assert schema["properties"]["configuration"]["properties"]["incohfit_qmax"]["description"] in str(
                error_info.message
            )
        assert "incohfit_qmax" in str(error_info.value)

    def test_common_description(self, mock_runNumber_exists):
        parameters = deepcopy(self.parameters_EQSANS)
        parameters["configuration"]["timeSliceInterval"] = -45
        schema = load_schema("EQSANS")
        with pytest.raises(ReductionParameterError) as error_info:
            reduction_parameters(parameters, validate=True)
            assert schema["properties"]["configuration"]["properties"]["timeSliceInterval"]["description"] in str(
                error_info.message
            )
            assert schema["properties"]["configuration"]["properties"]["timeSliceInterval"]["default"] in str(
                error_info.message
            )
        assert "timeSliceInterval" in str(error_info.value)

    def test_missing_required(self, mock_runNumber_exists):
        parameters = deepcopy(self.parameters_EQSANS)
        del parameters["configuration"]["numQBins"]
        with pytest.raises(ReductionParameterError) as error_info:
            reduction_parameters(parameters, validate=True)
        assert "numQBins" in str(error_info.value)

    def test_eval_condition(self, mock_runNumber_exists):
        parameters = deepcopy(self.parameters_EQSANS)
        parameters["configuration"]["numQBins"] = None
        with pytest.raises(ReductionParameterError) as error_info:
            reduction_parameters(parameters, validate=True)
        assert "numQBins" in str(error_info.value)

    def test_print_error_property(self, mock_runNumber_exists):
        parameters = deepcopy(self.parameters_EQSANS)
        parameters["configuration"]["numQBins"] = None
        with pytest.raises(ReductionParameterError) as error_info:
            reduction_parameters(parameters, validate=True)
            assert "numQBins" in error_info.error_user_friendly


if __name__ == "__main__":
    pytest.main([__file__])
