import json
import jsonschema
import sys


schema = {
    "$id": "https://shaman.ornl.gov/drtsans.schema.json",
    "$schema": "http://json-schema.org/draft-07/schema#",
    "title": "drtsans input",
    "type": "object",
    "properties": {
        "runNumber": {
            "type": "string",
            "minLength": 2,
            "description": "the run number"
        },
        "outputFilename": {
            "type": "string",
            "description": "the output filename"
        },
        "instrumentName": {
            "type": "string",
            "enum": ["EQSANS", "CG2", "CG3"],
            "description": "the instrument name"
        },
        "iptsNumber": {
            "type": "string",
            "description": "the IPTS of the experiment"
        },
        "thickness": {
            "type": "string",
            "description": "sample thickness in cm"
        },
        "transmission": {
            "type": "object",
            "properties": {
                "runNumber": {
                    "type": "string",
                    "description": "the run number for transmission"
                }
            },
            "required": ["runNumber"],
            "additionalProperties": False
        },
        "background": {
            "type": "object",
            "properties": {
                "runNumber": {
                    "type": "string",
                    "description": "the run number for the background"
                },
                "transmission": {
                    "type": "object",
                    "properties": {
                        "runNumber": {
                            "type": "string",
                            "description": "the run number for transmission"
                        }
                    },
                    "required": ["runNumber"],
                    "additionalProperties": False
                }
            },
            "required": ["runNumber", "transmission"],
            "additionalProperties": False
        },
        "empty": {
            "type": "object",
            "properties": {
                "runNumber": {
                    "type": "string",
                    "description": "the run number for empty"
                }
            },
            "required": ["runNumber"],
            "additionalProperties": False
        },
        "configuration": {
            "type": "object",
            "properties": {
                "outputDir": {
                    "type": "string",
                    "description": "destination folder",
                    "minLength": 3
                },
                "sampleApertureSize": {
                    "type": "string",
                    "description": "sample aperture size (will override the one from the file if set)"
                },
                "maskFileName": {
                    "type": "string",
                    "description": "filename for a mask"
                },
                "useMaskFileName": {
                    "type": "boolean",
                    "description": "can choose not to use the mask file, even if given"
                },
                "useDefaultMask": {
                    "type": "boolean",
                    "description": "flag to use default mask"
                },
                "beamFluxFileName": {
                    "type": "string",
                    "description": "filename for beam flux"
                },
                "useBeamFluxFileName": {
                    "type": "boolean",
                    "description": "can choose not to use the beam flux file, even if given"
                },
                "darkFileName": {
                    "type": "string",
                    "description": "filename for the dark current measurement"
                },
                "useDarkFileName": {
                    "type": "boolean",
                    "description": "can choose not to use the dark current file, even if given"
                },
                "sensitivityFileName": {
                    "type": "string",
                    "description": "filename for processed sensitivity"
                },
                "useSensitivityFileName": {
                    "type": "boolean",
                    "description": "can choose not to use the sensitivity file, even if given"
                },
                "absoluteScale": {
                    "type": "string",
                    "description": "absolute scale factor (should be able to convert to a floating point number)"
                },
                "normalization": {
                    "type": "string",
                    "enum": ["Total charge", "Monitor", "Time", "None"],
                    "description": "selector for normalization type"
                },
                "sampleOffset": {
                    "type": "string",
                    "description": "sample offset (should be able to convert to a floating point number)"
                },
                "useSampleOffset": {
                    "type": "boolean",
                    "description": "override the sample offset value in the file with the one in the configuration"
                },
                "detectorOffset": {
                    "type": "string",
                    "description": "detector offset (should be able to convert to a floating point number)"
                },
                "useDetectorOffset": {
                    "type": "boolean",
                    "description": "override the detector offset value in the file with the one in the configuration"
                },
                "useSolidAngleCorrection": {
                    "type": "boolean",
                    "description": "flag to apply solid angle correction"
                },
                "useDetectorTubeType": {
                    "type": "boolean",
                    "description": "flag to use the tube type when calculating solid angle correction"
                },
                "useFlightPathCorrection": {
                    "type": "boolean",
                    "description": "flag to ???"
                },
                "useThetaDepTransCorrection": {
                    "type": "boolean",
                    "description": "flag to ???"
                },
                "nPixelsRadiusForTransmission": {
                    "type": "string",
                    "description": "??? (should be able to convert to an integer)"
                },
                "numQxQyBins": {
                    "type": "integer",
                    "description": "number of bins in qx and qy"
                },
                "QbinType": {
                    "type": "string",
                    "enum": ["linear", "log"],
                    "description": "1D binning type, either 'linear' or 'log'",
                },
                "numQBins": {
                    "type": "integer",
                    "description": "number of bins in |q|"
                },
                "useTOFcuts": {
                    "type": "boolean",
                    "description": "flag to use TOF limits in the data"
                },
                "TOFmin": {
                    "type": "string",
                    "description": "minimum TOF (should be able to convert to a floating point number)"
                },
                "TOFmax": {
                    "type": "string",
                    "description": "maximum TOF (should be able to convert to a floating point number)"
                },
                "useMaskBackTubes": {
                    "type": "boolean",
                    "description": "flag to mask the tubes in the back"
                },
                "wavelenStepType": {
                    "type": "string",
                    "enum": ["constant lambda", "constant dlambda/lambda"],
                    "description": "wether the lambda binning is linear or logarithmic"
                },
                "wavelenStep": {
                    "type": "number",
                    "description": "wavelength step"
                }
            },
            "minProperties": 31,
            "additionalProperties": False
        }
    },
    "minProperties": 9,
    "additionalProperties": False
}


def validate_file(filename):
    r"""
    Function to validate a json file acording to the schema above

    Parameters
    ----------

    filename: str
        the name of the jkson file to be validated
    """
    with open(filename, 'r') as jsf:
        json_parameters = json.load(jsf)
    validator = jsonschema.Draft4Validator(schema)
    if validator.is_valid(json_parameters):
        print("Valid drtsans json")
    else:
        error_message = ''
        for error in sorted(validator.iter_errors(json_parameters), key=str):
            error_message += '\\'.join(error.absolute_path) + ':\t' + error.message + '\n'
        print(error_message)


if __name__ == "__main__":
    """
    usage: python3 drtsans_schema.py reduction.json
    """
    filename = sys.argv[1]
    validate_file(filename)
