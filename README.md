# New SANS Command interface

Here we are going to split workflow algorithms into functions.

# Enviroment file

```sh
cp env.base .env

# Edit the .env file according to your paths

```

# Use a virtual environment

```sh
# Create the environment with python 2 (for Mantid it has to be python 2....)
virtualenv -p $(which python2) venv

# Activate it
source venv/bin/activate

# Install the requirements
pip install -r requirements.txt

```

# Run the tests


```sh

# Run one single test
pytest -s -v  tests/legacy/test_eqsansload.py::test_get_config_file

# Run all tests within the same file
pytest -s -v  tests/legacy/test_eqsansload.py

# Run all tests in the folder tests
pytest -s -v

```

# How to use the PropertyManager

This might be useful to debug old algorithms

```python
from __future__ import print_function

from mantid.kernel import PropertyManagerDataService, PropertyManager

# Create it
pm = PropertyManager()
pmds = PropertyManagerDataService.add("pm_name", pm)

# Add properties
pm.declareProperty("p1", "v1")
pm.declareProperty("p2", "v2")

# print all the properties
for k, v in zip(pm.keys(), pm.values()):
    print("{} -> ".format(k), end="")
    try:
        print(v.value)
    except:
        print(v)
```

# List to split

**Python Algorithms**

```sh
Framework/PythonInterface/plugins/algorithms/WorkflowAlgorithms/EQSANSAzimuthalAverage1D.py
Framework/PythonInterface/plugins/algorithms/WorkflowAlgorithms/EQSANSDirectBeamTransmission.py
Framework/PythonInterface/plugins/algorithms/WorkflowAlgorithms/EQSANSNormalise.py
Framework/PythonInterface/plugins/algorithms/WorkflowAlgorithms/HFIRSANSReduction.py
Framework/PythonInterface/plugins/algorithms/WorkflowAlgorithms/NormaliseByThickness.py
Framework/PythonInterface/plugins/algorithms/WorkflowAlgorithms/ReactorSANSResolution.py
Framework/PythonInterface/plugins/algorithms/WorkflowAlgorithms/SANSAbsoluteScale.py
Framework/PythonInterface/plugins/algorithms/WorkflowAlgorithms/SANSAzimuthalAverage1D.py
Framework/PythonInterface/plugins/algorithms/WorkflowAlgorithms/SANSBeamSpreaderTransmission.py
Framework/PythonInterface/plugins/algorithms/WorkflowAlgorithms/SANSDirectBeamTransmission.py

# This one probably not!!!
Framework/PythonInterface/plugins/algorithms/WorkflowAlgorithms/SANSMask.py

Framework/PythonInterface/plugins/algorithms/WorkflowAlgorithms/SANSReduction.py
Framework/PythonInterface/plugins/algorithms/WorkflowAlgorithms/TransmissionUtils.py
```


**C++ Algorihtms**

```sh
Framework/WorkflowAlgorithms/src/EQSANSDarkCurrentSubtraction.cpp

# issue #2
Framework/WorkflowAlgorithms/src/EQSANSLoad.cpp

Framework/WorkflowAlgorithms/src/EQSANSMonitorTOF.cpp
Framework/WorkflowAlgorithms/src/EQSANSQ2D.cpp
Framework/WorkflowAlgorithms/src/HFIRDarkCurrentSubtraction.cpp
Framework/WorkflowAlgorithms/src/HFIRLoad.cpp
Framework/WorkflowAlgorithms/src/HFIRSANSNormalise.cpp

# issue #1
Framework/WorkflowAlgorithms/src/SANSBeamFinder.cpp

Framework/WorkflowAlgorithms/src/SANSBeamFluxCorrection.cpp
Framework/WorkflowAlgorithms/src/SANSSensitivityCorrection.cpp
Framework/WorkflowAlgorithms/src/SANSSolidAngleCorrection.cpp
Framework/WorkflowAlgorithms/src/SetupEQSANSReduction.cpp
Framework/WorkflowAlgorithms/src/SetupHFIRReduction.cpp
```
