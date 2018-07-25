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

# One test
pytest -s -v  tests/test_eqsansload.py::test_get_config_file

# One file
pytest -s -v  tests/test_eqsansload.py

# All
pytest -s -v

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
(issue #2) Framework/WorkflowAlgorithms/src/EQSANSLoad.cpp
Framework/WorkflowAlgorithms/src/EQSANSMonitorTOF.cpp
Framework/WorkflowAlgorithms/src/EQSANSQ2D.cpp
Framework/WorkflowAlgorithms/src/HFIRDarkCurrentSubtraction.cpp
Framework/WorkflowAlgorithms/src/HFIRLoad.cpp
Framework/WorkflowAlgorithms/src/HFIRSANSNormalise.cpp
(issue #1) Framework/WorkflowAlgorithms/src/SANSBeamFinder.cpp
Framework/WorkflowAlgorithms/src/SANSBeamFluxCorrection.cpp
Framework/WorkflowAlgorithms/src/SANSSensitivityCorrection.cpp
Framework/WorkflowAlgorithms/src/SANSSolidAngleCorrection.cpp
Framework/WorkflowAlgorithms/src/SetupEQSANSReduction.cpp
Framework/WorkflowAlgorithms/src/SetupHFIRReduction.cpp
```
