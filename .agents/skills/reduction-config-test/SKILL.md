---
name: reduction-config-test
description: Use when the user asks to create a drtsans integration test from, with, or using a JSON file path, especially prompts like "create a test with JSON file /path/to/config.json"; generates a no-assertion reduction test that loads all files and reduces one configuration.
---

# Reduction Config Test

Use this skill to turn a JSON reduction configuration into a focused integration test that loads
all files and reduces that one configuration. The generated test reproduces behavior only; it must
not assert on outputs or compare to expected data.

## Workflow

1. Read the input JSON file and parse it as JSON.
   - Do not manually transcribe JSON by hand when a parser is available.
   - Embed the parsed JSON in the generated test as a Python dictionary literal.

2. Determine the instrument from the configuration.
   - Prefer `instrumentName`.
   - Treat `EQSANS` or `EQ-SANS` as `tests/integration/drtsans/tof/eqsans/`.
   - Treat `GPSANS` or `CG2` as `tests/integration/drtsans/mono/gpsans/`.
   - Treat `BIOSANS` or `CG3` as `tests/integration/drtsans/mono/biosans/`.
   - If the instrument is missing or unsupported, ask the user which supported instrument path to use.

3. Propose the test file path.
   - Default filename: `test_reduction_config.py`.
   - Default full path: the instrument directory plus `test_reduction_config.py`.
   - Show the user the full path and ask them to confirm or provide a different filename.
   - If the user provides only a filename, keep it in the selected instrument directory.
   - Use a pytest-compatible filename matching this repository's `python_files = ["test*.py"]` configuration.

4. Propose the output directory amendment.
   - Default output directory: `/tmp/reduction_config_test`.
   - Ask the user to confirm this path or provide another output directory.
   - Ensure the generated test sets `reduction_input["configuration"]["outputDir"]` to the confirmed path.
   - If the parsed JSON has no `configuration` dictionary, create one before assigning `outputDir`.

5. Generate the test.
   - Add imports for the selected instrument:
     - EQSANS: `from drtsans.tof.eqsans import load_all_files, reduction_parameters, reduce_single_configuration`
     - GPSANS/CG2: `from drtsans.mono.gpsans import load_all_files, reduction_parameters, reduce_single_configuration`
     - BIOSANS/CG3: `from drtsans.mono.biosans import load_all_files, reduction_parameters, reduce_single_configuration`
   - Include `import pytest`.
   - Mark the test with `@pytest.mark.datarepo`.
   - Name the function `test_reduction_config`.
   - The body must contain, in this order:
     1. `reduction_input = {...}` using the parsed JSON dictionary literal.
     2. The `configuration.outputDir` amendment.
     3. `reduction_input = reduction_parameters(parameters_particular=reduction_input, validate=True)`.
     4. `loaded = load_all_files(reduction_input)`.
     5. `reduce_single_configuration(loaded, reduction_input)`.
   - Do not add assertions.
   - Do not add expected-data checks, gold-file comparisons, or output-file existence checks.
   - Do not add cleanup unless the user explicitly asks for it.

6. Do not run the generated test.
   - Do not run integration tests, unit tests, pytest collection, or import checks for the generated file.
   - If verification is needed, limit it to static inspection of the generated file content.

## Template

Use this structure, selecting the import path for the instrument:

```python
import pytest

from drtsans.tof.eqsans import load_all_files, reduction_parameters, reduce_single_configuration


@pytest.mark.datarepo
def test_reduction_config():
    reduction_input = {
        # parsed JSON dictionary literal
    }
    reduction_input.setdefault("configuration", {})
    reduction_input["configuration"]["outputDir"] = "/tmp/reduction_config_test"
    reduction_input = reduction_parameters(parameters_particular=reduction_input, validate=True)

    loaded = load_all_files(reduction_input)
    reduce_single_configuration(loaded, reduction_input)
```
