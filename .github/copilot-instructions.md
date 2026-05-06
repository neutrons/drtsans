# GitHub Copilot Instructions for drtsans

This file configures how GitHub Copilot should assist with development in the
**drtsans** repository. This project is the ORNL SANS (Small-Angle Neutron
Scattering) instruments reduction package, used by scientists and researchers
for neutron scattering data analysis.

## 🎯 Core Principles

1. **Always assess before acting** - Understand the current state before proposing changes
2. **Always provide itemized plans** - Break work into clear, testable steps
3. **Focus on progress tracking** - Users should always know what's happening and what's next
4. **Test incrementally** - Each step should be testable
5. **Review after major changes** - Use sub-agents for code review when appropriate
6. **Document findings** - Keep track of important discoveries and patterns in
   code comments and documentation

---

## 🚫 What to Avoid

- **Do not refactor code that was not requested** — only change what is needed.
- **Do not add unrequested features** (extra parameters, new abstractions, etc.).
- **Do not use `print()`** — always use the Mantid logger for messages.
- **Do not guess at file paths or workspace names** — read the relevant source files first.
- **Do not skip tests** — every implementation step must be paired with a test.
- **Do not over-engineer** — prefer simple, readable code that matches the existing style.

---

## 🔄 Standard Workflow

### Step 1: Assess
Before responding to any request, examine:
- Current implementation (read relevant files)
- Existing tests (unit and integration tests in `tests/`)
- Related code that might be affected
- Project structure and conventions
- Dependencies on Mantid framework

**Output**: Provide a brief summary of what exists now and what needs to change.

**When the request is ambiguous**: check the codebase first — the answer is often in
existing patterns or related functions. If still unclear, state your assumption explicitly
(e.g., *"I'm assuming `q_max` defaults to `0.5` based on the existing API — let me know
if you intended otherwise"*) and proceed. Ask at most one focused question only if the
ambiguity would cause materially different implementations.

**Example Assessment**:
```
Current state:
- src/drtsans/iq.py has binning functions for I(Q) profiles
- Unit tests exist in tests/unit/
- Integration tests use pytest fixtures from tests/conftest.py
- Code uses Mantid framework for data manipulation
- Ruff is configured for linting

Need to change:
- Add new binning method for annular integration
- Add type hints (already partially present)
- Add comprehensive tests
- Update documentation in docs/
```

### Step 2: Plan
Create an itemized, actionable plan with:
- Clear, numbered steps
- Each step is independently testable
- Dependencies between steps noted
- Expected test coverage outlined

**Format**:
```markdown
## Itemized Plan

### Core Implementation (3 steps)
1. **Add new binning function** - Implement in src/drtsans/iq.py
2. **Add input validation** - Check for valid Q-range and intensity arrays
3. **Integrate with existing workflow** - Update dataobjects and API

### Testing (2 steps)
4. **Add unit tests** - Test normal cases, edge cases, error conditions in tests/unit/
5. **Add integration test** - Test end-to-end workflow in tests/integration/

### Documentation (1 step)
6. **Update docstrings and docs** - Add usage examples and update API reference

**Testing approach**: After steps 1-3, run step 4 to verify. Step 5 verifies the complete feature.
```

### Step 3: Implement Incrementally
- Complete ONE item from the plan at a time
- Show the user what you're doing: "Implementing step 1: Add new binning function"
- After each significant step, pause if appropriate

### Step 4: Test
After implementation:
- Run relevant tests with pytest
- Show test results to the user
- Fix any failures before proceeding
- Check test coverage with pytest-cov

### Step 5: Review (For Major Changes)
After completing a significant feature or refactor:
- Invoke the **Plan** agent (`run_subagent` with `agentName: "Plan"`) to review the changes
- Provide the review context: what changed and why
- Address any issues found
- Update documentation as needed

## 🛠️ Technology Stack

This project uses specific technologies and frameworks:

### Core Framework
- **Mantid**: Framework for neutron and muon data reduction and analysis
  - Use `mantid.simpleapi` for algorithms
  - Use `mantid.kernel` for logging and configuration
  - Use `mantid.api` for workspace types

### Data Processing
- **numpy**: Array operations and numerical computing
- **pandas**: Data manipulation and analysis
- **h5py**: HDF5 file handling for NeXus files
- **lmfit**: Curve fitting (version 1.3.3)
- **matplotlib**: Plotting and visualization
- **mpld3**: Interactive matplotlib plots

### Development Tools
- **Testing**: pytest with plugins (pytest-cov, pytest-qt, pytest-mock, pytest-xvfb, pytest-xdist)
- **Linting/Formatting**: ruff (configured in pyproject.toml)
- **Documentation**: Sphinx with RTD theme
- **Pre-commit hooks**: Configured in .pre-commit-config.yaml
- **Package management**: pixi (conda-based environment management)
- **Version control**: versioningit for automatic versioning from git tags

### Project Structure
```
drtsans/
├── src/drtsans/          # Main package source
├── tests/                # Test suite
│   ├── unit/            # Unit tests
│   ├── integration/     # Integration tests
│   └── conftest.py      # Pytest fixtures
├── docs/                 # Sphinx documentation
├── notebooks/            # Jupyter notebooks for examples
├── scripts/              # Utility scripts
└── data/                 # Test data references
```

## 📝 Code Quality Standards

### Always Include:
1. **Type hints** - For function parameters and return values (already in use)
2. **Docstrings** - For all public functions/classes (numpy/sphinx style preferred)
3. **Error handling** - Use Mantid logger for warnings and errors
4. **Input validation** - Check assumptions about data shapes and types
5. **Comments** - Explain scientific methodology and "why", not "what"

### SANS-Specific Considerations:
- **Units**: Always document physical units (e.g., Q in 1/Å, I in 1/cm)
- **Mantid workspaces**: Use proper workspace types (MatrixWorkspace, IEventWorkspace)
- **Data formats**: Support NeXus (.nxs.h5), processed data formats
- **Instrument specifics**: Code may handle EQSANS, GPSANS, BIOSANS

### Code Example Template:
```python
from typing import Optional, Union
import numpy as np
from mantid.api import MatrixWorkspace
from mantid.kernel import Logger

logger = Logger("drtsans.module_name")


def process_sans_data(
    input_workspace: Union[str, MatrixWorkspace],
    q_min: float = 0.001,
    q_max: float = 0.5,
    mask_edges: bool = True
) -> MatrixWorkspace:
    """
    Process SANS data and calculate I(Q) profile.

    This function bins scattered intensity as a function of momentum transfer Q,
    applying proper error propagation and optional edge masking for detector panels.

    Parameters
    ----------
    input_workspace : str or MatrixWorkspace
        Input workspace containing reduced SANS data
    q_min : float, optional
        Minimum Q value in 1/Å (default: 0.001)
    q_max : float, optional
        Maximum Q value in 1/Å (default: 0.5)
    mask_edges : bool, optional
        Whether to mask detector edges (default: True)

    Returns
    -------
    MatrixWorkspace
        Binned I(Q) workspace with proper units

    Raises
    ------
    ValueError
        If Q range is invalid or workspace has incompatible units

    Example
    -------
    >>> ws = process_sans_data("EQSANS_68200", q_min=0.005, q_max=0.3)
    >>> print(ws.readY(0))  # Get intensity values

    Notes
    -----
    - Q values are momentum transfer: Q = 4π sin(θ)/λ
    - Intensities are in absolute units (1/cm) after proper normalization
    """
    # Validate inputs
    if q_min <= 0 or q_max <= q_min:
        raise ValueError(f"Invalid Q range: [{q_min}, {q_max}]")

    # Convert to workspace if string
    if isinstance(input_workspace, str):
        input_workspace = mtd[input_workspace]

    logger.information(f"Processing SANS data: Q range [{q_min}, {q_max}] 1/Å")

    # Implementation...

    return output_workspace
```

## 🧪 Testing Guidelines

### Test Structure
Every test should follow Arrange-Act-Assert:
```python
def test_bin_intensity_into_q1d_normal_case(clean_workspace):
    """Test I(Q) binning with typical SANS data."""
    # Arrange: Set up test workspace and expected results
    input_ws = clean_workspace("EQSANS_68200")
    q_min, q_max = 0.005, 0.3
    expected_bins = 50

    # Act: Call the function
    result = bin_intensity_into_q1d(input_ws, qmin=q_min, qmax=q_max, bins=expected_bins)

    # Assert: Verify expectations
    assert result.intensity.shape[0] == expected_bins
    assert np.all(result.mod_q >= q_min)
    assert np.all(result.mod_q <= q_max)
    assert not np.any(np.isnan(result.intensity))
```

### Test Organization
- **Unit tests** (`tests/unit/`) - Test individual functions
- **Integration tests** (`tests/integration/`) - Test complete workflows
- **Test data** - Referenced from `/SNS/EQSANS/shared/sans-backend/data/` or `tests/data/`
- **Fixtures** - Use pytest fixtures defined in `tests/conftest.py`

### Test Markers
Use pytest markers for special test categories:
```python
@pytest.mark.datarepo
def test_with_real_data():
    """Test using drtsans-data repository."""
    pass

@pytest.mark.mount_eqsans
def test_with_sns_mount():
    """Test requiring /SNS/EQSANS/ data mount."""
    pass
```

### Running Tests
```bash
# All tests
pixi run test

# Unit tests only
pixi run unit-test

# Integration tests
pixi run integration-test

# With coverage
pytest --cov=src --cov-report=term-missing
```

## 🔍 Code Review Process

### When to Trigger Review:
- After implementing a new reduction algorithm
- After significant refactoring of core modules
- Before marking major work as complete
- When requested by user

### Review Checklist:
- [ ] Code follows ruff linting rules
- [ ] All functions have type hints and docstrings
- [ ] Tests cover normal, edge, and error cases
- [ ] Mantid algorithms are used correctly
- [ ] Physical units are documented and correct
- [ ] Error messages use Mantid logger
- [ ] No performance regressions for large datasets
- [ ] Documentation is updated in docs/
- [ ] Works with all supported instruments (EQSANS, GPSANS, BIOSANS)

## 📚 Documentation Standards

### Module Docstrings:
```python
"""
Module for SANS data I(Q) binning operations.

This module provides functions for binning scattered neutron intensity as a
function of momentum transfer Q, including:

- 1D radial binning: bin_intensity_into_q1d()
- 2D Cartesian binning: bin_intensity_into_q2d()
- Wedge selection: select_i_of_q_by_wedge()
- Annular binning: bin_annular_into_q1d()

Typical workflow:
    1. Load reduced SANS data
    2. Choose binning method and parameters
    3. Apply binning with proper error propagation
    4. Save results in CanSAS or ASCII format

Example:
    >>> from drtsans import iq
    >>> result = iq.bin_intensity_into_q1d(workspace, qmin=0.005, qmax=0.3)
"""
```

### Function Docstrings:
Use numpy/sphinx style (already used in the codebase):
```python
def bin_intensity_into_q1d(input_workspace, qmin, qmax, bins=100):
    """
    Bin SANS intensity into 1D I(Q) profile.

    Parameters
    ----------
    input_workspace : str or MatrixWorkspace
        Input workspace containing Q values and intensities
    qmin : float
        Minimum Q value in 1/Å
    qmax : float
        Maximum Q value in 1/Å
    bins : int, optional
        Number of Q bins (default: 100)

    Returns
    -------
    IQmod
        Named tuple containing mod_q, intensity, error, delta_mod_q
    """
```

## 🚨 Common Scenarios and Responses

### Scenario: User asks to add a new reduction feature

**Response Pattern:**
```markdown
Let me assess the current reduction pipeline.

[Read relevant files in src/drtsans/]

**Current state:**
- Main reduction API is in src/drtsans/api.py
- I(Q) calculation is in src/drtsans/iq.py
- Uses Mantid algorithms for low-level operations
- Tests exist in tests/unit/ and tests/integration/

**Itemized Plan:**
1. Implement new function in appropriate module
2. Add Mantid algorithm calls as needed
3. Add proper error handling with Mantid logger
4. Add unit tests in tests/unit/
5. Add integration test with realistic data
6. Update API documentation in docs/

Proceeding with implementation...

**Step 1: [Description]**
[Make change, show code with proper docstrings]

**Step 2: [Description]**
[Make change, ensure Mantid compatibility]

**Testing:**
[Run tests with pytest]

Completed! Would you like me to:
- Review the code with the Plan agent?
- Add example notebook in notebooks/?
- Update the user documentation?
```

### Scenario: User reports incorrect I(Q) calculation

**Response Pattern:**
```markdown
Let me investigate the I(Q) calculation.

[Read src/drtsans/iq.py and related modules]

**Assessment:**
- Issue location: src/drtsans/iq.py, line XXX
- Root cause: [Explanation of physics/algorithm issue]
- Impact: Affects 1D radial binning for all instruments

**Itemized Fix Plan:**
1. Correct the momentum transfer calculation
2. Update error propagation
3. Add regression test with known good data
4. Verify against reference results

**Implementing fix:**
[Make changes with proper scientific justification]

**Testing:**
[Run unit and integration tests]

Fixed! The calculation now correctly accounts for [scientific detail].
Added test case to prevent regression.
```

## ⚡ Efficiency Guidelines

1. **Parallel operations** - Use pytest-xdist for parallel test execution
2. **Reuse fixtures** - Leverage existing fixtures in tests/conftest.py
3. **Mantid workspace management** - Clean up temporary workspaces
4. **Avoid redundant I/O** - Cache test data when possible
5. **Use pixi tasks** - Leverage pre-configured pixi tasks for common operations

## 🔗 SANS-Specific Integration Patterns

### Loading Data Pattern:
```python
from drtsans.load import load_events
from drtsans.path import abspath

# Load SANS data with proper path resolution
data_file = abspath("EQSANS_68200", instrument="EQSANS")
workspace = load_events(data_file, output_workspace="raw_data")
```

### Reduction Workflow Pattern:
```python
from drtsans import iq, sensitivity, solid_angle_correction
from mantid.simpleapi import mtd

# Typical SANS reduction workflow
ws = load_events("EQSANS_68200")
ws = solid_angle_correction(ws)
ws = apply_sensitivity_correction(ws, sensitivity_file="sensitivity.nxs")
result = iq.bin_intensity_into_q1d(ws, qmin=0.005, qmax=0.3)

# Save results
save_ascii_binned_1D("output.dat", "EQSANS 68200", result)
```

### Working with Mantid Pattern:
```python
from mantid.simpleapi import mtd, DeleteWorkspace
from mantid.kernel import Logger

logger = Logger("drtsans.mymodule")

def process_data(workspace_name):
    """Process workspace and clean up."""
    try:
        ws = mtd[workspace_name]
        # Process workspace
        result = do_processing(ws)
        logger.information(f"Processed {workspace_name}")
        return result
    finally:
        # Clean up temporary workspaces
        if "__temp" in workspace_name and mtd.doesExist(workspace_name):
            DeleteWorkspace(workspace_name)
```

## 📋 Quick Reference

### Assessment Checklist:
- [ ] Read relevant source files in src/drtsans/
- [ ] Check existing tests in tests/unit/ and tests/integration/
- [ ] Review Mantid algorithm usage
- [ ] Check for instrument-specific code
- [ ] Review documentation in docs/
- [ ] Identify affected components


### Testing Checklist:
- [ ] Add unit tests in tests/unit/
- [ ] Add integration tests if needed
- [ ] Use appropriate pytest markers
- [ ] Test with realistic SANS data
- [ ] Run tests with pixi run test
- [ ] Check coverage with pytest-cov
- [ ] Verify with all instruments if applicable

### Review Checklist:
- [ ] Invoke the Plan agent for major changes
- [ ] Provide clear context
- [ ] Address found issues
- [ ] Update documentation in docs/
- [ ] Verify Mantid compatibility
- [ ] Check scientific correctness
- [ ] Confirm with user

## 🎓 SANS-Specific Guidance

Remember: Users are scientists working with neutron scattering data.

### Always:
- Explain scientific concepts when relevant (Q-space, scattering intensity, etc.)
- Reference physical units (Q in 1/Å, I in 1/cm, wavelength in Å)
- Consider instrument-specific behavior (EQSANS, GPSANS, BIOSANS)
- Use proper Mantid workspace handling
- Explain data reduction steps and their purpose

### Example Explanations:
```
"I'm adding binning in Q-space, where Q is the momentum transfer
(Q = 4π sin(θ)/λ). This converts the 2D detector image into a 1D
I(Q) profile that's easier to analyze and fit."

"Using Mantid's Integration algorithm here to sum counts across
wavelength bins while properly propagating uncertainties."

"The sensitivity correction accounts for pixel-to-pixel detector
efficiency variations, ensuring accurate absolute intensities."
```

---

**Remember**: Every interaction should leave the user with working, tested,
scientifically correct code that properly handles SANS data reduction
workflows using the Mantid framework.
