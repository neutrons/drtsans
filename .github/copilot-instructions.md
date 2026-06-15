## Workflow

- Assess before acting: read relevant source, tests, conventions, and Mantid dependencies.
- For non-trivial work, provide a concise itemized plan with testable steps.
- Implement incrementally and keep the user informed about progress and next steps.
- When a request is ambiguous, inspect existing patterns first; if still unclear, state the assumption and proceed.
- Test the smallest relevant behavior first, then broaden testing when needed.
- Ask permission before invoking a review agent for major features or refactors.

## Scope Rules

- Do not refactor unrelated code without permission.
- Do not add unrequested features, parameters, or abstractions; propose them separately when useful.
- Do not guess file paths, workspace names, or APIs.
- Do not use `print()` for runtime messages; use the Mantid logger.

## Project Context

- Users are scientists working with neutron scattering data.
- Explain scientific concepts when relevant, such as Q-space, scattering intensity, and data reduction steps.
- Account for EQSANS, GPSANS, and BIOSANS behavior.
- Document physical units, including Q in 1/Å, intensity in 1/cm, and wavelength in Å.

## Technology Stack

- **Mantid**: use `mantid.simpleapi` for algorithms, `mantid.kernel` for logging/configuration, and `mantid.api` for workspace types.
- **Testing**: pytest with pytest-cov, pytest-qt, pytest-mock, pytest-xvfb, and pytest-xdist.
- **Linting/formatting**: ruff, configured in `pyproject.toml`.
- **Documentation**: Sphinx with RTD theme.
- **Package management**: pixi.
- **Versioning**: versioningit from git tags.

## Project Structure

```text
drtsans2/
├── src/drtsans/          # Main package source
├── tests/                # Test suite
│   ├── unit/             # Unit tests
│   ├── integration/      # Integration tests
│   └── conftest.py       # Pytest fixtures
├── docs/                 # Sphinx documentation
├── notebooks/            # Jupyter notebooks for examples
├── scripts/              # Utility scripts
└── data/                 # Test data references
```

## Code Standards

- Add type hints for function parameters and return values.
- Add numpy/sphinx-style docstrings for public functions and classes.
- Validate assumptions about data shapes, types, units, and workspace types.
- Use Mantid logging for warnings and errors.
- Explain scientific methodology or intent in comments, not obvious mechanics.
- Use proper Mantid workspace types, such as `MatrixWorkspace` and `IEventWorkspace`.
- Support NeXus (`.nxs.h5`) and processed data formats where relevant.

## Pixi
Run Pixi tasks with permission to write to the repository and Pixi/cache directories.
In Codex read-only sandboxes, request sandbox escalation with `prefix_rule: ["pixi", "run"]`.

## Testing

- Follow Arrange-Act-Assert.
- Place unit tests in `tests/unit/` and integration tests in `tests/integration/`.
- Use fixtures from `tests/conftest.py` and existing pytest markers.
- Reference test data from `/SNS/EQSANS/shared/sans-backend/data/` or `tests/data/`.
- Clean up temporary Mantid workspaces.
- Cache test data where practical to avoid redundant I/O.

Run tests with:

```bash
pixi run unit-test
pixi run integration-test
pixi run test
```

## Review

Ask permission before review. Use review for new reduction algorithms, significant core refactors, major work before completion, or explicit user requests.

When approved, invoke the **Plan** agent with `run_subagent` and `agentName: "Plan"`, include what changed and why, then address findings.

Review checklist:

- Ruff rules pass.
- Public functions have type hints and docstrings.
- Tests cover normal, edge, and error cases.
- Mantid algorithms and workspace types are used correctly.
- Physical units are documented and correct.
- Warnings and errors use Mantid logging.
- Large datasets do not suffer obvious performance regressions.
- Documentation is updated when behavior changes.
- Supported instruments still work: EQSANS, GPSANS, and BIOSANS.
