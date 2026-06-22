---
name: pixi-python
description: Use whenever running Python tooling in a Pixi-managed repository, including python, pytest, ruff, mypy, pip, or similar commands; this environment skill applies alongside task-specific skills when pixi.toml, pixi.lock, .pixi/, or [tool.pixi] sections in pyproject.toml are present.
---

# Pixi Python Repositories

Before running Python tooling, check whether the repository is managed by Pixi.

Treat the repo as Pixi-managed if any of these are present:

- `pixi.toml`
- `pixi.lock`
- `.pixi/`
- `pyproject.toml` containing `[tool.pixi...]`

Use Pixi for Python commands:

```bash
pixi run python ...
pixi run pytest ...
pixi run ruff ...
pixi run mypy ...
```

Do not first try bare `python`, `pytest`, `pip`, `ruff`,
or similar tools unless explicitly checking the system environment.

**Important**: If the the sandbox is read-only and `pixi run ...` fails because files must be written under `.pixi/`,
`.pytest_cache/`, or test result paths, rerun the same `pixi run ...` command outside the sandbox.
