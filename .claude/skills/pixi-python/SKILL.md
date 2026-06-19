---
name: pixi-python
description: Use when working in a Python repository managed by Pixi, including projects with pixi.toml, pixi.lock, .pixi/, or [tool.pixi] sections in pyproject.toml.
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
