# Pixi

Run Pixi tasks with permission to write to the repository and Pixi/cache directories.
In Codex read-only sandboxes, request sandbox escalation with `prefix_rule: ["pixi", "run"]`.

# Testing

Run unit tests with `pixi run unit-test`.
Run integration tests with `pixi run integration-test`.
Run the full test suite with `pixi run test`.
