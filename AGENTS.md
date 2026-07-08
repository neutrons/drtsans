@.github/copilot-instructions.md

## Personal Codex Context

At the start of work in this repository, inspect non-sensitive auxiliary Codex
context under `~/.codex/` when filesystem permissions allow it.

Prefer these locations:

- `~/.codex/rules/*.rules`
- `~/.codex/config.toml`
- `~/.codex/skills/*/SKILL.md`
- `~/.codex/agents/**`
- `~/.codex/plugins/**/skills/**/SKILL.md`

Do not read or summarize sensitive runtime state unless the user explicitly
asks:

- `~/.codex/auth.json`
- `~/.codex/history.jsonl`
- `~/.codex/sessions/**`
- `~/.codex/log/**`
- `~/.codex/*.sqlite*`
- `~/.codex/shell_snapshots/**`

Summarize which auxiliary files were loaded and apply only instructions relevant to the current task.
