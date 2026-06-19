---
name: review-comment-audit-report
description: Use when the user wants a detailed audit of whether pull request review comments were addressed, especially with one section per reviewer comment, before/after code snippets, current file and line references, and resolution status.
---

# Review Comment Audit Report

Use this skill when asked to verify whether a pull request author addressed review comments and the user wants a detailed written report.

## Workflow

1. Resolve the pull request.
   - If the user gives a pull request URL or number, use it.
   - Otherwise infer it from the current branch when possible.
   - Identify the reviewer whose comments should be audited. If unspecified, use the authenticated user or ask.

2. Fetch thread-aware review data.
   - Use GitHub review threads, not only flat comments.
   - Preserve the comment author, file path, original line, current line if available, `isResolved`, `isOutdated`, and replies in the thread.
   - Do not assume a resolved thread means the request was fully implemented.

3. Identify likely addressing commits.
   - Look for commits after the review comments.
   - Prefer commit messages such as `address comments`, `fix review comments`, or similar.
   - Use `git show --unified=8 <commit> -- <file>` to capture before/after code.
   - Use `nl -ba` or equivalent on the relevant before and after file revisions so snippets include line numbers:
     - `Before` snippets use line numbers from the revision before the addressing change.
     - `After` snippets use line numbers from the revision after the addressing change, or the final inspected PR/head state.
   - Inspect the current file with `nl -ba <file>` to provide current line-number links.

4. Inspect current code or docs.
   - Confirm whether the final branch state still contains the addressing change.
   - If the current state differs from the addressing commit, mention both.

5. Write one section per reviewer comment.
   - Include the comment number in the section title.
   - Use a short descriptive title, not only the file name.
   - Keep the explanation grounded in concrete code or documentation changes.

6. Output the report.
   Mandatory order:
   1. First, post the full report body in the chat.
   2. Then save the exact same report body to a Markdown file under `/tmp`.
   3. Finally, tell the user the saved file path.

   Do not replace the chat report with a summary, excerpt, or file link.
   Do not rely on a concise final answer when this skill is active; the full report is the final answer.


## Section Format

Use this structure for each comment:

````markdown
**N. Short Title**

Your comment:
Briefly quote or paraphrase the reviewer comment.

Status:
Addressed / Partially addressed / Not addressed / Resolved by discussion.

Where:
- [file.py](/abs/path/file.py:123)
- GitHub thread state: resolved/unresolved, outdated/current.

Before:
```text
123  old code
124  old code
```

After:
```text
145  new code
146  new code
```

Explanation:
Explain exactly how the change responds to the comment. If the thread is marked resolved but the code only partially implements the request, say that clearly.
````

## Reporting Rules

- The complete report must appear in the chat, even if it is long.
- A saved `/tmp/*.md` file is required but is not a substitute for posting the report.
- The final response should contain the full report first, followed by the saved file path.
- Do not use a table for detailed audits unless the user explicitly requests one.
- Use one numbered section per reviewer comment.
- Include file names and actively look up the current line number for every reference.
- Before/after code snippets must include line numbers inside the snippet.
- `Before` snippet line numbers must come from the file revision before the addressing change.
- `After` snippet line numbers must come from the file revision after the addressing change, or from the final inspected PR/head state when that is what the report is evaluating.
- Do not use `diff` fences for numbered snippets unless the line numbers remain clear; prefer `text` fences for numbered before/after excerpts.
- Show before/after code snippets when a code or documentation change exists.
- If the author addressed the comment only in discussion, say `resolved by discussion, not by code change`.
- If a thread is resolved but the implementation is incomplete, say `partially addressed` and explain the gap.
- If no code change is visible for a comment, state that directly.
- Do not resolve threads, reply on GitHub, submit reviews, or mutate the pull request unless explicitly asked.

## Good Output Characteristics

- The report should be easy to scan but not compressed into a table.
- Each section should stand alone: a reader should understand the original request, the implementation response, and any remaining gap.
- Snippets should be small enough to focus on the relevant change but complete enough to prove the conclusion.
- Prefer precise language such as `addressed`, `partially addressed`, `not addressed`, or `resolved by discussion`.
