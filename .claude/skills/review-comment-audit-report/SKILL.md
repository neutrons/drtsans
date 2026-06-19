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
   - Inspect the current file with `nl -ba <file>` to provide current line-number links.

4. Inspect current code or docs.
   - Confirm whether the final branch state still contains the addressing change.
   - If the current state differs from the addressing commit, mention both.

5. Write one section per reviewer comment.
   - Include the comment number in the section title.
   - Use a short descriptive title, not only the file name.
   - Keep the explanation grounded in concrete code or documentation changes.

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
```diff
- old code
```

After:
```diff
+ new code
```

Explanation:
Explain exactly how the change responds to the comment. If the thread is marked resolved but the code only partially implements the request, say that clearly.
````

## Reporting Rules

- Do not use a table for detailed audits unless the user explicitly requests one.
- Use one numbered section per reviewer comment.
- Include file names and current line numbers for every code reference.
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
