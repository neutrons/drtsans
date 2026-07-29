---
name: assess-reviewer-comments
description: Assess unresolved GitHub pull request review comments from a specified reviewer and PR number. Use when the user wants Codex to collect unresolved PR review threads, judge each comment's merit and current relevance against the PR code, propose resolution actions when needed, and rank comments by importance before implementation.
---

# Assess Reviewer Comments

## Overview

Use this skill to turn unresolved review comments from a specific GitHub reviewer into an ordered action plan. The result should distinguish valid findings from stale, low-impact, or mistaken comments, and should not mutate GitHub state unless the user separately asks for that.

## Workflow

1. Resolve the repository, pull request, and reviewer.
   - Require a PR number or URL and a reviewer GitHub handle.
   - If the repository is not explicit, infer it from the current Git remote.
   - Normalize the reviewer handle by ignoring a leading `@`.
   - If either the PR or reviewer is ambiguous, ask before continuing.

2. Fetch unresolved thread-aware review data.
   - Use GitHub review threads, not only flat PR review comments.
   - Include only threads with `isResolved: false`.
   - Include a thread when the specified reviewer authored the thread's first comment or an unresolved actionable comment within the thread.
   - Preserve each thread's URL, file path, original line, current line if available, `isOutdated`, diff hunk, reviewer comment body, replies, and timestamps.
   - Also read relevant top-level PR conversation comments from that reviewer if they appear to refer to unresolved review work, but separate them from review threads because GitHub does not mark top-level conversation comments as resolved.

3. Inspect the current PR code before judging comments.
   - Read the current file and surrounding code for each comment.
   - Check related tests, docs, or call sites when needed to determine whether the comment is valid.
   - For outdated threads, map the concern to the current code before deciding relevance.
   - Do not assume a comment is correct because it is unresolved; verify it against the code.
   - Do not assume a thread is irrelevant because it is outdated; decide whether the underlying issue still exists.

4. Assess merit and relevance.
   - Merit is whether the reviewer's concern is technically correct, consistent with project conventions, and supported by the current code.
   - Relevance is whether addressing it would improve correctness, maintainability, user-facing behavior, tests, docs, or review acceptance for this PR.
   - Prefer concrete evidence from source, tests, docs, and existing patterns.
   - Mark comments as one of: `High relevance`, `Medium relevance`, `Low relevance`, `Not relevant`, or `Needs clarification`.
   - State when a comment is valid but out of scope for the PR.
   - State when a comment appears mistaken, already handled, superseded, or only a preference.

5. Offer actions.
   - For valid and relevant comments, propose the smallest concrete code, test, docs, or discussion action that would resolve the concern.
   - For partially valid comments, split the actionable part from the non-actionable part.
   - For low-relevance or mistaken comments, propose a concise reply explaining why no code change is recommended.
   - If more information is needed, propose the exact clarification question to ask.
   - Do not implement changes, post replies, submit reviews, or resolve GitHub threads unless explicitly asked.

6. Rank the comments.
   - Rank by expected PR impact, not by chronological order.
   - Prioritize correctness, data loss, failing CI, public API or CLI behavior, security, and reproducibility.
   - Next prioritize maintainability, project conventions, tests, docs, and reviewer-blocking concerns.
   - Deprioritize style preferences, already-obsolete comments, comments outside PR scope, and comments whose requested change would add risk without clear benefit.
   - For ties, put easier low-risk fixes before invasive changes.

## Recommended Data Collection

Use the GitHub app when it exposes the needed PR and comment data. Use `gh` GraphQL when thread-level resolution state or inline review context is required.

Useful fields for review threads:

- `isResolved`
- `isOutdated`
- `path`
- `line`
- `originalLine`
- `diffSide`
- `comments.nodes.author.login`
- `comments.nodes.body`
- `comments.nodes.diffHunk`
- `comments.nodes.url`
- `comments.nodes.createdAt`
- `comments.nodes.replyTo`

Fetch enough PR metadata to identify the base branch, head branch, latest commits, review decision, and whether the reviewer has requested changes.

## Output Format

Start with a compact summary:

```markdown
Found N unresolved review thread(s) from @reviewer on PR #number.
Recommended order: comment A, comment B, comment C.
```

Then provide one section per ranked comment:

```markdown
1. High relevance — Short title

Reviewer comment:
Brief paraphrase or short quote.

Where:
- GitHub: <thread/comment URL>
- File: path/to/file.ext:line

Assessment:
Explain whether the comment is correct and relevant, using concrete source or test evidence.

Recommended action:
Smallest action that resolves the concern, or a reply/clarification if no code change is recommended.
```

End with an implementation plan only if the user asked for next steps or fixes. Otherwise stop after the assessment and action recommendations.

## Rules

- Keep the work read-only by default.
- Do not collapse multiple unrelated comments into one assessment unless they are duplicate concerns.
- Do not rely solely on GitHub resolution state, reactions, or reviewer status.
- Do not characterize a comment as invalid without citing the specific code, tests, docs, or project convention that supports that conclusion.
- Prefer direct, implementable actions over broad advice.
- If no unresolved comments from the reviewer are found, report that and mention whether resolved or outdated comments were seen.
