---
name: defect-review
description: Review a Defect (bug report) before Development Review and estimate the fix effort in resource hours. Use when the user describes a bug/Defect and wants a Development-Review-style scoping pass, a breakdown into Task items with hour estimates, or asks "how many hours to fix this" / "estimate this defect/bug". Not for Stories (use story-review instead — Stories are estimated in Story points, Defects are estimated directly in hours).
---

# Review Defect

## Overview

This skill prepares one reviewer-style comment for the Neutron Data Project
(NDP) "Development Review" process for a single **Defect**: read the Defect,
evaluate it against the Defect-specific criteria checklist, break it into
elemental **Task** items with their own hour estimates, and roll those up into
a total resource-hour estimate. It is for **Defects only** — Defects are
estimated directly in resource hours, never in Story points. For a new
Story/feature/enhancement, use the `story-review` skill instead.

Do not confuse a Task with a Defect. A Defect is broken down into Task items
the way a project is broken down into a to-do list. A Task is a concrete
implementation step; it is never itself a Story or a Defect.

## Input

A free-text description of the Defect to be fixed. This may include:

- The symptom/problem being reported.
- Reproduction steps, if known.
- The expected behavior or correct result.
- Any related Story, parent item, or context.
- Additional files, such as logs, error output, screenshots, or code
  snippets, that are relevant to the Defect.

If the input is missing information needed to complete a criteria item,
**do not invent it** — surface it as an open question instead (see Workflow).

## Workflow

1. **Review related context.** Read any parent items, related Stories or
   Defects, and associated requirements supplied with or linked from the
   Defect. If relevant context is missing, note it as an open question.

2. **Restate understanding.** Summarize the Defect in 1–2 sentences, in your
   own words, so the reviewer/owner can confirm the problem was understood
   correctly.

3. **Walk the Defect Review Criteria checklist.** For each item below, mark
   it `Clear`, `Unclear`, or `N/A`, with a one-line justification:
   - **Clarity** — is the Defect clear and understood?
   - **Reproduction** — are there instructions for reproducing the problem?
   - **Expected behavior** — is it clear how the software should behave, or
     what result it should produce?
   - **Equations & formulae** — are all equations/formulae given explicitly,
     or referenced to a citable source (e.g., a paper)? Existing code is
     *not* an acceptable reference for a formula.
   - **Data for verification** — is data provided to verify the fix against?
     Output from other software is *not* generally acceptable as
     verification data.
   - **Inputs & outputs** — is input test data available where applicable,
     and is example output data available to compare against?
   - **APIs** — are explicit parameters and return values identified and
     characterized, where applicable?
   - **UIs** — are fields, file characteristics, and visual/interaction
     features clear, where applicable? Wireframes, mock-ups, or sketches are
     a good sign here, not a requirement.

   Anything marked `Unclear` becomes an **open question** — a specific
   question you would ask the Defect owner, not a guess you make on their
   behalf.

4. **State whether the fix is known.** Reproducibility is not the same as
   having a known fix. Explicitly note: is the fix approach already
   understood, or does it first require investigation? If the fix is
   unknown, add a dedicated investigation Task ("determine the fix
   approach") as its own Task, estimated separately. In that case, make clear
   that the estimate may cover investigation only; the actual fix Task(s)
   should generally not be scheduled or estimated with confidence until the
   approach is understood.

5. **Decompose into Task items.** List the elemental Task items required to
   address the Defect — concrete, individually completable steps, similar
   to a to-do list. Per NDP guidance, size each Task to roughly 4–16 hours;
   if a Task looks larger than that, split it into smaller Tasks rather than
   leaving it oversized.

6. **Estimate each Task in hours.** Give every individual Task item its own
   man-hour estimate (a single number or a tight range) — not just a total
   for the whole Defect. This is what lets a reviewer sanity-check the
   estimate task by task.

7. **Roll up to a total.** Sum the per-Task hour estimates into a total
   resource-hour estimate for the Defect. Defects are **never** converted to
   Story points — report the total in hours only.

8. **Save the review to a file.** Write the full `REVIEW` output (see Output
   format below) to a markdown file:
   - If the Defect number is known (given explicitly, or evident from the
     input), save to `/tmp/review_<DefectNumber>.md`.
   - If no Defect number is available, save to `/tmp/review_defect.md`.
   - Wrap all lines in the saved file so none exceeds 120 characters.

## Output format

```
REVIEW

Understanding: <1-2 sentence restatement of the Defect>

Criteria:
- Clarity: Clear/Unclear — <justification>
- Reproduction: Clear/Unclear — <justification>
- Expected behavior: Clear/Unclear — <justification>
- Equations & formulae: Clear/Unclear/N/A — <justification>
- Data for verification: Clear/Unclear/N/A — <justification>
- Inputs & outputs: Clear/Unclear/N/A — <justification>
- APIs: Clear/Unclear/N/A — <justification>
- UIs: Clear/Unclear/N/A — <justification>

Fix known: Yes/No — <note; if No, an investigation Task is included below>

Open questions:
- <question for the Defect owner, if any>

Tasks:
- <Task 1 description> (est. <H1> hrs)
- <Task 2 description> (est. <H2> hrs)
- ...

Estimate: <sum of task hours> resource hours
```

The `REVIEW` prefix matches the NDP convention for review comments left in a
work item's discussion field, so this output can be pasted there directly.
This exact content is also what gets written to the markdown file described
in Workflow step 8.

## Rules

- Do not fabricate reproduction steps, expected behavior, equations, or
  inputs/outputs that the prompt did not supply — raise them as open
  questions instead.
- This skill is for Defects, not Stories. Stories are estimated in Story
  points using a different criteria set — use `story-review` for those.
- Never convert a Defect's hour estimate into Story points.
- If the fix approach is unknown, call that out explicitly and estimate the
  investigation separately. Do not present a confident final fix estimate
  until the fix approach is understood.
- Keep the criteria list above as the source of truth for this skill; if it
  changes, update this file deliberately rather than drifting from it
  silently.
- No line in the saved review file should exceed 120 characters; wrap prose
  to fit.
