---
name: story-review
description: Review a new User Story before Development Review and estimate its effort in Story points. Use when the user describes a new Story/feature/enhancement and wants a Development-Review-style scoping pass, a breakdown into Task items with hour estimates, or asks "how many story points is this" / "estimate this story". Not for Defects (defects are estimated in resource hours, not points, and follow a different criteria set).
---

# Story Review

## Overview

This skill reproduces the Neutron Data Project (NDP) "Development Review"
process for a single **User Story**: read the Story, evaluate it against a
fixed criteria checklist, break it into elemental **Task** items with their
own hour estimates, and roll those up into a Story point estimate. It is for
**Stories only** — Defects are estimated in resource hours directly and use a
different (though similar) criteria set; do not use this skill to estimate a
Defect.

Do not confuse a Task with a Story. A Story (or Defect) is broken down into
Task items the way a project is broken down into a to-do list. A Task is a
concrete implementation step; it is never itself a Story or a Defect, and it
is never re-estimated in Story points — only in hours.

## Input

A free-text description of the Story to be implemented. This may include:

- The Story text itself (as a user-facing requirement).
- Any parent Story, Capability, or related requirement given for context.
- Existing acceptance criteria, equations, sample data, API shapes, or
  wireframes the requester has already provided.
- Additional files, such as images, PDFs, or code snippets, that are relevant to the
  Story.

If the input is missing information needed to complete a criteria item,
**do not invent it** — surface it as an open question instead (see Workflow).

## Workflow

1. **Restate understanding.** Summarize the Story in 1–2 sentences, in your
   own words, so the reviewer/owner can confirm scope was understood
   correctly.

2. **Walk the User Story Criteria checklist.** For each item below, mark it
   `Clear`, `Unclear`, or `N/A`, with a one-line justification:
   - **Clarity** — is the Story clear and understood?
   - **Acceptance criteria** — is it clear enough that a test could be
     written to verify the implementation passes it?
   - **Equations & formulae** — are all equations/formulae given explicitly,
     or referenced to a citable source (e.g., a paper)? Existing code is
     *not* an acceptable reference for a formula.
   - **Data for verification** — is data provided to verify the
     implementation against? Output from other software is *not* generally
     acceptable as verification data.
   - **Inputs & outputs** — is input test data available where applicable,
     and is example output data available to compare against?
   - **APIs** — are explicit parameters and return values identified and
     characterized, where applicable?
   - **UIs** — are fields, file characteristics, and visual/interaction
     features clear, where applicable? Wireframes, mock-ups, or sketches are
     a good sign here, not a requirement.

   Anything marked `Unclear` becomes an **open question** — a specific
   question you would ask the Story owner, not a guess you make on their
   behalf.

3. **Decompose into Task items.** List the elemental Task items required to
   implement the Story — concrete, individually completable steps, similar
   to a to-do list. Per NDP guidance, size each Task to roughly 4–16 hours;
   if a Task looks larger than that, split it into smaller Tasks rather than
   leaving it oversized.

4. **Estimate each Task in hours.** Give every individual Task item its own
   man-hour estimate (a single number or a tight range) — not just a total
   for the whole Story. This is what lets a reviewer sanity-check the
   estimate task by task.

5. **Roll up to a total and convert to Story points.** Sum the per-Task hour
   estimates into a total-hours estimate for the Story, then map that total
   to Story points using the table below.

6. **Flag Epics.** If the resulting point estimate exceeds 8, mark the Story
   as an **EPIC** and propose a concrete decomposition into child Stories
   (each targeted at ≤ 8 points), following the NDP pattern: the parent
   Story carries no functionality of its own — its points are simply the
   sum of its children's points. This Story-level split (into smaller
   Stories) is a distinct step from the Task-level breakdown in step 3, and
   is only needed when the total estimate is too large for a single Story.

7. **Save the review to a file.** Write the full `REVIEW` output (see Output
   format below) to a markdown file:
   - If the Story number is known (given explicitly, or evident from the
     input), save to `/tmp/review_<StoryNumber>.md`.
   - If no Story number is available, save to `/tmp/review_story.md`.
   - Wrap all lines in the saved file so none exceeds 120 characters.

## Story point conversion table

| Story points | Hours range |
|---|---|
| 1 | under 2 hours |
| 2 | 2 to 8 hours |
| 3 | 8 to 20 hours |
| 5 | 20 to 40 hours |
| 8 | 40 to 80 hours |
| > 8 | **EPIC** — recommend splitting into sub-stories, each estimated at or under 8 points |

## Output format

```
REVIEW

Understanding: <1-2 sentence restatement of the Story>

Criteria:
- Clarity: Clear/Unclear — <justification>
- Acceptance criteria: Clear/Unclear — <justification>
- Equations & formulae: Clear/Unclear/N/A — <justification>
- Data for verification: Clear/Unclear/N/A — <justification>
- Inputs & outputs: Clear/Unclear/N/A — <justification>
- APIs: Clear/Unclear/N/A — <justification>
- UIs: Clear/Unclear/N/A — <justification>

Open questions:
- <question for the Story owner, if any>

Tasks:
- <Task 1 description> (est. <H1> hrs)
- <Task 2 description> (est. <H2> hrs)
- ...

[EPIC — recommend splitting into sub-stories: <proposed child stories>]  (only if > 8 pts)

Total task estimate: <sum of task hours> hours
Estimate: <N> pts
```

The `REVIEW` prefix matches the NDP convention for review comments left in a
work item's discussion field, so this output can be pasted there directly.
This exact content is also what gets written to the markdown file described
in Workflow step 7.

## Rules

- Do not fabricate equations, inputs/outputs, or API contracts that the
  prompt did not supply — raise them as open questions instead.
- This skill is for Stories, not Defects. Defects are estimated directly in
  resource hours based on reviewer experience, not Story points, and use a
  Defect-specific criteria set (reproduction steps, expected behavior).
- Keep the criteria list and the point/hour table above as the source of
  truth for this skill; if either changes, update this file deliberately
  rather than drifting from it silently.
- No line in the saved review file should exceed 120 characters; wrap prose
  to fit.
