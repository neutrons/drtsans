## Description of work:

**Check all that apply:**
- [ ] updated documentation and checked that it looks correct in the [pull request preview](https://docs.readthedocs.com/platform/stable/pull-requests.html)
- [ ] Source added/refactored
- [ ] Unit tests added/refactored
- [ ] Integration tests added/refactored
- [ ] Included a manual test for the reviewer
- [ ] Verified that tests requiring the /SNS and /HFIR filesystems pass without fail

**References:**
- Links to IBM EWM items:
- Links to related issues or pull requests:

## :warning: Manual test for the reviewer
<!-- Instructions for testing here. -->

## Check list for the reviewer
- [ ] best software practices
    + [ ] clearly named variables (better to be verbose in variable names)
    + [ ] code comments explaining the intent of code blocks
- [ ] All the tests are passing
- [ ] The documentation is up to date and looks correct in the [pull request preview](https://docs.readthedocs.com/platform/stable/pull-requests.html)
- [ ] code comments added when explaining intent

### Execution of tests requiring the /SNS and /HFIR filesystems
It is strongly encouraged that the reviewer runs the "manual" tests in their local machine
because these are not run by the GitLab CI. It is necessary that the remote /SNS and /HFIR filesystems
are mounted in the machine running the tests.
In the below code block, substitute `<MERGE_REQUEST_NUMBER>` for the actual merge request number`


```bash
cd /path/to/my/drtsans/
git fetch origin merge-requests/<MERGE_REQUEST_NUMBER>/head:mr<MERGE_REQUEST_NUMBER>
git switch mr<MERGE_REQUEST_NUMBER>
pixi run manual-test
```
