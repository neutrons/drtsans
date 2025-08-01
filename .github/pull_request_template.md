## Description of work:

Check all that apply:
- [ ] updated documentation and checked that it looks correct in the [pull request preview](https://docs.readthedocs.com/platform/stable/pull-requests.html)
- [ ] Source added/refactored
- [ ] Added unit tests
- [ ] Added integration tests
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
It is strongly encouraged that the reviewer runs the following tests in their local machine
because these tests are not run by the GitLab CI. It is assumed that the reviewer has the /SNS and /HFIR filesystems
remotely mounted in their machine.

```bash
cd /path/to/my/local/drtsans/repo/
git fetch origin merge-requests/<MERGE_REQUEST_NUMBER>/head:mr<MERGE_REQUEST_NUMBER>
git switch mr<MERGE_REQUEST_NUMBER>
conda activate <my_drtsans_dev_environment>
pytest -m mount_eqsans ./tests/unit/ ./tests/integration/
```
In the above code snippet, substitute `<MERGE_REQUEST_NUMBER>` for the actual merge request number. Also substitute
`<my_drtsans_dev_environment>` with the name of the conda environment you use for development. It is critical that
you have installed the repo in this conda environment in editable mode with `pip install -e .` or `conda develop .`
