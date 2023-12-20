## Description of work:

- [ ] Source added/refactored?
- [ ] Added unit tests?
- [ ] Added integration tests?
- [ ] Verified that tests requiring the /SNS and /HFIR filesystems do pass?

## To Test by Reviewer:

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

### Manual test for this merge request

(here the author of the merge request can include a test specific to this merge request)

### Review best practices

The reviewer should check for the discussed software practices

- [ ] all internal functions should have an underbar, as is python standard
- [ ] need more comments in meat of code explaining intent
- [ ] clearly named variables (better to be verbose in variable names)
