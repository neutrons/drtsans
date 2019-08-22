drt-sans
========

For end users go to [next version](http://scse-sans-demo.ornl.gov/)
[QA version](http://scse-ui.ornl.gov:8080/)

Use [jupyter](https://jupyter.sns.gov/) to have a play with the
code. The kernel to select is `sans at ...`.


**This is a python3 only package.**

Using the Docker packaged environment
-----------------------------------------------

This the instructions for someone who wants to use the Docker container
created through the automated build pipeline to develop drt-sans, use
drt-sans to develop reduction scripts, or test existing drt-sans
functionality. The SNS analysis cluster does not have Docker installed
and Docker is required to follow these instructions.

1. (If not installed) Install Docker: https://docs.docker.com/install/
2. Download the latest sans-backend-run.sh script from the feature, release, or master branch for which you are testing: http://scse-mantid-demo.ornl.gov/sans-backend
3. Run the script with `sudo bash sans-backend-run.sh -h` to see the help menu.

Current options include:
* -i) launches a bash shell
* -u) forces an update of the application.
* -h) prints this message.

You must download the wrapper script from the above link as the build process modifies the copy in version control.

Set-up for development in a virtual environment
-----------------------------------------------

This is the instructions for a person who is developing *both*
drt-sans and mantid in tandem. It assumes that you are using a virtual
environment and have a local build of mantid. As one caveat, for this
method, mantid must be build against the same version of python being
used in the virtual environment

1. Checkout the code.

```sh
$ git clone git@code.ornl.gov:sns-hfir-scse/sans/sans-backend.git
$ cd sans-backend
```

2. Create the virtual environment and activate it. The
   `--system-site-packages` argument lets it use things installed on
   the system as a whole for compiling mantid
```sh
$ virtualenv -p $(which python3) --system-site-packages .venv
$ source .venv/bin/activate
```

3. Add the built version of mantid to the python path in the virtual
   environment
```sh
$ ~/build/mantid/bin/AddPythonPath.py
```

4. Install the code in `develop` mode.
```sh
$ python setup.py develop
```

5. Try it out. Start `python` and try
```python
import mantid
import ornl
```
Verify you can run the unit tests:
```sh
python -m pytest tests/unit/new/
```

6. Be done. Deactivate the virtual environment using
```
$ deactivate
```

As an alternative, you can use [direnv](https://direnv.net) to do a
fair amount of the work. Unfortunately, it doesn't currently handle
`--system-site-packages` correctly so you'll have to work around
it. Create the virtual environment using
```sh
$ virtualenv -p $(which python3) --system-site-packages .direnv/python-$(python3 -c "import platform as p;print(p.python_version())")
```
Then you create a file `.envrc` in your source directory
```sh
$ echo "layout python3" > .envrc
```
Finally, tell direnv that you want it to work in this directory
```sh
$ direnv allow
```
Then follow steps 3-5 from above. After this, the virtual environment
with load when you enter the source tree, and unload when you leave.


# NEXT branch

Branched off of master on June 18, 2019 in preparation for version 1.0
