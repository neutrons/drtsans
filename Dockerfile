FROM code.ornl.gov:4567/sns-hfir-scse/docker-containers/mantid-framework-nightly/master

WORKDIR /tmp/input

COPY ornl /opt/sans-backend/ornl
COPY tests /opt/sans-backend/tests
COPY scripts /opt/sans-backend/scripts
COPY .gitattributes /opt/sans-backend/.gitattributes
COPY MANIFEST.in /opt/sans-backend/MANIFEST.in
COPY pytest.ini /opt/sans-backend/pytest.ini
COPY requirements.txt /opt/sans-backend/requirements.txt
COPY setup.cfg /opt/sans-backend/setup.cfg
COPY setup.py /opt/sans-backend/setup.py
COPY test_job.sh /opt/sans-backend/test_job.sh
COPY versioneer.py /opt/sans-backend/versioneer.py

ENV PYTHONPATH $PYTHONPATH:/opt/conda/lib/python3.6/site-packages:/opt/conda/envs/mantid/bin/:/opt/conda/bin:/opt/sans-backend
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/opt/conda/lib/mantid/plugins/