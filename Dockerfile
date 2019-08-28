FROM code.ornl.gov:4567/sns-hfir-scse/docker-containers/mantid-framework-nightly/master

WORKDIR /tmp/input
USER snsdata

COPY ornl tests scripts .gitattributes MANIFEST.in pytest.ini setup.cfg setup.py test_job.sh versioneer.py /opt/sans-backend/

RUN ls /opt/sans-backend/scripts/
ENV PYTHONPATH $PYTHONPATH:/opt/conda/lib/python3.6/site-packages:/opt/conda/envs/mantid/bin/:/opt/conda/bin:/opt/sans-backend
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/opt/conda/lib/mantid/plugins/