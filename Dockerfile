FROM code.ornl.gov:4567/sns-hfir-scse/docker-containers/mantid-framework-nightly/master

COPY . /opt/sans-backend

ENV PYTHONPATH $PYTHONPATH:/opt/conda/lib/python3.6/site-packages:/opt/conda/envs/mantid/bin/:/opt/conda/bin:/opt/sans-backend
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/opt/conda/lib/mantid/plugins/