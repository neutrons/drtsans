FROM code.ornl.gov:4567/sns-hfir-scse/docker-containers/mantid-framework-nightly/master

WORKDIR /opt/sns-backend
COPY . /opt/sns-backend
RUN ./scripts/conda-setup.sh

ENV PYTHONPATH $PYTHONPATH:/opt/conda/lib/python3.6/site-packages:/opt/conda/envs/mantid/bin/:/opt/conda/bin:/opt/sans-backend
