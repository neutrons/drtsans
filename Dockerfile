FROM code.ornl.gov:4567/sns-hfir-scse/docker-containers/mantid-framework-nightly/master

WORKDIR /opt/sns-backend
COPY ornl /opt/sns-backend/ornl
COPY scripts/conda-setup.sh /opt/sns-backend/scripts/conda-setup.sh
COPY scripts/run-inputs.sh /opt/sns-backend/scripts/run-inputs.sh
RUN ./scripts/conda-setup.sh

ENV PYTHONPATH $PYTHONPATH:/opt/conda/lib/python3.6/site-packages:/opt/conda/envs/mantid/bin/:/opt/conda/bin:/opt/sans-backend
