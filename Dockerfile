FROM continuumio/miniconda3

RUN apt update
RUN apt install -y freeglut3-dev libglu1-mesa
RUN conda create -n mantid python=3.6
RUN echo "source activate mantid" > ~/.bashrc
ENV PATH /opt/conda/envs/mantid/bin:$PATH

COPY . /opt/sns-backend
WORKDIR /opt/sns-backend
RUN .scripts/conda-setup.sh

ENV PYTHONPATH $PYTHONPATH:/opt/conda/lib/python3.6/site-packages:/opt/conda/envs/mantid/bin/:/opt/conda/bin:/opt/sans-backend
