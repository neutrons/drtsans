FROM ubuntu:latest

WORKDIR /opt/sns-backend
COPY . /opt/sns-backend

RUN apt update && apt dist-upgrade -y
RUN apt install -y wget python3-pip python3-dev freeglut3-dev libglu1-mesa 
RUN wget https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-x86_64.sh
RUN chmod +x Anaconda3-2019.03-Linux-x86_64.sh
RUN ./Anaconda3-2019.03-Linux-x86_64.sh -b -p /opt/anaconda
RUN /opt/anaconda/bin/conda init bash
RUN ./conda-setup.sh

#RUN pip3 install --upgrade pip
#RUN pip3 install virtualenv


#ENV /opt/anaconda/lib/python3.6/site-packages/
#RUN virtualenv -p $(which python3) venv && . venv/bin/activate && pip3 install -r requirements_test.txt

