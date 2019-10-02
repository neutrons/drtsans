FROM code.ornl.gov:4567/sns-hfir-scse/docker-containers/mantid-framework-nightly/master

RUN apt-get install -y sssd
COPY sssd.conf /etc/sssd/sssd.conf
RUN chmod 600 /etc/sssd/sssd.conf
RUN chown root:root /etc/sssd/sssd.conf
COPY nsswitch.conf /etc/nsswitch.conf
COPY SNSCA.cert.pem /etc/pki/tls/certs/

WORKDIR /tmp/input
USER snsdata

COPY drtsans /opt/sans-backend/drtsans
COPY tests /opt/sans-backend/tests
COPY scripts /opt/sans-backend/scripts
COPY docs /opt/sans-backend/docs
COPY .gitattributes MANIFEST.in pytest.ini setup.cfg setup.py test_job.sh versioneer.py requirements.txt requirements_dev.txt /opt/sans-backend/

ENV PYTHONPATH $PYTHONPATH:/opt/conda/lib/python3.6/site-packages:/opt/conda/envs/mantid/bin/:/opt/conda/bin:/opt/sans-backend
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/opt/conda/lib/mantid/plugins/

RUN bash /opt/sans-backend/scripts/conda-setup.sh