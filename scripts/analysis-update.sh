cat /etc/idmapd.conf

export PATH=/SNS/software/miniconda2/bin:$PATH
source activate sans
conda install -q -y  -c mantid/label/nightly mantid-framework
pip install https://$CI_USER:$CI_PASS@code.ornl.gov/sns-hfir-scse/sans/sans-backend.git