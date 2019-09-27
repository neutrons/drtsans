groupadd -g 52321 sns_ndav_team
useradd -s /bin/bash -u 13972 -m -g sns_ndav_team lj7
runuser -l lj7 -c "bash /opt/sans-backend/scripts/analysis-update.sh"