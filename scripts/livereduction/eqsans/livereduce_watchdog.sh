#!/usr/bin/env bash

# --- CONFIG ---
WATCHDOG_TARGET="/var/log/SNS_applications/livereduce.log"
SERVICE="livereduce.service"
INTERVAL=60         # how often we check the file, in seconds
THRESHOLD=300       # inactivity threshold in seconds

WATCHDOG_LOG="/var/log/SNS_applications/livereduce_watchdog.log"
# Track when we last issued a restart so we don't keep restarting every loop
last_restart=0

# Infinite loop
while true; do

  if [[ ! -e "$WATCHDOG_TARGET" ]]; then
    echo "$(date '+%F %T') ERROR: '$WATCHDOG_TARGET' not found." >> "$WATCHDOG_LOG"

  else
    # Get file mtime (epoch seconds) and current time
    last_mod=$(stat -c %Y "$WATCHDOG_TARGET")
    now=$(date +%s)
    age=$(( now - last_mod ))

    if (( age >= THRESHOLD )); then
      # Only restart if we haven't already done so in this inactivity window
      since_restart=$(( now - last_restart ))
      if (( since_restart >= THRESHOLD )); then
        echo -e "\n#############################################################################" >> "$WATCHDOG_LOG"
        echo "$(date '+%F %T') No change for $age s in $WATCHDOG_TARGET" >> "$WATCHDOG_LOG"
        echo "---- Last 20 lines of $WATCHDOG_TARGET before restart:" >> "$WATCHDOG_LOG"
        tail -n 20 "$WATCHDOG_TARGET" >> "$WATCHDOG_LOG"
        echo -e "\nrestarting $SERVICE." >> "$WATCHDOG_LOG"
        # Restart the service (use systemctl or service as appropriate)
        if command -v systemctl &>/dev/null; then
          systemctl stop "$SERVICE"
          systemctl start "$SERVICE"
          systemctl status "$SERVICE" >> "$WATCHDOG_LOG"
        else
          service "$SERVICE" stop
          service "$SERVICE" start >> "$WATCHDOG_LOG"
          service "$SERVICE" status >> "$WATCHDOG_LOG"
        fi
        last_restart=$now
      fi
    fi
  fi

  sleep "$INTERVAL"
done
