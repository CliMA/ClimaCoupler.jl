#!/bin/bash
# Hang watchdog for calibration_relay.sh. It handles exactly ONE pathology:
# the driver process is alive but its workers died (walltime), so it will
# hang on ClimaCalibrate's empty_pool_timeout for up to 11 hours with no
# exit for the relay to observe. Everything else - clean exits, crashes,
# NaN'd members - is the relay's job, so this script only ever KILLS the
# driver's process group; the relay's normal loop then relaunches.
#
#   calibration_watchdog.sh <run_dir>
#
# Started (and killed) automatically by calibration_relay.sh; also runnable
# standalone against a manually launched driver if <run_dir>/relay_driver.pid
# holds the driver's process-group id.
#
# Kill condition: THREE consecutive wakes (default 30 min apart) where all of
#   1. the driver process group is alive (pid file valid),
#   2. no file under run_dir has been written in WATCHDOG_SILENCE_MIN
#      (default 45; relay bookkeeping files excluded, since the watchdog's
#      own logging must not reset the silence clock),
#   3. qstat succeeds and shows no matching worker jobs (default pattern
#      "julia", ClimaCalibrate's worker job prefix) - a driver waiting on a
#      QUEUED worker is healthy, and a failed qstat is inconclusive (login
#      qstat is a laggy cache; the single-miss false positive is a measured
#      failure mode).
# Any contrary evidence resets the strike counter. Three strikes at 30 min
# on top of 45 min silence detects a hang in ~2 h - slow, but the
# alternative cost is an 11 h empty_pool_timeout stall, and dead workers
# bill nothing while we wait.
set -u

[ $# -eq 1 ] || { echo "usage: calibration_watchdog.sh <run_dir>" >&2; exit 64; }
RUN_DIR=$1
PERIOD_MIN=${WATCHDOG_PERIOD_MIN:-30}
SILENCE_MIN=${WATCHDOG_SILENCE_MIN:-45}
JOB_PATTERN=${WATCHDOG_JOB_PATTERN:-julia}
RELAY_LOG=$RUN_DIR/relay.log
PID_FILE=$RUN_DIR/relay_driver.pid

log() { echo "=== $(date '+%F %T') watchdog: $*" | tee -a "$RELAY_LOG"; }

log "start: period ${PERIOD_MIN}m, silence threshold ${SILENCE_MIN}m, job pattern '$JOB_PATTERN'"

strikes=0
while :; do
    sleep $((PERIOD_MIN * 60))

    pgid=$(cat "$PID_FILE" 2>/dev/null) || { strikes=0; continue; }
    case $pgid in ''|*[!0-9]*) strikes=0; continue;; esac
    # Driver group already gone -> the relay is handling it.
    kill -0 -- "-$pgid" 2>/dev/null || { strikes=0; continue; }

    # Any recent write anywhere in the run dir counts as life. Exclude the
    # relay's own bookkeeping so logging does not defeat the check.
    recent=$(find "$RUN_DIR" -type f \
                  ! -name relay.log ! -name relay_driver.pid \
                  -mmin "-$SILENCE_MIN" -print -quit 2>/dev/null)
    [ -n "$recent" ] && { strikes=0; continue; }

    # qstat must SUCCEED and show no workers; a failed qstat is inconclusive.
    if ! qout=$(qstat -u "$USER" 2>/dev/null); then
        log "qstat unavailable; inconclusive, no strike"
        continue
    fi
    if echo "$qout" | grep -q "$JOB_PATTERN"; then
        strikes=0
        continue
    fi

    strikes=$((strikes + 1))
    log "driver alive (pgid $pgid), run dir silent >${SILENCE_MIN}m, no '$JOB_PATTERN' jobs - strike $strikes/3"
    if [ "$strikes" -ge 3 ]; then
        log "KILLING hung driver process group $pgid"
        kill -TERM -- "-$pgid" 2>/dev/null
        sleep 60
        kill -KILL -- "-$pgid" 2>/dev/null
        strikes=0
    fi
done
