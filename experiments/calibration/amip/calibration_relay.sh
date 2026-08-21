#!/bin/bash
# Relay wrapper for long calibrations: run a driver command repeatedly until
# the target number of iterations exists on disk, relying on ClimaCalibrate's
# resume (completed iterations are skipped on relaunch). Generic over
# calibrations - it knows nothing about configs; it only watches the run
# directory and the driver's exit.
#
#   calibration_relay.sh <run_dir> <target_iterations> <driver command...>
#
# Example (inside tmux on a login node):
#
#   bash experiments/calibration/amip/calibration_relay.sh \
#       /glade/derecho/scratch/$USER/my_calibration_out 10 \
#       bash /glade/derecho/scratch/$USER/my_calibration_out/driver_my_run.sh
#
# The config's n_iterations must be set to the SAME total target; the relay
# never edits configs, it just relaunches until iteration_<target> is done.
#
# Progress is measured as the number of iteration_*/G_ensemble.jld2 files: a
# completed iteration always has one, and the trailing posterior-only
# directory (iteration_<N+1>) never does.
#
# Failure policy:
#   - An attempt that completes ZERO new iterations is a "fast failure"
#     (config error, model NaN'ing every member, broken environment).
#     RELAY_FASTFAIL_LIMIT consecutive fast failures (default 2) stop the
#     relay rather than rerunning the same failure forever.
#   - RELAY_MAX_RELAUNCHES (default target/2 + 2, from the ~2-8 iterations
#     that fit one 12 h worker walltime) bounds the total attempts.
#
# A hang watchdog (calibration_watchdog.sh, same directory) is started
# automatically and killed on exit; it handles the one case the relay cannot
# see - a driver alive but stalled on a dead worker pool - by killing the
# driver's process group so the relay's normal relaunch path takes over.
# Disable it with RELAY_NO_WATCHDOG=1.
#
# Every decision is appended to <run_dir>/relay.log.
set -u

[ $# -ge 3 ] || {
    echo "usage: calibration_relay.sh <run_dir> <target_iterations> <driver command...>" >&2
    exit 64
}
RUN_DIR=$1
TARGET=$2
shift 2

[ -d "$RUN_DIR" ] || { echo "run_dir $RUN_DIR does not exist" >&2; exit 64; }
case $TARGET in ''|*[!0-9]*) echo "target_iterations must be a number" >&2; exit 64;; esac

RELAY_LOG=$RUN_DIR/relay.log
PID_FILE=$RUN_DIR/relay_driver.pid
FASTFAIL_LIMIT=${RELAY_FASTFAIL_LIMIT:-2}
MAX_RELAUNCHES=${RELAY_MAX_RELAUNCHES:-$((TARGET / 2 + 2))}

log() { echo "=== $(date '+%F %T') relay: $*" | tee -a "$RELAY_LOG"; }

count_done() { find "$RUN_DIR" -maxdepth 2 -name G_ensemble.jld2 2>/dev/null | wc -l; }

# --- watchdog -----------------------------------------------------------
WATCHDOG_PID=""
if [ "${RELAY_NO_WATCHDOG:-0}" != "1" ]; then
    here=$(cd "$(dirname "$0")" && pwd)
    bash "$here/calibration_watchdog.sh" "$RUN_DIR" &
    WATCHDOG_PID=$!
    log "watchdog started (pid $WATCHDOG_PID)"
fi
cleanup() {
    [ -n "$WATCHDOG_PID" ] && kill "$WATCHDOG_PID" 2>/dev/null
    rm -f "$PID_FILE"
}
trap cleanup EXIT

# --- relay loop ----------------------------------------------------------
log "start: target $TARGET iterations in $RUN_DIR; driver: $*"
log "limits: max $MAX_RELAUNCHES attempts, stop after $FASTFAIL_LIMIT zero-progress attempts"

attempt=0
fastfail=0
while :; do
    done_before=$(count_done)
    if [ "$done_before" -ge "$TARGET" ]; then
        log "COMPLETE: $done_before/$TARGET iterations on disk"
        exit 0
    fi
    attempt=$((attempt + 1))
    if [ "$attempt" -gt "$MAX_RELAUNCHES" ]; then
        log "STOP: relaunch budget exhausted ($MAX_RELAUNCHES attempts, $done_before/$TARGET done)"
        exit 1
    fi

    log "attempt $attempt/$MAX_RELAUNCHES: $done_before/$TARGET done, launching driver"
    # Byte offset of driver.log before the attempt, for the worker-error scan
    # below (driver.log is appended across attempts).
    dlog=$RUN_DIR/driver.log
    dlog_off=$(stat -c %s "$dlog" 2>/dev/null || echo 0)
    # New session/process group so the watchdog (and we) can kill driver +
    # children together. setsid forks iff the caller gave the background job
    # its own group (job-control shells) and execs in place otherwise, so
    # $! is NOT a reliable pgid; instead the driver process writes its own
    # pid (== the new group's pgid under either behavior) to the pid file,
    # and --wait keeps $! alive to propagate the driver's exit status.
    setsid --wait bash -c 'echo "$$" > "$1"; shift; exec "$@"' _ "$PID_FILE" "$@" &
    wait $!
    rc=$?
    rm -f "$PID_FILE"

    done_after=$(count_done)
    log "attempt $attempt: driver exited rc=$rc, iterations $done_before -> $done_after"

    # Case-B advisory: workers dying LATE in an iteration (< 50% of members
    # still in flight) do not trip ClimaCalibrate's failure-rate halt, so the
    # iteration completes with imputed members and stands on disk as normal
    # progress. Detect the signature - iterations advanced while this
    # attempt's slice of driver.log grew worker-error lines - and flag it for
    # the post-run review; the relay cannot (and should not) undo it.
    if [ "$done_after" -gt "$done_before" ] && [ -f "$dlog" ]; then
        nerr=$(tail -c +"$((dlog_off + 1))" "$dlog" 2>/dev/null |
               grep -c "${RELAY_WORKER_ERROR_PATTERN:-Error running on worker}")
        [ "$nerr" -gt 0 ] && log "ADVISORY: attempt $attempt advanced iterations" \
            "but logged $nerr worker error(s) - an iteration completed in this" \
            "attempt may contain imputed (failed) members; verify its" \
            "G_ensemble before trusting the update"
    fi

    if [ "$done_after" -le "$done_before" ]; then
        fastfail=$((fastfail + 1))
        log "zero progress ($fastfail/$FASTFAIL_LIMIT consecutive)"
        if [ "$fastfail" -ge "$FASTFAIL_LIMIT" ]; then
            log "STOP: $fastfail consecutive zero-progress attempts - not a walltime" \
                "pattern; inspect driver.log before relaunching"
            exit 2
        fi
    else
        fastfail=0
    fi
done
