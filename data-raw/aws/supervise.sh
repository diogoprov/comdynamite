#!/bin/bash
# (Adapted from crywolf/dev/dyn-sim/aws and fpcalib/aws.) Keep a spot run alive: until s3://<bucket>/fpcalib/runs/<run>/ has a DONE or
# FAILED marker, relaunch the run (resuming from its job files) whenever its
# instance ends, trying the instance types given in order (give <workers> as "auto" to use all
# vCPUs but two of whichever type starts). When no type has spot capacity, wait 5 min and try
# again; that does not count as a relaunch. Prints one line per event.
# Usage: data-raw/aws/supervise.sh <bucket> <run> <script> <workers> <instance_profile> <max_hours> <instance_id> <type> [<type> ...]
#   e.g. data-raw/aws/supervise.sh my-bucket dyn data-raw/s1_type1_shifts.R 190 pcps-sim-ec2 12 i-0abc c7a.48xlarge c7i.48xlarge m7a.48xlarge
# MAXR (environment, default 5): number of relaunches before giving up.
# MAXWAIT (environment, default 24): attempts in a row without capacity before giving up.
set -uo pipefail
BUCKET=$1; RUN=$2; SCRIPT=$3; WORKERS=$4; ROLE=$5; HOURS=$6; ID=$7; shift 7; TYPES=("$@")
MAXR=${MAXR:-5}
MAXWAIT=${MAXWAIT:-24}
wait_n=0
cd "$(git rev-parse --show-toplevel)"
S3=s3://${BUCKET}/fpcalib/runs/${RUN}
marker() { aws s3 ls "${S3}/$1" >/dev/null 2>&1; }
state() { aws ec2 describe-instances --instance-ids "$1" \
  --query 'Reservations[].Instances[].State.Name' --output text 2>/dev/null || echo gone; }
n=0
while true; do
  if marker DONE; then echo "${RUN}: finished ($(date))"; exit 0; fi
  if marker FAILED; then echo "${RUN}: script failed, see ${S3}/logs/ ($(date))"; exit 1; fi
  s=$(state "$ID")
  case "$s" in
    running|pending) sleep 300; continue;;
  esac
  [ "$n" -ge "$MAXR" ] && { echo "${RUN}: gave up after ${n} relaunches ($(date))"; exit 1; }
  launched=""
  for T in "${TYPES[@]}"; do
    out=$(data-raw/aws/launch_spot.sh "$BUCKET" "$RUN" "$SCRIPT" "$T" "$WORKERS" "$ROLE" "$HOURS" 2>&1 | tail -1)
    case "$out" in i-*) launched=$T; ID=$out; break;; esac
  done
  if [ -n "$launched" ]; then
    n=$((n + 1)); wait_n=0
    echo "${RUN}: instance ended (${s}); relaunch ${n} on ${launched} ${ID} ($(date))"
  else
    wait_n=$((wait_n + 1))
    echo "${RUN}: no spot capacity for ${TYPES[*]} (${wait_n}/${MAXWAIT}): ${out} ($(date))"
    [ "$wait_n" -ge "$MAXWAIT" ] && { echo "${RUN}: gave up waiting for capacity ($(date))"; exit 1; }
  fi
  sleep 300
done
