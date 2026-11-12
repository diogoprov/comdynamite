#!/bin/bash
# Launch one spot instance that runs a simulation script of data-raw/ (e.g.
# data-raw/s1_type1_shifts.R) and terminates at the end (adapted from the AWS recipe of
# crywolf/dev/dyn-sim/aws and the compendium fpcalib).
# Usage: data-raw/aws/launch_spot.sh <bucket> <run> <script> <instance_type> <workers> <instance_profile> [max_hours] [reps]
#   e.g. data-raw/aws/launch_spot.sh my-bucket dyn-test data-raw/s1_type1_shifts.R c7a.2xlarge 8 pcps-sim-ec2 2 1
# <run>: label of the run; code, job files, results and logs go to
#   s3://<bucket>/fpcalib/runs/<run>/ (upload the code first with pack_code.sh;
#   the instance profile pcps-sim-ec2 writes only under fpcalib/, hence the prefix).
#   Launching again with the same <run> resumes from the job files saved there.
# max_hours: hard limit after which the instance terminates (default 12).
# reps: sets REPS for the script (default: the script's own default).
# <instance_profile>: IAM instance profile with read/write access to the prefix.
# JOBS (environment, optional): job folder of the script (default
#   data-raw/s1_jobs; other scripts set their own).
# EXTRA_ENV (environment, optional): more variables for the script, e.g. "PARTITION=0"
#   (space-separated NAME=value pairs).
# AZ (environment, optional): availability zone, e.g. us-east-1c, to pick the
#   cheapest spot price; otherwise EC2 chooses.
# MARKET (environment, optional): "ondemand" for an on-demand instance, which is not
#   interrupted (for jobs longer than the usual time between spot interruptions); default spot.
set -euo pipefail
BUCKET=$1; RUN=$2; SCRIPT=$3; TYPE=$4; WORKERS=$5; ROLE=$6; MAX_HOURS=${7:-12}; REPS=${8:-}
cd "$(git rev-parse --show-toplevel)"
JOBS=${JOBS:-data-raw/s1_jobs}
aws s3 ls "s3://${BUCKET}/fpcalib/runs/${RUN}/code.tar.gz" >/dev/null || \
  { echo "No code for run ${RUN}; run data-raw/aws/pack_code.sh ${BUCKET} ${RUN} first." >&2; exit 1; }
UD=$(mktemp)
sed -e "s|__BUCKET__|${BUCKET}|" -e "s|__RUN__|${RUN}|" -e "s|__SCRIPT__|${SCRIPT}|" \
  -e "s|__JOBS__|${JOBS}|" -e "s|__WORKERS__|${WORKERS}|" -e "s|__REPS__|${REPS}|" \
  -e "s|__MAX_HOURS__|${MAX_HOURS}|" -e "s|__EXTRA_ENV__|${EXTRA_ENV:-}|" \
  data-raw/aws/user_data.sh > "${UD}"

# latest Amazon Linux 2023 AMI (x86_64), looked up without SSM permissions
AMI=$(aws ec2 describe-images --owners amazon \
  --filters "Name=name,Values=al2023-ami-2023.*-kernel-*-x86_64" "Name=state,Values=available" \
  --query 'sort_by(Images, &CreationDate)[-1].ImageId' --output text)
echo "AMI ${AMI}; job folder ${JOBS}; market ${MARKET:-spot}" >&2
SPOT='{"MarketType":"spot","SpotOptions":{"SpotInstanceType":"one-time","InstanceInterruptionBehavior":"terminate"}}'
if [ "${MARKET:-spot}" = "ondemand" ]; then SPOT=""; fi

aws ec2 run-instances \
  --image-id "${AMI}" \
  --instance-type "${TYPE}" \
  ${AZ:+--placement AvailabilityZone=${AZ}} \
  ${SPOT:+--instance-market-options "${SPOT}"} \
  --instance-initiated-shutdown-behavior terminate \
  --iam-instance-profile Name="${ROLE}" \
  --block-device-mappings '[{"DeviceName":"/dev/xvda","Ebs":{"VolumeSize":40,"VolumeType":"gp3","DeleteOnTermination":true}}]' \
  --user-data "file://${UD}" \
  --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=comdynamite-${RUN}},{Key=project,Value=comdynamite}]" \
                       "ResourceType=volume,Tags=[{Key=Name,Value=comdynamite-${RUN}},{Key=project,Value=comdynamite}]" \
  --query 'Instances[0].InstanceId' --output text
rm -f "${UD}"
