#!/bin/zsh
# Overnight batch: v4 dynamic occupancy for the NAAMP species with longevity data.
# Run from anywhere; works from the repository root. Waits until START (HH:MM; empty = now),
# runs 4 species at a time and skips species whose fit already exists (rerun to resume).
# Launch e.g.: nohup caffeinate -i -s data-raw/persistence_proxy/run_overnight.sh &
cd "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
P=data-raw/persistence_proxy; L=$P/fits/logs; mkdir -p $L
START=${START:-}
if [[ -n "$START" ]]; then
  echo "scheduled $(date) for $START" >> $L/overnight.log
  while [[ "$(date +%H:%M)" != "$START" ]]; do sleep 20; done
fi
echo "start $(date)" >> $L/overnight.log
Rscript -e 'x <- readRDS("data-raw/persistence_proxy/species_table.rds"); x <- x$Species[!is.na(x$lon)];
  first <- c("Anaxyrus americanus", "Gastrophryne carolinensis"); cat(c(first, setdiff(x, first)), sep = "\n")' > $L/queue.txt
P=$P L=$L xargs -P 4 -I{} zsh -c '
  sp="$1"; f=$P/fits/pilot4_${sp// /_}.rds
  if [[ -f $f ]]; then echo "skip $sp" >> $L/overnight.log; exit 0; fi
  it=1000; [[ "$sp" == "Anaxyrus americanus" || "$sp" == "Gastrophryne carolinensis" ]] && it=2000
  SPECIES="$sp" ITER=$it Rscript $P/pilot4.R > $L/v4_${sp// /_}.log 2>&1
  echo "done $sp (iter $it, exit $?) $(date)" >> $L/overnight.log' _ {} < $L/queue.txt
echo "finished $(date)" >> $L/overnight.log
