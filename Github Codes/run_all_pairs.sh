#!/usr/bin/env bash
#
# run_all_pairs.sh
#
# Sweep the MR-CCC primary analysis over all 20 ordered cell-type pairs,
# one fresh R process per pair.
#
# Usage, FROM THE PROJECT ROOT (not from Github Codes/):
#
#   bash "Github Codes/run_all_pairs.sh" 2>&1 | tee Results/run_all_pairs.log
#
# Expect roughly 60-120 minutes per pair (100,000 iterations x 4 chains, with
# pathway scoring dominating), so about twenty-six hours in total. The log file
# lets you check progress without disturbing the run.
#
# RESUMABLE. run_pair.R skips any pair whose Results/<C1>_<C2>_MR_CCC.rds
# already exists, so if the sweep is interrupted, re-running this script
# picks up where it stopped. Delete a pair's .rds to force a re-run.
#
# The pairs are ordered so that the main-text axis runs first and the
# remaining nineteen follow; if the sweep has to be cut short, the pairs
# that matter most are already done.

set -u

CELLS=(NKCells MonocytesCells BCells CD4Cells CD8Cells)

# Main-text axis first.
ORDER=("NKCells MonocytesCells")
for c1 in "${CELLS[@]}"; do
  for c2 in "${CELLS[@]}"; do
    [ "$c1" = "$c2" ] && continue
    [ "$c1 $c2" = "NKCells MonocytesCells" ] && continue
    ORDER+=("$c1 $c2")
  done
done

echo "Sweeping ${#ORDER[@]} ordered pairs."
START=$(date +%s)

for pair in "${ORDER[@]}"; do
  # shellcheck disable=SC2086
  Rscript "Github Codes/run_pair.R" $pair
  status=$?
  if [ $status -ne 0 ]; then
    # Do not abort the whole sweep for one bad pair: record it and continue,
    # so that a single failure late at night does not cost the remaining
    # pairs. The .rds for the failed pair simply will not exist, and
    # re-running this script will retry it.
    echo "WARNING: pair '$pair' exited with status $status -- continuing."
  fi
done

ELAPSED=$(( ($(date +%s) - START) / 60 ))
echo "Sweep finished in ${ELAPSED} minutes."
echo "Pairs with output:"
ls -1 Results/*_MR_CCC.rds 2>/dev/null | wc -l
