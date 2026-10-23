#!/usr/bin/env bash
# Compare the final logged connectivity objective for the two readers.
set -euo pipefail

# Settings: edit these defaults or supply them as environment variables.
ZDRIVE=${ZDRIVE:-./build/src/driver/zdrive.exe}
HYPERGRAPH=${HYPERGRAPH:-./webbase-1M.mtx.zoltan.hg}
CORES=${CORES:-8}       # One MPI rank per core.
BLOCKS=${BLOCKS:-2}     # Number of partition blocks, independent of core count.
OUTPUT_DIR=${OUTPUT_DIR:-reader-comparison}

fail() {
  printf 'FAIL: %s\n' "$1" >&2
  exit 1
}

extract_objective() {
  # Ignore the initial evaluation; select CUTL only after load balancing.
  sed -n '/^AFTER load balancing$/,$p' "$1" |
    grep -F 'Zoltan_LB_Eval_HG  CUTL (' |
    tail -n 1 |
    awk '{print $NF}'
}

if [[ $# -eq 1 || $# -gt 3 ]]; then
  printf 'Usage: bash %s [ZDRIVE GRAPH.hg [new-output-directory]]\n' "$0" >&2
  exit 2
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
zdrive=$(realpath -- "${1:-$ZDRIVE}")
graph=$(realpath -- "${2:-$HYPERGRAPH}")
output=${3:-$OUTPUT_DIR}

[[ "$CORES" =~ ^[1-9][0-9]*$ ]] || fail 'CORES must be a positive integer'
[[ "$BLOCKS" =~ ^[1-9][0-9]*$ ]] || fail 'BLOCKS must be a positive integer'
[[ -x "$zdrive" ]] || fail "Binary is not executable: $zdrive"
[[ -f "$graph" ]] || fail "Graph does not exist: $graph"

# Keep the exact inputs and logs together in a fresh directory.
mkdir -- "$output"
for reader in old parallel; do
  sed -E "s/NUM_GLOBAL_PARTS=.*/NUM_GLOBAL_PARTS=$BLOCKS/" \
    "$script_dir/$reader.inp" > "$output/$reader.inp"
done
ln -s -- "$graph" "$output/graph.hg"
cd -- "$output"

printf '\nReader comparison\n  Graph:  %s\n  Cores:  %s\n  Blocks: %s\n  Logs:   %s\n\n' \
  "$graph" "$CORES" "$BLOCKS" "$PWD"

for reader in old parallel; do
  printf 'Running %-8s reader ...\n' "$reader"
  "${MPIEXEC:-mpirun}" -np "$CORES" "$zdrive" "$reader.inp" \
    > "$reader.log" 2>&1 || fail "$reader run failed; see $PWD/$reader.log"
done

grep -qF 'Read hypergraph in parallel using MPI-IO' parallel.log ||
  fail 'Parallel reader was not activated; see parallel.log'

old_objective=$(extract_objective old.log) ||
  fail 'Final objective missing from old.log'
parallel_objective=$(extract_objective parallel.log) ||
  fail 'Final objective missing from parallel.log'

printf '\n%-12s %s\n' 'Reader' 'CUTL (km1)'
printf '%-12s %s\n' '------------' '----------'
printf '%-12s %s\n' 'Old' "$old_objective" 'Parallel' "$parallel_objective"

[[ -n "$old_objective" && "$old_objective" == "$parallel_objective" ]] ||
  fail 'Final printed objectives differ or are missing'

printf '\nPASS: final printed objectives are equal\n'
