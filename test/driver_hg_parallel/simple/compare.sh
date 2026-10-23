#!/usr/bin/env bash
# One graph, eight MPI ranks, one run per reader; compare final logged CUTL (km1).
set -euo pipefail
if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: bash $0 /path/to/zdrive.exe /path/to/graph.hg [new-output-directory]" >&2
  exit 2
fi
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
zdrive=$(realpath -- "$1")
graph=$(realpath -- "$2")
output=${3:-reader-comparison}
[[ -x "$zdrive" && -f "$graph" ]] || { echo "Missing binary or graph" >&2; exit 2; }
mkdir -- "$output"
cp -- "$here/old.inp" "$here/parallel.inp" "$output/"
ln -s -- "$graph" "$output/graph.hg"
cd -- "$output"
for reader in old parallel; do
  echo "Running $reader reader on 8 ranks; log: $PWD/$reader.log"
  "${MPIEXEC:-mpirun}" -np 8 "$zdrive" "$reader.inp" > "$reader.log" 2>&1
done
# Confirm the new reader ran, then select CUTL only AFTER load balancing.
grep -qF 'Read hypergraph in parallel using MPI-IO' parallel.log
objective() {
  sed -n '/^AFTER load balancing$/,$p' "$1" |
    grep -F 'Zoltan_LB_Eval_HG  CUTL (' | tail -n 1 | awk '{print $NF}'
}
old=$(objective old.log) || { echo 'FAIL: final objective missing from old.log' >&2; exit 1; }
parallel=$(objective parallel.log) || { echo 'FAIL: final objective missing from parallel.log' >&2; exit 1; }
printf 'Old CUTL (km1): %s\nNew CUTL (km1): %s\n' "$old" "$parallel"
if [[ -n "$old" && "$old" == "$parallel" ]]; then
  echo 'PASS: final printed objectives are equal'
else
  echo 'FAIL: final printed objectives differ or are missing' >&2
  exit 1
fi
