# Parallel native hypergraph input

## Small suite you can inspect and run

After building the MPI driver, run this from the source root. `--build` is the
configured **Autotools build directory**, using the same MPI installation as
`mpicc` and `mpirun`:

```sh
python3 test/driver_hg_parallel/check_reader.py \
  --build /path/to/build --ranks 1,2,4 \
  --output /tmp/zoltan-reader-check
```

Use a new output directory each time. The script compiles the small C comparison
probe automatically and keeps all files for inspection. On a cluster, run inside
an allocation with enough cores. `--mpicc`, `--mpiexec`, and `--mpi-args` let you
select the matching MPI compiler and launcher. You can use `--ranks 1,8` to
include ranks with no initial vertices or edges. No Python packages are required.

The fixtures contain six vertices and four edges:

- `fixtures/weighted-zero.hg`: zero-based PaToH input with vertex/edge weights.
- `fixtures/unweighted-one.hg`: one-based IBM input without weights.
- `fixtures/expected.json`: supplied assignments with hand-computed objectives
  and balance, including a three-part example where cut-net and km1 differ.

The checks are:

1. Validate the Python objective evaluator against those hand-worked examples.
2. Compare each rank's old/new reader arrays byte for byte under linear,
   cyclic, and owner distribution: offsets, IDs, ordered pins, owners and weights.
3. Run both readers with identical PHG parameters and seed for both cut objectives.
   Require identical final assignments; independently calculate cut-net, km1,
   partition weights and imbalance from the original text file. Also compare
   the calculated objectives with the driver's final evaluation.

`expected.json` describes the supplied example assignments; actual PHG output is
evaluated separately. Tiny weighted cases use an imbalance tolerance of 1.5 to
accommodate coarse vertex weights. Both readers use exactly the same settings.

For `--ranks 1,2,4`, successful output ends with:

```text
PASS: 18 exact input comparisons and 12 final partition comparisons. Results: ...
```

Inspect `summary.json`, each case's `zdrive.inp`, `run.log`, `assignments.json`,
`metrics.json`, and the native `graph.out.*` files in the output directory.
Probe logs end with `Exact per-rank reader equality: PASS`. Any mismatch or
command failure causes a nonzero exit status. The larger suite described below
also checks malformed files, empty edges, long records and vector weights.

## Reader behavior

The C driver (`zdrive`) can load an **uncompressed native `.hg` file** with
parallel MPI-IO and parallel text parsing. Select it explicitly in `zdrive.inp`:

```text
File Type = hypergraph, initial read = parallel
File Name = /path/to/hypergraph_without_the_final_dot_hg
```

Omitting `initial read = parallel` retains the legacy reader. `initial read =
chunks` retains its existing MatrixMarket meaning. This feature uses the native
Zoltan IBM/PaToH header conventions, not hMETIS headers. It does not alter PHG,
its callbacks, its parameters, or any other partitioning algorithm.

Each MPI rank reads a contiguous byte range, adjusted to complete newline
records. A prefix sum of record counts recovers the legacy implicit edge and
vertex IDs. Pins and weights are parsed from memory, avoiding the legacy
per-character wrapper and temporary token copies. Complete edges are exchanged
with `MPI_Alltoallv` according to the existing initial distribution rules.
Source rank order preserves the exact legacy local edge/pin order. No global
edge or pin array is assembled on rank zero. One rank needs no edge exchange.

Vertex weights still pass through rank zero and the existing vertex distributor;
this deliberately reuses its ownership and ordering rules. The supported initial
distributions are linear, cyclic, owner, and assignment-file distribution, including
`initial procs` subsets. As in the old reader, this does not turn graph connectivity
sidecars or geometry into a new file format. Sidecar handling remains in the driver.

The reader retains the driver's existing signed `int` counts. It checks MPI
message-count overflow rather than silently overflowing. This is not a solution
for inputs exceeding the driver's existing global count limits. IEEE finite
floating-point weights are converted with `strtod` and cast to `float`, matching
legacy `atof` conversion. Duplicate pins, empty edges under linear/cyclic
ownership, and CRLF are preserved. With multiple ranks, owner distribution of an
empty edge is rejected because the old distributor dereferences its first pin
and has undefined behavior. Records must be newline-terminated; trailing blank
lines are allowed. Compressed files, invalid IDs, missing/excess weights, invalid
numeric tokens, and nonfinite weights are rejected collectively. Malformed-input
acceptance quirks of the old `atoi`/`atof` parser are not reproduced.

`DRIVER: Input loading time` measures `read_mesh()` wall time, including reading,
parsing, initial distribution and mesh construction. It reports the maximum across
ranks, with its starting barrier and final reduction outside the measured region.
Distributed directory construction and partitioning happen afterwards. It is not
pure filesystem I/O time.

## Build and validation

Build the normal MPI C driver using either supported build system. For example,
from an empty build directory outside the source tree:

```sh
/path/to/Zoltan/configure --enable-mpi --with-mpi-compilers \
  --disable-f90interface --disable-zoltan-cppdriver --with-gnumake \
  CC=mpicc CXX=mpicxx CFLAGS='-O2 -g'
make -j 6 everything
```

The standalone probe compares actual local arrays from the serial reader and
legacy distributor against the parallel implementation, including all weight
bytes, indices, IDs, pins, and pin owners. It includes `dr_hg_io.c` in the test
translation unit to access its static legacy distributor without exporting that
helper from production code. Compile it from the source root (GNU-compatible
linker options below):

```sh
mpicc -O2 -g -ffunction-sections -fdata-sections \
  -I/path/to/build/src/include -Isrc/include -Isrc/driver -Isrc/ch \
  test/driver_hg_parallel/reader_compare.c src/driver/dr_hg_readfile.c \
  src/driver/dr_err.c src/driver/dr_util.c \
  src/ch/ch_dist_graph.c src/ch/ch_init_dist.c \
  -Wl,--gc-sections -L/path/to/build/src -lzoltan -lm \
  -o /path/to/build/reader_compare

python3 test/driver_hg_parallel/run.py \
  --driver /path/to/build/src/driver/zdrive.exe \
  --probe /path/to/build/reader_compare
```

Tests cover IBM/PaToH input, both ID bases, vector vertex/edge weights,
duplicates, empty edges, empty ranks, long records spanning byte partitions,
CRLF, trailing blank lines, initial distributions and process subsets, and
collective rejection of malformed input. Small PHG runs require identical final
vertex-to-part assignments. A separate Python evaluator recomputes cut-net,
connectivity (`km1`) and part weights from the original complete graph; it also
checks the configured balance bound. These metrics must equal across readers.
Both configured cut objectives are exercised. Test weights are exactly
representable binary fractions, so comparisons require exact equality.

## Benchmark the supplied webbase graph

```sh
python3 test/driver_hg_parallel/run.py \
  --driver /path/to/build/src/driver/zdrive.exe \
  --probe /path/to/build/reader_compare \
  --benchmark webbase-1M.mtx.zoltan.hg --ranks 1,8 \
  --check-partitions --report /tmp/zoltan-webbase-io-benchmark.json
```

This first compares the actual graph's rank-local arrays. Loading benchmarks use
one MPI rank per bound core, disable partitioning and text/debug output, discard
one warm-up per reader/rank count, alternate reader order, and record five
repetitions. The input remains on its original filesystem through a symlink.
These are cached end-to-end loading times, not cold storage bandwidth numbers.
`--check-partitions` additionally runs PHG with eight parts and identical
parameters/seeds through each reader, requires identical assignments, and
independently recomputes both objectives and balance from the original file.
Read the JSON report for individual timings and the compared metrics. Run this
inside a compute-node allocation to measure SuperMUC-NG rather than the local
machine. MPI launch arguments may need adapting to the site's launcher.
