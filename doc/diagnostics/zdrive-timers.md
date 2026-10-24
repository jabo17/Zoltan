# Diagnostic zdrive timers

Branch `diagnostic/zdrive-timers` starts from `feature/parallel-io`.
Rebuild zdrive from this branch using your existing build configuration.
Timers are always enabled on this diagnostic branch; main is unchanged.

Use the usual input, including:

```
File Type = hypergraph, initial read = parallel
Compression = uncompressed
Zoltan Parameters = PHG_USE_TIMERS=2
```

The last parameter enables existing PHG internal timers, including its Build
stage (construction of the algorithm's hypergraph from driver callbacks).
It is optional and separate from the new driver timers. No partitioning
algorithm or parameter defaults are changed by this branch.

Capture both output streams:

```
mpiexec --timeout 600 --map-by ppr:48:node:PE=1 --bind-to core -n 48 \
  /path/to/zdrive.exe zdrive.inp > run.log 2>&1
```

The configuration summary explicitly confirms `Parallel native hypergraph
reader (MPI-IO) enabled`. Timing BEGIN/END lines are printed on rank zero
to stderr and flushed immediately. Scheduler/filesystem forwarding may still
delay visibility. A BEGIN without its END identifies the active section.
END reports minimum, mean, and maximum elapsed seconds across ranks.

| Section | Includes |
|---|---|
| driver_initialization | Initialization after MPI_Init, command-file reading, callback setup |
| input_total | Complete read_mesh call; includes all nested input sections below |
| input_read_and_parse | Reader and auxiliary input files, before input distribution |
| parallel_header | Rank-zero header read and broadcast |
| parallel_file_read | MPI file open, record boundaries, buffer allocation, byte reads, close |
| parallel_parse_and_validate | Record counting, parsing, allocation, validation and associated collectives |
| parallel_vertex_weight_gather | Gathering vertex weights for the legacy vertex distributor |
| input_distribution | Vertex/edge distribution and pin-owner construction |
| driver_mesh_construction | Conversion into the driver's mesh/element arrays and cleanup of temporary input arrays |
| initial_distributed_directory | Creation and update of the directory mapping vertices to ranks |
| initial_evaluation | Driver and Zoltan evaluations before partitioning; disabled at zdrive debug level 0 |
| partition_including_algorithm_setup | Zoltan_LB_Partition, including PHG Build and partitioning stages |
| partition_result_processing | Driver handling of partition results before migration |
| migration_including_barrier | Existing migration section, including its entry barrier |
| postmigration_checks_and_zoltan_copy | Fixed-object checks and the driver's Zoltan copy tests |
| final_evaluation | Driver and Zoltan evaluations after migration |
| cleanup_and_output | Destruction, result output and freeing input data |
| total_after_mpi_init | Application runtime after MPI_Init and before MPI_Finalize |

Sections are nested: do not add parent times to their children. Reading times
include allocation and MPI synchronization; `parallel_file_read` is not a
measurement of storage-device throughput alone. There are no added entry
barriers, but END uses collective reductions and therefore can affect timing
and synchronization slightly. Different rank maxima need not add up exactly.
No completed timing is promised on error/abort paths.

MPI launch, MPI_Init and MPI_Finalize are excluded from total_after_mpi_init.
Compare that total with your surrounding mpiexec timer to identify time spent
outside the measured driver. A delay before the first BEGIN is before or in
MPI_Init (or delayed log delivery).

Validation: rebuilt the MPI Autotools driver and ran the serial and parallel
readers on webbase with two ranks, two parts, seed 2, and no refinement.
Both produced final CUTL 447491.000; all expected section reports completed.
