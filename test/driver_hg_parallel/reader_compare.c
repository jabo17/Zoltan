/* Standalone integration probe: compile with function sections so the unused
 * driver entry points are discarded.  Including this file exposes the legacy
 * static distributor to the test without changing its production visibility. */
#include "../../src/driver/dr_hg_io.c"

int Debug_Driver = 0;
int Debug_Chaco_Input = 0;
int Chaco_In_Assign_Inv = 0;

struct graph {
  int nv, gnv, ne, gne, pins, first, vd, ed, base;
  int *index, *vertex, *gid, *owner;
  float *vw, *ew;
};

static void load(struct graph *g, const char *file, int parallel,
                 PARIO_INFO_PTR pio, short **assignments)
{
  int rank, dim=0, graphwd=0, *xadj=NULL, *adj=NULL;
  float *x=NULL, *y=NULL, *z=NULL, *graphw=NULL;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (parallel) {
    if (!HG_readfile_parallel(file, &g->nv, &g->gne, &g->pins,
         &g->first, &g->ne, &g->index, &g->vertex, &g->vd, &g->vw,
         &g->ed, &g->ew, &g->base)) MPI_Abort(MPI_COMM_WORLD, 2);
  } else if (!rank) {
    ZOLTAN_FILE *f=ZOLTAN_FILE_open((char *)file, "r", STANDARD);
    if (!f || HG_readfile(rank, f, &g->nv, &g->ne, &g->pins,
         &g->index, &g->vertex, &g->vd, &g->vw, &g->ed, &g->ew,
         &g->base) != ZOLTAN_OK) MPI_Abort(MPI_COMM_WORLD, 2);
    ZOLTAN_FILE_close(f);
  }
  MPI_Bcast(&g->base, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (!chaco_dist_graph(MPI_COMM_WORLD, pio, 0, &g->gnv, &g->nv,
       &xadj, &adj, &g->vd, &g->vw, &graphwd, &graphw, &dim,
       &x, &y, &z, assignments)) MPI_Abort(MPI_COMM_WORLD, 2);
  if (parallel) {
    if (!HG_dist_parallel(pio, g->base, g->gne, &g->ne, g->first,
          &g->gid, &g->index, &g->vertex, &g->owner, g->ed,
          &g->ew, *assignments)) MPI_Abort(MPI_COMM_WORLD, 2);
  } else if (!dist_hyperedges(MPI_COMM_WORLD, pio, 0, g->base, g->gnv,
       &g->gne, &g->ne, &g->gid, &g->index, &g->vertex, &g->owner,
       &g->ed, &g->ew, *assignments)) MPI_Abort(MPI_COMM_WORLD, 2);
}

static int same(const void *a, const void *b, size_t count, size_t size)
{ return !count || !memcmp(a, b, count*size); }

int main(int argc, char **argv)
{
  struct graph old={0}, parallel={0};
  struct Parallel_IO pio={0};
  int rank, ok, all, nranks, i;
  short *assignments=NULL;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);
  if (argc != 4) MPI_Abort(MPI_COMM_WORLD, 2);
  pio.init_dist_type=atoi(argv[2]);
  pio.init_dist_procs=atoi(argv[3]);
  if (pio.init_dist_type==INITIAL_FILE && !rank) {
    /* The fixture uses 17 vertices; exercise assignment-file ownership too. */
    assignments=(short *)malloc(17*sizeof(short));
    for (i=0; i<17; ++i) assignments[i]=(short)((i*7)%nranks);
  }
  load(&old, argv[1], 0, &pio, &assignments);
  load(&parallel, argv[1], 1, &pio, &assignments);
  ok=old.nv==parallel.nv && old.gnv==parallel.gnv && old.ne==parallel.ne &&
     old.gne==parallel.gne && old.vd==parallel.vd && old.ed==parallel.ed &&
     old.base==parallel.base;
  if (ok) ok=same(old.index, parallel.index, old.ne+1, sizeof(int)) &&
     same(old.gid, parallel.gid, old.ne, sizeof(int)) &&
     same(old.vertex, parallel.vertex, old.index[old.ne], sizeof(int)) &&
     same(old.owner, parallel.owner, old.index[old.ne], sizeof(int)) &&
     same(old.vw, parallel.vw, (size_t)old.nv*old.vd, sizeof(float)) &&
     same(old.ew, parallel.ew, (size_t)old.ne*old.ed, sizeof(float));
  if (!ok) fprintf(stderr, "Rank %d: reader mismatch\n", rank);
  MPI_Allreduce(&ok, &all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (!rank && all) puts("Exact per-rank reader equality: PASS");
  /* One-shot probe; explicitly release output arrays for sanitizer use. */
  free(old.index); free(old.vertex); free(old.gid); free(old.owner); free(old.vw); free(old.ew);
  free(parallel.index); free(parallel.vertex); free(parallel.gid); free(parallel.owner); free(parallel.vw); free(parallel.ew);
  free(assignments);
  MPI_Finalize();
  return all ? 0 : 1;
}
