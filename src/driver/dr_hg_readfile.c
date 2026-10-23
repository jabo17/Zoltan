/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */

#include <stdlib.h>
#include <ctype.h>
#include "zoltan.h"

#include "dr_const.h"
#include "dr_hg_readfile.h"
#include "dr_util_const.h"
#include "dr_compress_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#define BUF_LEN 1000
#define ERROR(proc,yo,msg,err) \
 {printf("Proc(%d) From(%s): %s\n",proc,yo,msg);return(err);}

static int   readfile  (int, ZOLTAN_FILE*, int*, int*, int*, int**, int**, int*,
 float**, int*, float**, int*);

static int nextstr (ZOLTAN_FILE* f, char *string);


/* These routines work with multiple (vector) vertex and hyperedge weights */

/*****************************************************************************/

int HG_readfile (
  int Proc,
  ZOLTAN_FILE* f,
  int *nVtx, int *nEdge, int *nPins,
  int **hindex,  int **hvertex,
  int *vwgt_dim, float **vwgt,
  int *ewgt_dim, float **ewgt,
  int *base)
{
  char string[BUF_LEN];
  const char *yo = "HG_readfile";
  int i, code;

  /* Initialize return values in case of error. */
  *nVtx  = *nEdge  = *nPins = *vwgt_dim = *ewgt_dim = *base = 0;
  *hindex = *hvertex = NULL;
  *vwgt = *ewgt = NULL;

  /* Skip comments. */
  while (1) {
    ZOLTAN_FILE_gets (string, BUF_LEN-1, f);

    if (*string != '%')
      break;
  }

  /* Extended format: nVtx nEdge nPins NM where both N&M are [0,9]. N is the
   * number of dimension for vertex weights and M is number of dimensions for
   * the edge weights. */
  if (sscanf(string, "%d %d %d %d", nVtx,nEdge,nPins,&code) != 4)
    ERROR(Proc, yo, "Unrecognized file format.", ZOLTAN_FATAL);

  if (*nVtx <= 1) {                  /* wrong guess if nVtx really is one! */
    i = sscanf(string, "%d %d %d %d %d %d", 
                       base, nVtx, nEdge, nPins, &code, vwgt_dim);
    if (i == 4)
      code = 0;
    else if (i < 4)
      ERROR(Proc, yo, "Error in PaToH file format.", ZOLTAN_FATAL);

    switch (code) {
    case 0:    /* No weights */
      *vwgt_dim = 0;
      *ewgt_dim = 0;
      break;
    case 1: /* Vwgts only; default is 1 vwgt unless 6 ints on input line */
      if (i == 5) *vwgt_dim = 1;
      *ewgt_dim = 0;
      break;
    case 2: /* Ewgts only; PaToH supports only one ewgt. */
      *vwgt_dim = 0;
      *ewgt_dim = 1;
      break;
    case 3: /* Ewgts & vwgts; default is 1 vwgt unless 6 ints on input line */
      if (i == 5) *vwgt_dim = 1;
      *ewgt_dim = 1;
      break;
    default: /* Invalid code */
      ERROR(Proc, yo, "Error in PaToH file format:  Illegal PaToH code.",
            ZOLTAN_FATAL);
      break;
    }

    return readfile(Proc, f, nVtx, nEdge, nPins, hindex, hvertex, vwgt_dim,
                    vwgt, ewgt_dim, ewgt, base);          /* patoh format */
  }
  else {
    *base = 1;
    if (sscanf(string, "%d %d %d %1d%1d",
                        nVtx,nEdge,nPins,vwgt_dim,ewgt_dim) != 5)
      ERROR(Proc, yo, "Error in HG file format.", ZOLTAN_FATAL);
    return readfile(Proc, f, nVtx, nEdge, nPins, hindex, hvertex, vwgt_dim,
                     vwgt, ewgt_dim, ewgt, base);          /* IBM format */
  }
}


/*****************************************************************************/
/* Read IBM & Patoh formats */

static int readfile (
 int Proc,
 ZOLTAN_FILE* f,
 int *nVtx, int *nEdge, int *nPins,
 int **index,   int **vertex,
 int *vwgt_dim, float **vwgt,
 int *ewgt_dim, float **ewgt,
 int *base)
{
  int i, j, pin;
  char string[BUF_LEN];
  float *p;
  const char *yo = "readfile";

  if (!(*index  = (int*) malloc ((*nEdge+1) * sizeof(int)))
   || !(*vertex = (int*) malloc ((*nPins    * sizeof(int)))))
      ERROR (Proc, yo, "Insufficient memory", ZOLTAN_MEMERR);

  if (*vwgt_dim > 0
   && !(*vwgt = (float*) malloc (*nVtx * *vwgt_dim * sizeof (float))))
      ERROR (Proc, yo, "Insufficient memory", ZOLTAN_MEMERR)

  if (*ewgt_dim > 0
   && !(*ewgt = (float*) malloc (*nEdge * *ewgt_dim * sizeof (float))))
      ERROR (Proc, yo, "Insufficient memory", ZOLTAN_MEMERR)

  pin = 0;
  p = *ewgt;
  for (i = 0; i < *nEdge; i++) {
     (*index)[i] = pin;
     for (j = 0; j < *ewgt_dim; j++) {
        if (nextstr (f, string) != 1)
           ERROR (Proc, yo, "Wrong number of edge weights", ZOLTAN_FATAL)
        *p++ = atof (string);
     }
     while (nextstr (f,string) == 1)
        (*vertex)[pin++] = atoi(string);
  }
  (*index)[i] = pin;
  if (pin != *nPins)
     ERROR (Proc, yo, "Assertion Error in readfile.", ZOLTAN_FATAL);

  p = *vwgt;
  if (*vwgt_dim > 0)
     for (i = 0; i < *nVtx; i++) {
        for (j = 0; nextstr (f, string) == 1 && j < *vwgt_dim; j++)
           *p++ = atof (string);
        if (j != *vwgt_dim)
	  ERROR (Proc, yo, "Wrong number of vertex weights", ZOLTAN_FATAL)
     }

  if (nextstr (f, string) == 1)
     ERROR (Proc, yo, "Input file is longer than expected", ZOLTAN_FATAL)
  return ZOLTAN_OK;
}


/*****************************************************************************/
/* Uses character reads to read arbitrary long hyperedge definition lines.   */
static int nextstr (ZOLTAN_FILE* f, char *string)
{
  char ch;

  *string = 0;
  while (1) {
    ch = ZOLTAN_FILE_getc(f);
    if (ch == EOF)          return 0;
    if (ch == '\n')         return -1;
    if (isspace(ch) == 0)   break;
  }

  while (1)  {
    if (ch == EOF || isspace(ch)) {
      *string = 0;
      ZOLTAN_FILE_ungetc (ch, f);  /* ??? */
      return 1;
    }
    *string++ = ch;
    ch = ZOLTAN_FILE_getc(f);
  }
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

/* Parallel .hg input is intentionally confined to the driver.  The byte
 * partitions contain whole physical records; their prefix sums give the same
 * implicit edge/vertex IDs as the serial parser. */
#ifdef HAVE_MPI
#include <limits.h>
#include <errno.h>
#include <math.h>
#include "ch_init_dist_const.h"

static int hg_all_ok(int ok)
{
  int all;
  MPI_Allreduce(&ok, &all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  return all;
}

/* Independent, large buffered reads, including files larger than INT_MAX. */
static int hg_read_bytes(MPI_File f, MPI_Offset off, char *buf, size_t bytes)
{
  while (bytes) {
    MPI_Status status;
    int got, n = bytes > 4*1024*1024 ? 4*1024*1024 : (int)bytes;
    if (MPI_File_read_at(f, off, buf, n, MPI_BYTE, &status) != MPI_SUCCESS)
      return 0;
    MPI_Get_count(&status, MPI_BYTE, &got);
    if (got != n) return 0;
    off += n; buf += n; bytes -= n;
  }
  return 1;
}

static int hg_record_boundary(MPI_File f, MPI_Offset cut,
                              MPI_Offset end, MPI_Offset *result)
{
  char buf[65536];
  MPI_Offset pos = cut - 1;
  if (cut >= end) { *result = end; return 1; }
  while (pos < end) {
    size_t n = end-pos > (MPI_Offset)sizeof(buf) ? sizeof(buf) : (size_t)(end-pos);
    char *nl;
    if (!hg_read_bytes(f, pos, buf, n)) return 0;
    nl = (char *)memchr(buf, '\n', n);
    if (nl) { *result = pos + (nl-buf) + 1; return 1; }
    pos += n;
  }
  *result = end;
  return 1;
}

static char *hg_skip_space(char *s)
{
  while (*s && isspace((unsigned char)*s)) ++s;
  return s;
}

static int hg_weight(char **s, float *out)
{
  char *end;
  double value;
  *s = hg_skip_space(*s);
  value = strtod(*s, &end); /* Same double-to-float conversion as atof. */
  if (end == *s || (*end && !isspace((unsigned char)*end)) || !isfinite(value))
    return 0;
  *out = (float)value;
  if (!isfinite(*out)) return 0;
  *s = end;
  return 1;
}

int HG_readfile_parallel(const char *filename, int *nv, int *gne, int *np,
 int *first_edge, int *ne, int **index, int **vertex, int *vd, float **vw,
 int *ed, float **ew, int *base)
{
  int rank, nranks, header[6] = {0}, ok = 1, i, local_pins = 0;
  int local_edges = 0, local_weights = 0, *counts = NULL, *displs = NULL;
  int capacity = 0, *pins = NULL;
  long long lines = 0, first = 0, total_lines, total_pins, pin_count, expected, used;
  MPI_File f;
  MPI_Offset body = 0, size = 0, begin, end;
  char *buf = NULL, *line, *stop;
  float *local_vw = NULL;
  *index = *vertex = NULL; *vw = *ew = NULL;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);
  if (!rank) {
    char text[BUF_LEN];
    FILE *in = fopen(filename, "r");
    int code = 0, n;
    if (!in) ok = 0;
    if (ok) {
      do { if (!fgets(text, sizeof(text), in)) { ok = 0; break; } }
      while (text[0] == '%');
      if (ok && sscanf(text, "%d %d %d %d", nv, gne, np, &code) != 4) ok = 0;
      if (ok && *nv <= 1) {
        n = sscanf(text, "%d %d %d %d %d %d", base, nv, gne, np, &code, vd);
        if (n == 4) code = 0;
        *ed = 0;
        if (n < 4 || code < 0 || code > 3) ok = 0;
        if (code == 0 || code == 2) *vd = 0;
        else if (n == 5) *vd = 1;
        if (code == 2 || code == 3) *ed = 1;
      } else if (ok) {
        *base = 1;
        if (sscanf(text, "%d %d %d %1d%1d", nv, gne, np, vd, ed) != 5) ok = 0;
      }
      body = (MPI_Offset)ftell(in);
      fclose(in);
    }
    if (ok && (*nv < 0 || *gne < 0 || *np < 0 || *vd < 0 || *ed < 0 ||
               *vd > MAX_CPU_WGTS || (*vd && *nv > INT_MAX / *vd) ||
               (*base != 0 && *base != 1))) ok = 0;
    header[0]=*nv; header[1]=*gne; header[2]=*np;
    header[3]=*vd; header[4]=*ed; header[5]=*base;
  }
  MPI_Bcast(&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (!ok) goto fail;
  MPI_Bcast(header, 6, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&body, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);
  *nv=header[0]; *gne=header[1]; *np=header[2];
  *vd=header[3]; *ed=header[4]; *base=header[5];
  if (MPI_File_open(MPI_COMM_WORLD, (char *)filename, MPI_MODE_RDONLY,
                    MPI_INFO_NULL, &f) != MPI_SUCCESS) goto fail;
  MPI_File_set_errhandler(f, MPI_ERRORS_RETURN);
  ok = MPI_File_get_size(f, &size) == MPI_SUCCESS && size >= body;
  if (!hg_all_ok(ok)) { MPI_File_close(&f); goto fail; }
  begin = body + (size-body)/nranks*rank + (size-body)%nranks*rank/nranks;
  end = body + (size-body)/nranks*(rank+1) + (size-body)%nranks*(rank+1)/nranks;
  if (rank) ok = hg_record_boundary(f, begin, size, &begin);
  if (ok && rank+1 < nranks) ok = hg_record_boundary(f, end, size, &end);
  if (ok) {
    buf = (char *)malloc((size_t)(end-begin)+1);
    ok = buf != NULL;
  }
  if (ok) ok = hg_read_bytes(f, begin, buf, (size_t)(end-begin));
  MPI_File_close(&f);
  if (!hg_all_ok(ok)) goto fail;
  buf[end-begin] = 0;
  /* Require terminated records, as the legacy nextstr EOF path drops the
   * last unterminated token.  Reject rather than change that behavior. */
  ok = end == begin || buf[end-begin-1] == '\n';
  for (line=buf; line < buf+(end-begin); ++line) if (*line == '\n') ++lines;
  MPI_Exscan(&lines, &first, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  if (!rank) first = 0;
  MPI_Allreduce(&lines, &total_lines, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  expected = (long long)*gne + (*vd ? *nv : 0);
  if (total_lines < expected) ok = 0;
  *first_edge = first < *gne ? (int)first : *gne;
  local_edges = first >= *gne ? 0 : (int)((first+lines < *gne ? first+lines : *gne)-first);
  used = first >= expected ? 0 : (first+lines < expected ? lines : expected-first);
  if (used-local_edges > INT_MAX) ok = 0;
  local_weights = ok ? (int)(used-local_edges) : 0;
  if (*vd && local_weights > INT_MAX / *vd) ok = 0;
  if (*ed && local_edges > INT_MAX / *ed) ok = 0;
  if (!hg_all_ok(ok)) goto fail;
  *ne = local_edges;
  *index = (int *)calloc((size_t)local_edges+1, sizeof(int));
  *ew = (float *)malloc(((size_t)local_edges* *ed+1)*sizeof(float));
  local_vw = (float *)malloc(((size_t)local_weights* *vd+1)*sizeof(float));
  ok = *index && *ew && local_vw;
  if (!hg_all_ok(ok)) goto fail;
  line = buf;
  for (i=0; i<local_edges && ok; ++i) {
    char *s = line;
    int j;
    stop = strchr(line, '\n'); *stop = 0; line = stop+1;
    (*index)[i] = local_pins;
    for (j=0; j<*ed && ok; ++j) ok = hg_weight(&s, *ew+(size_t)i* *ed+j);
    while (ok && *(s=hg_skip_space(s))) {
      long pin;
      char *tail;
      errno = 0;
      pin = strtol(s, &tail, 10);
      if (errno || tail == s || (*tail && !isspace((unsigned char)*tail)) ||
          pin < *base || pin-(long)*base >= *nv || local_pins == INT_MAX) { ok=0; break; }
      if (local_pins == capacity) {
        int next = capacity > INT_MAX/2 ? INT_MAX : (capacity ? capacity*2 : 1024);
        int *grown = (int *)realloc(pins, (size_t)next*sizeof(int));
        if (!grown) { ok=0; break; }
        pins=grown; capacity=next;
      }
      pins[local_pins++] = (int)pin;
      s = tail;
    }
  }
  (*index)[local_edges] = local_pins;
  for (i=0; i<local_weights && ok; ++i) {
    char *s = line;
    int j;
    stop=strchr(line, '\n'); *stop=0; line=stop+1;
    for (j=0; j<*vd && ok; ++j) ok=hg_weight(&s, local_vw+(size_t)i* *vd+j);
    if (*hg_skip_space(s)) ok=0;
  }
  /* The serial reader permits trailing blank lines (including webbase's).
   * They are not vertex records and must not enter the weight gather. */
  if (ok && *hg_skip_space(line)) ok=0;
  pin_count = local_pins;
  MPI_Allreduce(&pin_count, &total_pins, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  if (total_pins != *np) ok=0;
  free(buf); buf=NULL;
  counts=(int *)malloc((size_t)nranks*sizeof(int));
  displs=(int *)malloc((size_t)nranks*sizeof(int));
  if (!rank && *vd) *vw=(float *)malloc((size_t)*nv* *vd*sizeof(float));
  ok = ok && counts && displs && (rank || !*vd || !*nv || *vw);
  if (!hg_all_ok(ok)) goto fail;
  /* Only the relatively small vertex-weight array passes through the existing
   * vertex distributor.  No edge or pin array is gathered onto rank zero. */
  local_weights *= *vd;
  MPI_Allgather(&local_weights, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD);
  displs[0]=0;
  for (i=1; i<nranks; ++i) {
    if (counts[i-1] > INT_MAX-displs[i-1]) ok=0;
    else displs[i]=displs[i-1]+counts[i-1];
  }
  if (!hg_all_ok(ok)) goto fail;
  MPI_Gatherv(local_vw, local_weights, MPI_FLOAT, *vw, counts, displs,
              MPI_FLOAT, 0, MPI_COMM_WORLD);
  free(local_vw); free(counts); free(displs);
  *vertex=pins;
  return 1;
fail:
  free(buf); free(pins); free(local_vw); free(counts); free(displs);
  free(*index); *index=NULL; free(*ew); *ew=NULL; free(*vw); *vw=NULL;
  if (!rank) fprintf(stderr, "Parallel HG reader: invalid/unsupported input, read error, or insufficient memory: %s\n", filename);
  return 0;
}

int HG_dist_parallel(PARIO_INFO_PTR pio, int base, int global_edges, int *ne,
 int first_edge, int **gid, int **index, int **vertex, int **owner,
 int ed, float **ew, short *assignments)
{
  int nranks, rank, ndistrib, type, h, p, i, total=0, received=0, ok=1;
  int *counts=NULL, *recvcounts=NULL, *displs=NULL, *rdispls=NULL, *cursor=NULL;
  int *send=NULL, *recv=NULL, *dest=NULL;
  int *newindex=NULL, *newgid=NULL, *newvertex=NULL, *newowner=NULL;
  float *newew=NULL;
  int nedge=0, npin=0;
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (nranks == 1) {
    /* No exchange or copy is needed; match dist_hyperedges' serial case. */
    int pins = (*index)[*ne];
    *gid = (int *)malloc(((size_t)*ne+1)*sizeof(int));
    *owner = (int *)calloc((size_t)pins+1, sizeof(int));
    if (!*gid || !*owner) { free(*gid); free(*owner); *gid=*owner=NULL; return 0; }
    for (h=0; h<*ne; ++h) (*gid)[h]=h+base;
    return 1;
  }
  ndistrib=pio->init_dist_procs;
  if (ndistrib <= 0 || ndistrib > nranks) ndistrib=nranks;
  type=pio->init_dist_type == INITIAL_FILE ? INITIAL_LINEAR : pio->init_dist_type;
  counts=(int *)calloc(nranks, sizeof(int)); recvcounts=(int *)calloc(nranks, sizeof(int));
  displs=(int *)calloc(nranks, sizeof(int)); rdispls=(int *)calloc(nranks, sizeof(int));
  cursor=(int *)calloc(nranks, sizeof(int)); dest=(int *)malloc(((size_t)*ne+1)*sizeof(int));
  ok=counts && recvcounts && displs && rdispls && cursor && dest && sizeof(float)==sizeof(int);
  if (!hg_all_ok(ok)) goto fail;
  for (h=0; h<*ne; ++h) {
    int n=(*index)[h+1]-(*index)[h];
    if (nranks==1) p=0;
    else if (type==INITIAL_CYCLIC) p=(first_edge+h)%ndistrib;
    else if (type==INITIAL_LINEAR)
      /* Match the legacy float expression, including its rounding. */
      p=(int)((float)(first_edge+h)*(float)ndistrib/(float)global_edges);
    else if (type==INITIAL_OWNER && n)
      p=ch_dist_proc((*vertex)[(*index)[h]], assignments, base);
    else { ok=0; break; } /* Legacy OWNER of an empty edge is undefined. */
    if (p<0 || p>=nranks || n > INT_MAX-2-ed || counts[p]>INT_MAX-2-ed-n) {ok=0; break;}
    dest[h]=p; counts[p]+=2+ed+n;
  }
  if (!hg_all_ok(ok)) goto fail;
  MPI_Alltoall(counts, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD);
  for (p=0; p<nranks; ++p) {
    displs[p]=total; rdispls[p]=received;
    if (counts[p]>INT_MAX-total || recvcounts[p]>INT_MAX-received) {ok=0; break;}
    total+=counts[p]; received+=recvcounts[p]; cursor[p]=displs[p];
  }
  if (!hg_all_ok(ok)) goto fail;
  send=(int *)malloc(((size_t)total+1)*sizeof(int));
  recv=(int *)malloc(((size_t)received+1)*sizeof(int));
  if (!hg_all_ok(send && recv)) goto fail;
  for (h=0; h<*ne; ++h) {
    int n=(*index)[h+1]-(*index)[h], *out=send+cursor[dest[h]];
    *out++=first_edge+h+base; *out++=n;
    if (ed) memcpy(out, *ew+(size_t)h*ed, (size_t)ed*sizeof(float));
    out+=ed;
    if (n) memcpy(out, *vertex+(*index)[h], (size_t)n*sizeof(int));
    cursor[dest[h]]+=2+ed+n;
  }
  MPI_Alltoallv(send, counts, displs, MPI_INT, recv, recvcounts, rdispls, MPI_INT, MPI_COMM_WORLD);
  free(send); send=NULL;
  /* Source ranks own consecutive input ranges.  MPI_Alltoallv concatenates
   * sources in rank order, retaining exactly the legacy ascending edge order. */
  for (i=0; i<received; i+=2+ed+recv[i+1]) { ++nedge; npin+=recv[i+1]; }
  newindex=(int *)malloc(((size_t)nedge+1)*sizeof(int));
  newgid=(int *)malloc(((size_t)nedge+1)*sizeof(int));
  newvertex=(int *)malloc(((size_t)npin+1)*sizeof(int));
  newowner=(int *)malloc(((size_t)npin+1)*sizeof(int));
  newew=(float *)malloc(((size_t)nedge*ed+1)*sizeof(float));
  if (!hg_all_ok(newindex && newgid && newvertex && newowner && newew)) goto fail;
  h=0; npin=0;
  for (i=0; i<received;) {
    int n=recv[i+1];
    newindex[h]=npin; newgid[h]=recv[i];
    if (ed) memcpy(newew+(size_t)h*ed, recv+i+2, (size_t)ed*sizeof(float));
    if (n) memcpy(newvertex+npin, recv+i+2+ed, (size_t)n*sizeof(int));
    npin+=n; ++h; i+=2+ed+n;
  }
  newindex[h]=npin;
  for (i=0; i<npin; ++i) newowner[i]=ch_dist_proc(newvertex[i], assignments, base);
  free(*index); free(*vertex); free(*ew);
  *index=newindex; *gid=newgid; *vertex=newvertex; *owner=newowner; *ew=newew; *ne=nedge;
  free(counts); free(recvcounts); free(displs); free(rdispls); free(cursor); free(dest); free(recv);
  return 1;
fail:
  free(counts); free(recvcounts); free(displs); free(rdispls); free(cursor); free(dest); free(send); free(recv);
  free(newindex); free(newgid); free(newvertex); free(newowner); free(newew);
  if (!rank) fprintf(stderr, "Parallel HG redistribution: unsupported distribution, count overflow, or insufficient memory\n");
  return 0;
}
#else
int HG_readfile_parallel(const char *filename, int *nv, int *gne, int *np,
 int *first_edge, int *ne, int **index, int **vertex, int *vd, float **vw,
 int *ed, float **ew, int *base)
{
  fprintf(stderr, "Parallel HG reader requires an MPI build\n");
  return 0;
}
int HG_dist_parallel(PARIO_INFO_PTR pio, int base, int global_edges, int *ne,
 int first_edge, int **gid, int **index, int **vertex, int **owner,
 int ed, float **ew, short *assignments)
{ return 0; }
#endif
