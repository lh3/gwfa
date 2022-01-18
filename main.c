#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <zlib.h>
#include "gwfa.h"
#include "kseq.h"
#include "ketopt.h"
KSEQ_INIT(gzFile, gzread)

int main_ed(int argc, char *argv[])
{
	gzFile fp1, fp2;
	kseq_t *ks1, *ks2;
	ketopt_t o = KETOPT_INIT;
	gwf_graph_t *g;
	int c, s;
	void *km = 0;

	while ((c = ketopt(&o, argc, argv, 1, "", 0)) >= 0) {
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: test-ed <in1.fa> <in2.fa>\n");
		return 1;
	}

	fp1 = gzopen(argv[o.ind+0], "r");
	fp2 = gzopen(argv[o.ind+1], "r");
	assert(fp1 && fp2);
	ks1 = kseq_init(fp1);
	ks2 = kseq_init(fp2);
	kseq_read(ks1);
	kseq_read(ks2);

	g = (gwf_graph_t*)calloc(1, sizeof(*g));
	g->n_vtx = 1, g->n_arc = 0;
	g->seq = calloc(g->n_vtx, sizeof(*g->seq));
	g->seq[0] = strdup(ks1->seq.s);
	g->len = (uint32_t*)calloc(g->n_vtx, sizeof(*g->len));
	g->len[0] = ks1->seq.l;

	gwf_ed_index(km, g);
	s = gwf_ed(km, g, ks2->seq.l, ks2->seq.s, 0, 0);
	gwf_cleanup(km, g);
	printf("%d\n", s);

	kseq_destroy(ks1);
	kseq_destroy(ks2);
	gzclose(fp1);
	gzclose(fp2);
	return 0;
}


int main_test(int argc, char *argv[])
{
	gwf_graph_t *g;
	char *q;
	int32_t ql, s, v, v0, v1;
	void *km = 0;

	g = (gwf_graph_t*)calloc(1, sizeof(*g));

#define XTEST 2

#if XTEST == 1
	g->n_vtx = 1, g->n_arc = 0;
	g->seq = calloc(g->n_vtx, sizeof(*g->seq));
	g->seq[0] = strdup("GCTGCGATAGACCCTT");
	q = strdup("AGCTGCAGACCCTT");
	v0 = v1 = 0;
#elif XTEST == 2
	g->n_vtx = 3, g->n_arc = 4;
	g->seq = calloc(g->n_vtx, sizeof(*g->seq));
	g->arc = calloc(g->n_arc, sizeof(*g->arc));
	g->seq[0] = strdup("CA");
	g->seq[1] = strdup("T");
	g->seq[2] = strdup("TA");
	g->arc[0] = 0ULL<<32 | 1;
	g->arc[1] = 1ULL<<32 | 2;
	g->arc[2] = 1ULL<<32 | 0;
	g->arc[3] = 2ULL<<32 | 1;
	q = strdup("TATTA");
	v0 = 2, v1 = 2;
#endif

	ql = strlen(q);
	g->len = (uint32_t*)calloc(g->n_vtx, sizeof(*g->len));
	for (v = 0; v < g->n_vtx; ++v)
		g->len[v] = strlen(g->seq[v]);

	gwf_ed_index(km, g);
	s = gwf_ed(km, g, ql, q, v0, v1);

	printf("%d\n", s);

	return 0;
}

int main(int argc, char *argv[])
{
	return main_ed(argc, argv);
}
