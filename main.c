#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <zlib.h>
#include "gfa.h"
#include "gfa-priv.h"
#include "gwfa.h"
#include "ketopt.h"
#include "kalloc.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

gwf_graph_t *gwf_gfa2gwf(const gfa_t *gfa, uint32_t v0)
{
	int32_t i, k;
	gwf_graph_t *g;
	gfa_sub_t *sub;
	sub = gfa_sub_from(0, gfa, v0, 1<<30);
	GFA_CALLOC(g, 1);
	g->n_vtx = sub->n_v;
	g->n_arc = sub->n_a;
	GFA_MALLOC(g->len, g->n_vtx);
	GFA_MALLOC(g->src, g->n_vtx);
	GFA_MALLOC(g->seq, g->n_vtx);
	GFA_MALLOC(g->arc, g->n_arc);
	GFA_MALLOC(g->ol,  g->n_arc);
	for (i = k = 0; i < sub->n_v; ++i) {
		uint32_t v = sub->v[i].v, len = gfa->seg[v>>1].len, j;
		const gfa_seg_t *s = &gfa->seg[v>>1];
		g->len[i] = len;
		g->src[i] = v;
		GFA_MALLOC(g->seq[i], len + 1);
		if (v&1) {
			for (j = 0; j < len; ++j)
				g->seq[i][j] = gfa_comp_table[(uint8_t)s->seq[len - j - 1]];
		} else memcpy(g->seq[i], s->seq, len);
		g->seq[i][len] = 0; // null terminated for convenience
		for (j = 0; j < sub->v[i].n; ++j) {
			uint64_t a = sub->a[sub->v[i].off + j];
			g->arc[k] = (uint64_t)i<<32 | a>>32;
			g->ol[k] = gfa->arc[(uint32_t)a].ow;
			++k;
		}
		assert(k <= g->n_arc);
	}
	return g;
}

void gwf_free(gwf_graph_t *g)
{
	int32_t i;
	for (i = 0; i < g->n_vtx; ++i) free(g->seq[i]);
	free(g->len); free(g->seq); free(g->arc); free(g->src); free(g->ol); free(g);
}

void gwf_graph_print(FILE *fp, const gwf_graph_t *g)
{
	int32_t i;
	for (i = 0; i < g->n_vtx; ++i)
		fprintf(fp, "S\t%d\t%s\tLN:i:%d\n", i, g->seq[i], g->len[i]);
	for (i = 0; i < g->n_arc; ++i)
		fprintf(fp, "L\t%d\t+\t%d\t+\t%dM\n", (uint32_t)(g->arc[i]>>32), (uint32_t)g->arc[i], g->ol[i]);
}

int main(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *ks;
	ketopt_t o = KETOPT_INIT;
	gfa_t *gfa;
	gwf_graph_t *g;
	gwf_path_t path;
	int c, print_graph = 0, traceback = 0;
	uint32_t v0 = 0<<1|0; // first segment, forward strand
	uint32_t max_lag = 0;
	void *km = 0;
	char *sname = 0;

	while ((c = ketopt(&o, argc, argv, 1, "ptl:s:", 0)) >= 0) {
		if (c == 'p') print_graph = 1;
		else if (c == 'l') max_lag = atoi(o.arg);
		else if (c == 's') sname = o.arg;
		else if (c == 't') traceback = 1;
	}
	if ((!print_graph && argc - o.ind < 2) || (print_graph && argc == o.ind)) {
		fprintf(stderr, "Usage: gwf-test [options] <target.gfa|fa> <query.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -l INT    max lag behind the furthest wavefront; 0 to disable [0]\n");
		fprintf(stderr, "  -s STR    starting segment name [first]\n");
		fprintf(stderr, "  -t        report the alignment path\n");
		fprintf(stderr, "  -p        output GFA in the forward strand\n");
		return 1;
	}

	km = km_init();

	gfa = gfa_read(argv[o.ind]);
	assert(gfa);
	if (sname) {
		int32_t sid;
		sid = gfa_name2id(gfa, sname);
		if (sid < 0) fprintf(stderr, "ERROR: failed to find segment '%s'\n", sname);
		else v0 = sid<<1 | 0; // TODO: also allow to change the orientation
	}
	g = gwf_gfa2gwf(gfa, v0);
	if (print_graph) {
		gwf_graph_print(stdout, g);
		return 0; // free memory
	}
	gwf_ed_index(km, g);

	fp = gzopen(argv[o.ind+1], "r");
	assert(fp);
	ks = kseq_init(fp);
	while (kseq_read(ks) >= 0) {
		int32_t s;
		s = gwf_ed(km, g, ks->seq.l, ks->seq.s, 0, -1, max_lag, traceback, &path);
		if (traceback) {
			int32_t i, last_len = -1, len = 0;
			printf("%s\t%ld\t0\t%ld\t+\t", ks->name.s, ks->seq.l, ks->seq.l);
			for (i = 0; i < path.nv; ++i) {
				uint32_t v = g->src[path.v[i]];
				printf("%c%s", "><"[v&1], gfa->seg[v>>1].name);
				last_len = gfa->seg[v>>1].len;
				len += last_len;
			}
			printf("\t%d\t0\t%d\t%d\n", len, len - (last_len - path.end_off) + 1, path.s);
		} else printf("%s\t%d\n", ks->name.s, s);
	}
	kseq_destroy(ks);
	gzclose(fp);
	gfa_destroy(gfa);

	gwf_cleanup(km, g);
	gwf_free(g);
	km_destroy(km);
	return 0;
}
