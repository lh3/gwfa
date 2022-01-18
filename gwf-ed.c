#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "gwfa.h"
#include "kalloc.h"
#include "ksort.h"

/**********************
 * Indexing the graph *
 **********************/

#define generic_key(x) (x)
KRADIX_SORT_INIT(gwf64, uint64_t, generic_key, 8)

void gwf_ed_index_arc_core(uint64_t *idx, uint32_t n_vtx, uint32_t n_arc, uint64_t *arc)
{
	uint32_t i, st;
	radix_sort_gwf64(arc, arc + n_arc);
	for (st = 0, i = 1; i <= n_arc; ++i) {
		if (i == n_arc || arc[i]>>32 != arc[st]>>32) {
			uint32_t v = arc[st]>>32;
			assert(v < n_vtx);
			idx[v] = (uint64_t)st << 32 | (i - st);
			st = i;
		}
	}
}

void gwf_ed_index(void *km, gwf_graph_t *g)
{
	KMALLOC(km, g->aux, g->n_vtx);
	gwf_ed_index_arc_core(g->aux, g->n_vtx, g->n_arc, g->arc);
}

void gwf_cleanup(void *km, gwf_graph_t *g)
{
	kfree(km, g->aux);
	g->aux = 0;
}

/**************************************
 * Graph WaveFront with edit distance *
 **************************************/

#include "khashl.h" // make it compatible with kalloc
#include "kdq.h"

KHASHL_INIT(KH_LOCAL, gwf_set64_t, gwf_set64, uint64_t, kh_hash_dummy, kh_eq_generic)

typedef struct { // a diagonal
	uint64_t vd; // NB: this wastes 4 bytes due to memory alignment
	int32_t k;
} gwf_diag_t;

#define ed_key(x) ((x).vd)
KRADIX_SORT_INIT(gwf_ed, gwf_diag_t, ed_key, 8)

KDQ_INIT(gwf_diag_t)

static inline void gwf_ed_push(kdq_t(gwf_diag_t) *a, uint32_t v, int32_t d, int32_t k)
{
	gwf_diag_t t;
	t.vd = (uint64_t)v<<32 | (80000000LL + d), t.k = k;
	*kdq_pushp(gwf_diag_t, a) = t;
}

static inline int32_t gwf_ed_update(gwf_diag_t *p, uint32_t v, int32_t d, int32_t k)
{
	uint64_t vd = (uint64_t)v<<32 | (80000000LL + d);
	if (p->vd == vd) {
		p->k = p->k > k? p->k : k;
		return 0;
	}
	return 1;
}

// for each diagonal, remove elements that are on on the wave front
static int32_t gwf_ed_dedup(void *km, int32_t n_a, gwf_diag_t *a)
{
	int32_t i, n, st;
	// TODO: if sorting too slow, add a linear pass to reduce the initial size and/or use stable radix sort
	radix_sort_gwf_ed(a, a + n_a);
	for (i = 1, st = 0, n = 0; i <= n_a; ++i) {
		if (i == n_a || a[i].vd != a[st].vd) {
			int32_t j, max_j = st;
			for (j = st + 1; j < i; ++j) // choose the far end (i.e. the wave front)
				if (a[max_j].k < a[j].k) max_j = j;
			a[n++] = a[max_j];
			st = i;
		}
	}
#ifdef GWF_DEBUG
	printf("[%s] %d -> %d\n", __func__, n_a, n);
#endif
	return n;
}

gwf_diag_t *gwf_ed_extend(void *km, const gwf_graph_t *g, int32_t ql, const char *q, int32_t v1, int32_t *is_end, int32_t *n_a_, gwf_diag_t *a, gwf_set64_t *h)
{
	int32_t i, x, n = *n_a_;
	kdq_t(gwf_diag_t) *A, *B;
	gwf_diag_t *b;

	*is_end = 0;
	for (i = 0, x = 1; i < 32; ++i, x <<= 1)
		if (x >= n) break;
	if (i < 4) i = 4;
	A = kdq_init2(gwf_diag_t, km, i); // $A is a queue. A stack also work, but a queue looks cleaner.
	for (i = 0; i < n; ++i) *kdq_pushp(gwf_diag_t, A) = a[i];
	kfree(km, a); // $a is not used as it has been copied to $A

	B = kdq_init2(gwf_diag_t, km, A->bits + 1); // $B is a stack.
	gwf_set64_clear(h);
	while (kdq_size(A)) {
		gwf_diag_t t = *kdq_shift(gwf_diag_t, A);
		uint32_t v = t.vd >> 32; // vertex
		int32_t d = (int64_t)((uint32_t)t.vd) - 80000000LL; // diagonal
		int32_t k = (int32_t)t.k; // wave front position
		int32_t i = k + d; // query position
		while (k + 1 < g->len[v] && i + 1 < ql && g->seq[v][k+1] == q[i+1])
			++i, ++k;
		if (k + 1 < g->len[v] && i + 1 < ql) {
			int32_t push1 = 1, push2 = 1;
			if (B->count >= 2) push1 = gwf_ed_update(&B->a[B->count - 2], v, d-1, k+1);
			if (B->count >= 1) push2 = gwf_ed_update(&B->a[B->count - 1], v, d,   k+1);
			if (push1) gwf_ed_push(B, v, d-1, k+1);
			if (push2) gwf_ed_push(B, v, d,   k+1);
			gwf_ed_push(B, v, d+1, k);
		} else if (i + 1 < ql) { // k + 1 == g->len[v]
			int32_t ov = g->aux[v]>>32, nv = (int32_t)g->aux[v], j, n_ext = 0;
			for (j = 0; j < nv; ++j) {
				uint32_t w = (uint32_t)g->arc[ov + j];
				int absent;
				khint_t itr;
				itr = gwf_set64_put(h, (uint64_t)w<<32 | (i + 1), &absent);
				if (q[i + 1] == g->seq[w][0]) {
					++n_ext;
					if (absent) gwf_ed_push(A, w, i+1, 0);
				} else if (absent) {
					gwf_ed_push(B, w, i,   0);
					gwf_ed_push(B, w, i+1, 0);
				}
			}
			if (nv == 0 || n_ext != nv)
				gwf_ed_push(B, v, d+1, k);
		} else if (k + 1 < g->len[v]) { // i + 1 == ql
			gwf_ed_push(B, v, d-1, k+1);
		} else if (v != v1) { // i + 1 == ql && k + 1 == g->len[v]
			int32_t ov = g->aux[v]>>32, nv = (int32_t)g->aux[v], j;
			for (j = 0; j < nv; ++j) {
				uint32_t w = (uint32_t)g->arc[ov + j];
				gwf_ed_push(B, w, i, 0);
			}
		} else {
			*is_end = 1, *n_a_ = 0;
			kdq_destroy(gwf_diag_t, A);
			kdq_destroy(gwf_diag_t, B);
			return 0;
		}
	}

	kdq_destroy(gwf_diag_t, A);
	*n_a_ = kdq_size(B), b = B->a, B->a = 0;
	kdq_destroy(gwf_diag_t, B);
	return b;
}

int32_t gwf_ed(void *km, const gwf_graph_t *g, int32_t ql, const char *q, int32_t v0, int32_t v1)
{
	int32_t s = 0, n_a = 1, is_end;
	gwf_diag_t *a;
	gwf_set64_t *h;
	h = gwf_set64_init();
	KCALLOC(km, a, 1);
	a[0].vd = (uint64_t)v0<<32 | 80000000LL, a[0].k = -1;
	while (n_a > 0) {
		a = gwf_ed_extend(km, g, ql, q, v1, &is_end, &n_a, a, h);
		if (g->n_vtx > 1 && ((s+1) & 0x7f) == 0) // only need to dedup for graphs
			n_a = gwf_ed_dedup(km, n_a, a);
		if (is_end || n_a == 0) break;
		++s;
#ifdef GWF_DEBUG
		printf("[%s] s=%d, n=%d\n", __func__, s, n_a);
#endif
	}
	gwf_set64_destroy(h);
	return is_end? s : -1;
}
