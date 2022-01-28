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

// index the graph such that we can quickly access the neighbors of a vertex
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

// free the index
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
#include "kvec.h"

static inline uint64_t gwf_gen_vd(uint32_t v, int32_t d)
{
	return (uint64_t)v<<32 | (80000000LL + d);
}

/*
 * Diagonal interval
 */
typedef struct {
	uint64_t vd0, vd1;
} gwf_intv_t;

typedef kvec_t(gwf_intv_t) gwf_intv_v;

#define intvd_key(x) ((x).vd0)
KRADIX_SORT_INIT(gwf_intv, gwf_intv_t, intvd_key, 8)

static int gwf_intv_is_sorted(int32_t n_a, const gwf_intv_t *a)
{
	int32_t i;
	for (i = 1; i < n_a; ++i)
		if (a[i-1].vd0 > a[i].vd0) break;
	return (i == n_a);
}

void gwf_ed_print_intv(size_t n, gwf_intv_t *a) // for debugging only
{
	size_t i;
	for (i = 0; i < n; ++i)
		printf("Z\t%d\t%d\t%d\n", (int32_t)(a[i].vd0>>32), (int32_t)a[i].vd0 - 80000000, (int32_t)a[i].vd1 - 80000000);
}

// merge overlapping intervals; input must be sorted
static size_t gwf_intv_merge_adj(size_t n, gwf_intv_t *a)
{
	size_t i, k;
	uint64_t st, en;
	if (n == 0) return 0;
	st = a[0].vd0, en = a[0].vd1;
	for (i = 1, k = 0; i < n; ++i) {
		if (a[i].vd0 > en) {
			a[k].vd0 = st, a[k++].vd1 = en;
			st = a[i].vd0, en = a[i].vd1;
		} else en = en > a[i].vd1? en : a[i].vd1;
	}
	a[k].vd0 = st, a[k++].vd1 = en;
	return k;
}

// merge two sorted interval lists
static size_t gwf_intv_merge2(gwf_intv_t *a, size_t n_b, const gwf_intv_t *b, size_t n_c, const gwf_intv_t *c)
{
	size_t i = 0, j = 0, k = 0;
	while (i < n_b && j < n_c) {
		if (b[i].vd0 <= c[j].vd0)
			a[k++] = b[i++];
		else a[k++] = c[j++];
	}
	while (i < n_b) a[k++] = b[i++];
	while (j < n_c) a[k++] = c[j++];
	return gwf_intv_merge_adj(k, a);
}

/*
 * Diagonal
 */
typedef struct { // a diagonal
	uint64_t vd; // NB: this wastes 4 bytes due to memory alignment
	int32_t k, ooo; // ooo = out of order
} gwf_diag_t;

typedef kvec_t(gwf_diag_t) gwf_diag_v;

#define ed_key(x) ((x).vd)
KRADIX_SORT_INIT(gwf_ed, gwf_diag_t, ed_key, 8)

KDQ_INIT(gwf_diag_t)

// push (v,d,k) to the end of the queue
static inline void gwf_diag_push(kdq_t(gwf_diag_t) *a, uint32_t v, int32_t d, int32_t k, int32_t ooo)
{
	gwf_diag_t t;
	t.vd = gwf_gen_vd(v, d), t.k = k, t.ooo = ooo;
	*kdq_pushp(gwf_diag_t, a) = t;
}

// determine the wavefront on diagonal (v,d)
static inline int32_t gwf_diag_update(gwf_diag_t *p, uint32_t v, int32_t d, int32_t k, int32_t ooo)
{
	uint64_t vd = gwf_gen_vd(v, d);
	if (p->vd == vd) {
		p->k = p->k > k? p->k : k, p->ooo = ooo;
		return 0;
	}
	return 1;
}

static int gwf_diag_is_sorted(int32_t n_a, const gwf_diag_t *a)
{
	int32_t i;
	for (i = 1; i < n_a; ++i)
		if (a[i-1].vd > a[i].vd) break;
	return (i == n_a);
}

// sort a[]. This uses the gwf_diag_t::ooo field to speed up sorting.
static void gwf_diag_sort(int32_t n_a, gwf_diag_t *a, void *km, gwf_diag_v *ooo)
{
#if 1
	int32_t i, j, k, n_b, n_c;
	gwf_diag_t *b, *c;

	kv_resize(gwf_diag_t, km, *ooo, n_a);
	for (i = 0, n_c = 0; i < n_a; ++i)
		if (a[i].ooo) ++n_c;
	n_b = n_a - n_c;
	b = ooo->a, c = b + n_b;
	for (i = j = k = 0; i < n_a; ++i) {
		if (a[i].ooo) c[k++] = a[i];
		else b[j++] = a[i];
	}
	radix_sort_gwf_ed(c, c + n_c);

	i = j = k = 0;
	while (i < n_b && j < n_c) {
		if (b[i].vd <= c[j].vd)
			a[k++] = b[i++];
		else a[k++] = c[j++];
	}
	while (i < n_b) a[k++] = b[i++];
	while (j < n_c) a[k++] = c[j++];
	for (i = 0; i < n_a; ++i) a[i].ooo = 0;
#else
	radix_sort_gwf_ed(a, a + n_a); // the whole function could just call this line but this would be much slower.
#endif
}

// remove diagonals not on the wavefront
static int32_t gwf_diag_dedup(int32_t n_a, gwf_diag_t *a, void *km, gwf_diag_v *ooo)
{
	int32_t i, n, st;
	if (!gwf_diag_is_sorted(n_a, a))
		gwf_diag_sort(n_a, a, km, ooo);
	for (i = 1, st = 0, n = 0; i <= n_a; ++i) {
		if (i == n_a || a[i].vd != a[st].vd) {
			int32_t j, max_j = st;
			if (st + 1 < i)
				for (j = st + 1; j < i; ++j) // choose the far end (i.e. the wavefront)
					if (a[max_j].k < a[j].k) max_j = j;
			a[n++] = a[max_j];
			st = i;
		}
	}
	return n;
}

// use forbidden bands to remove diagonals not on the wavefront
static int32_t gwf_mixed_dedup(int32_t n_a, gwf_diag_t *a, int32_t n_b, gwf_intv_t *b)
{
	int32_t i = 0, j = 0, k = 0;
	while (i < n_a && j < n_b) {
		//printf("Z\t%d\t%d\t%d\t%d\t%d\n", (uint32_t)(a[i].vd>>32), (int32_t)a[i].vd - 80000000, (uint32_t)(b[j].vd0>>32), (int32_t)b[j].vd0 - 80000000, (int32_t)b[j].vd1 - 80000000);
		if (a[i].vd >= b[j].vd0 && a[i].vd < b[j].vd1) ++i;
		else if (a[i].vd >= b[j].vd1) ++j;
		else a[k++] = a[i++];
	}
	while (i < n_a) a[k++] = a[i++];
	return k;
}

/*
 * Core GWFA routine
 */
KHASHL_INIT(KH_LOCAL, gwf_set64_t, gwf_set64, uint64_t, kh_hash_dummy, kh_eq_generic)

typedef struct {
	void *km;
	gwf_set64_t *h;
	gwf_intv_v intv;
	gwf_intv_v tmp, swap;
	gwf_diag_v ooo;
} gwf_edbuf_t;

// remove diagonals not on the wavefront
static int32_t gwf_dedup(gwf_edbuf_t *buf, int32_t n_a, gwf_diag_t *a)
{
	if (buf->intv.n + buf->tmp.n > 0) {
		if (!gwf_intv_is_sorted(buf->tmp.n, buf->tmp.a))
			radix_sort_gwf_intv(buf->tmp.a, buf->tmp.a + buf->tmp.n);
		kv_copy(gwf_intv_t, buf->km, buf->swap, buf->intv);
		kv_resize(gwf_intv_t, buf->km, buf->intv, buf->intv.n + buf->tmp.n);
		buf->intv.n = gwf_intv_merge2(buf->intv.a, buf->swap.n, buf->swap.a, buf->tmp.n, buf->tmp.a);
	}
	int32_t n0 = n_a, n1;
	n_a = gwf_diag_dedup(n_a, a, buf->km, &buf->ooo);
	n1 = n_a;
	if (buf->intv.n > 0)
		n_a = gwf_mixed_dedup(n_a, a, buf->intv.n, buf->intv.a);
#ifdef GWF_DEBUG
	printf("[%s] intv.n=%ld; dedup: %d -> %d -> %d\n", __func__, buf->intv.n, n0, n1, n_a);
#endif
	return n_a;
}

static gwf_diag_t *gwf_ed_extend(gwf_edbuf_t *buf, const gwf_graph_t *g, int32_t ql, const char *q, int32_t v1, int32_t *end_v, int32_t *end_off, int32_t *n_a_, gwf_diag_t *a)
{
	int32_t i, x, n = *n_a_;
	kdq_t(gwf_diag_t) *A, *B;
	gwf_diag_t *b;

	*end_v = -1, *end_off = -1;
	for (i = 0, x = 1; i < 32; ++i, x <<= 1)
		if (x >= n) break;
	if (i < 4) i = 4;
	A = kdq_init2(gwf_diag_t, buf->km, i); // $A is a queue
	for (i = 0; i < n; ++i) *kdq_pushp(gwf_diag_t, A) = a[i]; // copy $a to $A; only necessary for graphs
	kfree(buf->km, a); // $a is not used as it has been copied to $A

	B = kdq_init2(gwf_diag_t, buf->km, A->bits + 1); // $B is a stack; could be replaced with a simple vector
	gwf_set64_clear(buf->h); // hash table $h to avoid visiting a vertex twice
	buf->tmp.n = 0;
	while (kdq_size(A)) {
		gwf_diag_t t = *kdq_shift(gwf_diag_t, A);
		uint32_t v = t.vd >> 32; // vertex
		int32_t d = (int64_t)((uint32_t)t.vd) - 80000000LL; // diagonal
		int32_t k = (int32_t)t.k; // wave front position
		int32_t i = k + d; // query position
		while (k + 1 < g->len[v] && i + 1 < ql && g->seq[v][k+1] == q[i+1]) // extend the diagonal $d to the wavefront
			++i, ++k;
		if (k + 1 < g->len[v] && i + 1 < ql) { // the most common case: the wavefront is in the middle
			int32_t push1 = 1, push2 = 1;
			if (B->count >= 2) push1 = gwf_diag_update(&B->a[B->count - 2], v, d-1, k+1, t.ooo);
			if (B->count >= 1) push2 = gwf_diag_update(&B->a[B->count - 1], v, d,   k+1, t.ooo);
			if (push1) gwf_diag_push(B, v, d-1, k+1, 1);
			if (push2 || push1) gwf_diag_push(B, v, d, k+1, 1);
			gwf_diag_push(B, v, d+1, k, t.ooo);
		} else if (i + 1 < ql) { // k + 1 == g->len[v]; reaching the end of the vertex but not the end of query
			int32_t ov = g->aux[v]>>32, nv = (int32_t)g->aux[v], j, n_ext = 0;
			gwf_intv_t *p;
			kv_pushp(gwf_intv_t, buf->km, buf->tmp, &p);
			p->vd0 = gwf_gen_vd(v, d), p->vd1 = p->vd0 + 1;
			for (j = 0; j < nv; ++j) { // traverse $v's neighbors
				uint32_t w = (uint32_t)g->arc[ov + j]; // $w is next to $v
				int absent;
				gwf_set64_put(buf->h, (uint64_t)w<<32 | (i + 1), &absent); // test if ($w,$i) has been visited
				if (q[i + 1] == g->seq[w][0]) { // can be extended to the next vertex without a mismatch
					++n_ext;
					if (absent) gwf_diag_push(A, w, i+1, 0, 1);
				} else if (absent) {
					gwf_diag_push(B, w, i,   0, 1);
					gwf_diag_push(B, w, i+1, 0, 1);
				}
			}
			if (nv == 0 || n_ext != nv) // add an insertion to the target; this *might* cause a duplicate in corner cases
				gwf_diag_push(B, v, d+1, k, 1);
		} else if (v1 < 0 || (v == v1 && k + 1 == g->len[v])) { // i + 1 == ql
			*end_v = v, *end_off = k, *n_a_ = 0;
			kdq_destroy(gwf_diag_t, A);
			kdq_destroy(gwf_diag_t, B);
			return 0;
		} else if (k + 1 < g->len[v]) { // i + 1 == ql; reaching the end of the query but not the end of the vertex
			gwf_diag_push(B, v, d-1, k+1, t.ooo); // add an deletion; this *might* case a duplicate in corner cases
		} else if (v != v1) { // i + 1 == ql && k + 1 == g->len[v]; not reaching the last vertex $v1
			int32_t ov = g->aux[v]>>32, nv = (int32_t)g->aux[v], j;
			for (j = 0; j < nv; ++j) {
				uint32_t w = (uint32_t)g->arc[ov + j];
				gwf_diag_push(B, w, i, 0, 1); // deleting the first base on the next vertex
			}
		} else {
			assert(0); // should never come here
		}
	}

	kdq_destroy(gwf_diag_t, A);
	n = kdq_size(B), b = B->a, B->a = 0;
	kdq_destroy(gwf_diag_t, B);

	*n_a_ = n = gwf_dedup(buf, n, b);
	return b;
}

int32_t gwf_ed(void *km, const gwf_graph_t *g, int32_t ql, const char *q, int32_t v0, int32_t v1)
{
	int32_t s = 0, n_a = 1, end_v, end_off;
	gwf_diag_t *a;
	gwf_edbuf_t buf;

	memset(&buf, 0, sizeof(buf));
	buf.km = km;
	buf.h = gwf_set64_init();
	KCALLOC(km, a, 1);
	a[0].vd = gwf_gen_vd(v0, 0), a[0].k = -1, a[0].ooo = 0; // the initial state
	while (n_a > 0) {
		a = gwf_ed_extend(&buf, g, ql, q, v1, &end_v, &end_off, &n_a, a);
		if (end_off >= 0 || n_a == 0) break;
		++s;
#ifdef GWF_DEBUG
		printf("[%s] s=%d, n=%d\n", __func__, s, n_a);
#endif
	}
	gwf_set64_destroy(buf.h);
	kfree(km, buf.intv.a); kfree(km, buf.tmp.a); kfree(km, buf.swap.a);
	return end_v >= 0? s : -1; // end_v < 0 could happen if v0 can't reach v1
}
