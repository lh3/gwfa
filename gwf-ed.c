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

#define GWF_DIAG_SHIFT 0x40000000

static inline uint64_t gwf_gen_vd(uint32_t v, int32_t d)
{
	return (uint64_t)v<<32 | (GWF_DIAG_SHIFT + d);
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
		printf("Z\t%d\t%d\t%d\n", (int32_t)(a[i].vd0>>32), (int32_t)a[i].vd0 - GWF_DIAG_SHIFT, (int32_t)a[i].vd1 - GWF_DIAG_SHIFT);
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
	uint64_t vd; // higher 32 bits: vertex ID; lower 32 bits: diagonal+0x4000000
	int32_t k;
	uint32_t o;
} gwf_diag_t;

typedef kvec_t(gwf_diag_t) gwf_diag_v;

#define ed_key(x) ((x).vd)
KRADIX_SORT_INIT(gwf_ed, gwf_diag_t, ed_key, 8)

KDQ_INIT(gwf_diag_t)

void gwf_ed_print_diag(size_t n, gwf_diag_t *a) // for debugging only
{
	size_t i;
	for (i = 0; i < n; ++i) {
		int32_t d = (int32_t)a[i].vd - GWF_DIAG_SHIFT;
		printf("Z\t%d\t%d\t%d\t%d\t%d\n", (int32_t)(a[i].vd>>32), d, a[i].k, d + a[i].k, a[i].o);
	}
}

// push (v,d,k) to the end of the queue
static inline void gwf_diag_push(void *km, gwf_diag_v *a, uint32_t v, int32_t d, int32_t k, uint32_t ooo)
{
	gwf_diag_t *p;
	kv_pushp(gwf_diag_t, km, *a, &p);
	p->vd = gwf_gen_vd(v, d), p->k = k, p->o = ooo;
}

// determine the wavefront on diagonal (v,d)
static inline int32_t gwf_diag_update(gwf_diag_t *p, uint32_t v, int32_t d, int32_t k, uint32_t ooo)
{
	uint64_t vd = gwf_gen_vd(v, d);
	if (p->vd == vd) {
		p->o = p->k > k? p->o : ooo;
		p->k = p->k > k? p->k : k;
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
	int32_t i, j, k, n_b, n_c;
	gwf_diag_t *b, *c;

	kv_resize(gwf_diag_t, km, *ooo, n_a);
	for (i = 0, n_c = 0; i < n_a; ++i)
		if (a[i].o) ++n_c;
	n_b = n_a - n_c;
	b = ooo->a, c = b + n_b;
	for (i = j = k = 0; i < n_a; ++i) {
		if (a[i].o) c[k++] = a[i];
		else b[j++] = a[i];
	}
	radix_sort_gwf_ed(c, c + n_c);
	for (k = 0; k < n_c; ++k) c[k].o = 0;

	i = j = k = 0;
	while (i < n_b && j < n_c) {
		if (b[i].vd <= c[j].vd)
			a[k++] = b[i++];
		else a[k++] = c[j++];
	}
	while (i < n_b) a[k++] = b[i++];
	while (j < n_c) a[k++] = c[j++];
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
	gwf_set64_t *ha; // hash table for adjacency
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

// reach the wavefront
static inline int32_t gwf_extend1(int32_t d, int32_t k, int32_t vl, const char *ts, int32_t ql, const char *qs)
{
	int32_t max_k = (ql - d < vl? ql - d : vl) - 1;
	const char *ts_ = ts + 1, *qs_ = qs + d + 1;
#if 0
	// int32_t i = k + d; while (k + 1 < g->len[v] && i + 1 < ql && g->seq[v][k+1] == q[i+1]) ++k, ++i;
	while (k < max_k && *(ts_ + k) == *(qs_ + k))
		++k;
#else
	uint64_t cmp = 0;
	while (k + 7 < max_k) {
		uint64_t x = *(uint64_t*)(ts_ + k); // warning: unaligned memory access
		uint64_t y = *(uint64_t*)(qs_ + k);
		cmp = x ^ y;
		if (cmp == 0) k += 8;
		else break;
	}
	if (cmp)
		k += __builtin_ctzl(cmp) >> 3; // on x86, this is done via the BSR instruction: https://www.felixcloutier.com/x86/bsr
	else if (k + 7 >= max_k)
		while (k < max_k && *(ts_ + k) == *(qs_ + k)) // use this for generic CPUs. It is slightly faster than the unoptimized version
			++k;
#endif
	return k;
}

// This is essentially Landau-Vishkin for linear sequences. The function speeds up alignment to long vertices. Not really necessary.
static void gwf_ed_extend_batch(void *km, const gwf_graph_t *g, int32_t ql, const char *q, int32_t n, gwf_diag_t *a, gwf_diag_v *B,
								kdq_t(gwf_diag_t) *A, gwf_intv_v *tmp_intv)
{
	int32_t j, m;
	int32_t v = a->vd>>32;
	int32_t vl = g->len[v];
	const char *ts = g->seq[v];
	gwf_diag_t *b;

	// wfa_extend
	for (j = 0; j < n; ++j) {
		int32_t k;
		k = gwf_extend1((int32_t)a[j].vd - GWF_DIAG_SHIFT, a[j].k, vl, ts, ql, q);
		a[j].k = k;
	}

	// wfa_next
	kv_resize(gwf_diag_t, km, *B, B->n + n + 2);
	b = &B->a[B->n];
	b[0].vd = a[0].vd - 1;
	b[0].k = a[0].k + 1;
	b[0].o = 0;
	b[1].vd = a[0].vd;
	b[1].k  = (n == 1 || a[0].k > a[1].k? a[0].k : a[1].k) + 1;
	b[1].o  = 0;
	for (j = 1; j < n - 1; ++j) {
		int32_t k = a[j-1].k;
		k = k > a[j].k + 1? k : a[j].k + 1;
		k = k > a[j+1].k + 1? k : a[j+1].k + 1;
		b[j+1].vd = a[j].vd, b[j+1].k = k, b[j+1].o = 0;
	}
	if (n >= 2) {
		b[n].vd = a[n-1].vd;
		b[n].k  = a[n-2].k > a[n-1].k + 1? a[n-2].k : a[n-1].k + 1;
		b[n].o  = 0;
	}
	b[n+1].vd = a[n-1].vd + 1;
	b[n+1].k  = a[n-1].k;
	b[n+1].o  = 0;

	// drop out-of-bound cells
	for (j = 0; j < n; ++j) {
		gwf_diag_t *p = &a[j];
		if (p->k == vl - 1 || (int32_t)p->vd - GWF_DIAG_SHIFT + p->k == ql - 1)
			p->o |= 1, *kdq_pushp(gwf_diag_t, A) = *p;
	}
	for (j = 0, m = 0; j < n + 2; ++j) {
		gwf_diag_t *p = &b[j];
		int32_t d = (int32_t)p->vd - GWF_DIAG_SHIFT;
		if (d + p->k < ql && p->k < vl) {
			b[m++] = *p;
		} else if (p->k == vl) {
			gwf_intv_t *q;
			kv_pushp(gwf_intv_t, km, *tmp_intv, &q);
			q->vd0 = gwf_gen_vd(v, d), q->vd1 = q->vd0 + 1;
		}
	}
	B->n += m;
}

// wfa_extend and wfa_next combined
static gwf_diag_t *gwf_ed_extend(gwf_edbuf_t *buf, const gwf_graph_t *g, int32_t ql, const char *q, int32_t v1,
								 int32_t *end_v, int32_t *end_off, int32_t *n_a_, gwf_diag_t *a)
{
	int32_t i, x, n = *n_a_, do_dedup = 1;
	kdq_t(gwf_diag_t) *A;
	gwf_diag_v B = {0,0,0};
	gwf_diag_t *b;

	*end_v = *end_off = -1;
	buf->tmp.n = 0;
	gwf_set64_clear(buf->ha); // hash table $h to avoid visiting a vertex twice
	for (i = 0, x = 1; i < 32; ++i, x <<= 1)
		if (x >= n) break;
	if (i < 4) i = 4;
	A = kdq_init2(gwf_diag_t, buf->km, i); // $A is a queue
	kv_resize(gwf_diag_t, buf->km, B, n * 2);
#if 0 // unoptimized version without calling gwf_ed_extend_batch() at all. The final result will be the same.
	A->count = n;
	memcpy(A->a, a, n * sizeof(*a));
#else // optimized for long vertices.
	for (x = 0, i = 1; i <= n; ++i) {
		if (i == n || a[i].vd != a[i-1].vd + 1) {
			gwf_ed_extend_batch(buf->km, g, ql, q, i - x, &a[x], &B, A, &buf->tmp);
			x = i;
		}
	}
	if (kdq_size(A) == 0) do_dedup = 0;
#endif
	kfree(buf->km, a); // $a is not used as it has been copied to $A

	while (kdq_size(A)) {
		gwf_diag_t t;
		int32_t ooo, v, d, k, i, vl;

		t = *kdq_shift(gwf_diag_t, A);
		ooo = t.o, v = t.vd >> 32; // vertex
		d = (int32_t)t.vd - GWF_DIAG_SHIFT; // diagonal
		k = t.k; // wavefront position on the vertex
		vl = g->len[v]; // $vl is the vertex length
		k = gwf_extend1(d, k, vl, g->seq[v], ql, q);
		i = k + d; // query position

		if (k + 1 < vl && i + 1 < ql) { // the most common case: the wavefront is in the middle
			int32_t push1 = 1, push2 = 1;
			if (B.n >= 2) push1 = gwf_diag_update(&B.a[B.n - 2], v, d-1, k+1, ooo);
			if (B.n >= 1) push2 = gwf_diag_update(&B.a[B.n - 1], v, d,   k+1, ooo);
			if (push1)          gwf_diag_push(buf->km, &B, v, d-1, k+1, 1);
			if (push2 || push1) gwf_diag_push(buf->km, &B, v, d,   k+1, 1);
			gwf_diag_push(buf->km, &B, v, d+1, k, ooo);
		} else if (i + 1 < ql) { // k + 1 == g->len[v]; reaching the end of the vertex but not the end of query
			int32_t ov = g->aux[v]>>32, nv = (int32_t)g->aux[v], j, n_ext = 0;
			gwf_intv_t *p;
			kv_pushp(gwf_intv_t, buf->km, buf->tmp, &p);
			p->vd0 = gwf_gen_vd(v, d), p->vd1 = p->vd0 + 1;
			for (j = 0; j < nv; ++j) { // traverse $v's neighbors
				uint32_t w = (uint32_t)g->arc[ov + j]; // $w is next to $v
				int absent;
				gwf_set64_put(buf->ha, (uint64_t)w<<32 | (i + 1), &absent); // test if ($w,$i) has been visited
				if (q[i + 1] == g->seq[w][0]) { // can be extended to the next vertex without a mismatch
					++n_ext;
					if (absent) {
						gwf_diag_t *p;
						p = kdq_pushp(gwf_diag_t, A);
						p->vd = gwf_gen_vd(w, i+1), p->k = 0, p->o = 1;
					}
				} else if (absent) {
					gwf_diag_push(buf->km, &B, w, i,   0, 1);
					gwf_diag_push(buf->km, &B, w, i+1, 0, 1);
				}
			}
			if (nv == 0 || n_ext != nv) // add an insertion to the target; this *might* cause a duplicate in corner cases
				gwf_diag_push(buf->km, &B, v, d+1, k, 1);
		} else if (v1 < 0 || (v == v1 && k + 1 == vl)) { // i + 1 == ql
			*end_v = v, *end_off = k, *n_a_ = 0;
			kdq_destroy(gwf_diag_t, A);
			kfree(buf->km, B.a);
			return 0;
		} else if (k + 1 < vl) { // i + 1 == ql; reaching the end of the query but not the end of the vertex
			gwf_diag_push(buf->km, &B, v, d-1, k+1, ooo); // add an deletion; this *might* case a duplicate in corner cases
		} else if (v != v1) { // i + 1 == ql && k + 1 == g->len[v]; not reaching the last vertex $v1
			int32_t ov = g->aux[v]>>32, nv = (int32_t)g->aux[v], j;
			for (j = 0; j < nv; ++j) {
				uint32_t w = (uint32_t)g->arc[ov + j];
				gwf_diag_push(buf->km, &B, w, i, 0, 1); // deleting the first base on the next vertex
			}
		} else assert(0); // should never come here
	}

	kdq_destroy(gwf_diag_t, A);
	*n_a_ = n = B.n, b = B.a;

	if (do_dedup) *n_a_ = n = gwf_dedup(buf, n, b);
	return b;
}

int32_t gwf_ed(void *km, const gwf_graph_t *g, int32_t ql, const char *q, int32_t v0, int32_t v1, gwf_path_t *path)
{
	int32_t s = 0, n_a = 1;
	gwf_diag_t *a;
	gwf_edbuf_t buf;

	memset(&buf, 0, sizeof(buf));
	buf.km = km;
	buf.ha = gwf_set64_init2(km);
	KCALLOC(km, a, 1);
	a[0].vd = gwf_gen_vd(v0, 0), a[0].k = -1, a[0].o = 0; // the initial state
	while (n_a > 0) {
		a = gwf_ed_extend(&buf, g, ql, q, v1, &path->end_v, &path->end_off, &n_a, a);
		if (path->end_off >= 0 || n_a == 0) break;
		++s;
#ifdef GWF_DEBUG
		printf("[%s] s=%d, n=%d\n", __func__, s, n_a);
#endif
	}
	gwf_set64_destroy(buf.ha);
	kfree(km, buf.intv.a); kfree(km, buf.tmp.a); kfree(km, buf.swap.a);
	path->s = path->end_v >= 0? s : -1;
	return path->s; // end_v < 0 could happen if v0 can't reach v1
}
