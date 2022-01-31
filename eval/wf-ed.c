#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "kalloc.h"
#include "kvec.h"

typedef struct {
	int32_t d, k;
} wf_diag_t;

typedef kvec_t(wf_diag_t) wf_diag_v;

static int32_t wf_step(void *km, int32_t tl, const char *ts, int32_t ql, const char *qs, wf_diag_v *A, wf_diag_v *B)
{
	int32_t j, n = A->n;
	wf_diag_t *a = A->a, *b;

	// extend
#if 0
	for (j = 0; j < n; ++j) {
		wf_diag_t *p = &a[j];
		int32_t k = p->k;
		int32_t max_k = (ql - p->d < tl? ql - p->d : tl) - 1;
		const char *ts_ = ts + 1, *qs_ = qs + p->d + 1;
		while (k < max_k && *(ts_ + k) == *(qs_ + k))
			++k;
		if (k + p->d == ql - 1) return 1; // found semi glocal
		p->k = k;
	}
#else
	for (j = 0; j < n; ++j) {
		wf_diag_t *p = &a[j];
		int32_t k = p->k, i = k + p->d;
		while (k + 1 < tl && i + 1 < ql && ts[k+1] == qs[i+1])
			++k, ++i;
		if (i == ql - 1) return 1; // found semi glocal
		p->k = k;
	}
#endif

	// next
	kv_resize(wf_diag_t, km, *B, n + 2);
	b = B->a, B->n = 0;
	b[0].d = a[0].d - 1;
	b[0].k = a[0].k + 1;
	b[1].d = a[0].d;
	b[1].k = (n == 1 || a[0].k > a[1].k? a[0].k : a[1].k) + 1;
	for (j = 1; j < n - 1; ++j) {
		int32_t k;
		k = a[j-1].k;
		k = k > a[j+1].k + 1? k : a[j+1].k + 1;
		k = k > a[j].k + 1? k : a[j].k + 1;
		b[j+1].d = a[j].d, b[j+1].k = k;
	}
	if (n >= 2) {
		b[n].d = a[n-1].d;
		b[n].k = a[n-2].k > a[n-1].k + 1? a[n-2].k : a[n-1].k + 1;
	}
	b[n+1].d = a[n-1].d + 1;
	b[n+1].k = a[n-1].k;

	// out-of-bound cells
	kv_resize(wf_diag_t, km, *A, n + 2);
	A->n = 0, a = A->a;
	for (j = 0; j < n + 2; ++j) {
		wf_diag_t *p = &b[j];
		int32_t i = p->d + p->k;
		if (i >= ql || p->k >= tl || i < -1 || p->k < -1) continue;
		a[A->n++] = *p;
	}
	return 0;
}

int32_t wf_ed(void *km, int32_t tl, const char *ts, int32_t ql, const char *qs)
{
	int32_t s = 0;
	wf_diag_v A = {0,0,0}, B = {0,0,0};
	kv_resize(wf_diag_t, km, A, 16);
	A.a[A.n].d = 0, A.a[A.n++].k = -1;
	while (1) {
		int32_t ret = wf_step(km, tl, ts, ql, qs, &A, &B);
		if (ret) break;
		++s;
		printf("[%s] s=%d, n=%ld\n", __func__, s, A.n);
	}
	kfree(km, A.a);
	kfree(km, B.a);
	return s;
}

#ifdef WF_ED_MAIN
#include <zlib.h>
#include "ketopt.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	gzFile fp1, fp2;
	kseq_t *ks1, *ks2;
	ketopt_t o = KETOPT_INIT;
	int c, s;

	while ((c = ketopt(&o, argc, argv, 1, "", 0)) >= 0) {
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: test-wf-ed <in1.fa> <in2.fa>\n");
		return 1;
	}

	fp1 = gzopen(argv[o.ind+0], "r");
	fp2 = gzopen(argv[o.ind+1], "r");
	assert(fp1 && fp2);
	ks1 = kseq_init(fp1);
	ks2 = kseq_init(fp2);
	kseq_read(ks1);
	kseq_read(ks2);

	s = wf_ed(0, ks1->seq.l, ks1->seq.s, ks2->seq.l, ks2->seq.s);
	printf("%s\t%s\t%d\n", ks1->name.s, ks2->name.s, s);

	kseq_destroy(ks1);
	kseq_destroy(ks2);
	gzclose(fp1);
	gzclose(fp2);
	return 0;
}
#endif
