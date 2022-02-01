#include <stdint.h>
#include <stdio.h>

typedef struct {
	int32_t d, k;
} wf_diag_t;

static int32_t wf_step(int32_t tl, const char *ts, int32_t ql, const char *qs, int32_t n, wf_diag_t *a)
{
	int32_t j, m;
	wf_diag_t *b = a + n + 2; // temporary array

	// wfa_extend
#if 0 // unoptimized original version
	for (j = 0; j < n; ++j) {
		wf_diag_t *p = &a[j];
		int32_t k = p->k, i = k + p->d;
		while (k + 1 < tl && i + 1 < ql && ts[k+1] == qs[i+1])
			++k, ++i;
		if (i == ql - 1) return 1; // found semi glocal
		p->k = k;
	}
#else // optimized version learned from WFA
	for (j = 0; j < n; ++j) {
		wf_diag_t *p = &a[j];
		int32_t k = p->k;
		int32_t max_k = (ql - p->d < tl? ql - p->d : tl) - 1;
		const char *ts_ = ts + 1, *qs_ = qs + p->d + 1;
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
		if (k + p->d == ql - 1) return -1; // found semi glocal
		p->k = k;
	}
#endif

	// wfa_next
	b[0].d = a[0].d - 1;
	b[0].k = a[0].k + 1;
	b[1].d = a[0].d;
	b[1].k = (n == 1 || a[0].k > a[1].k? a[0].k : a[1].k) + 1;
	for (j = 1; j < n - 1; ++j) {
		int32_t k = a[j-1].k;
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

	// drop out-of-bound cells
	for (j = 0, m = 0; j < n + 2; ++j)
		if (b[j].d + b[j].k < ql && b[j].k < tl)
			a[m++] = b[j];
	return m;
}

// mem should be at least (tl+ql)*16 long
int32_t wf_ed(int32_t tl, const char *ts, int32_t ql, const char *qs, uint8_t *mem)
{
	int32_t s = 0, n = 1;
	wf_diag_t *a;
	a = (wf_diag_t*)mem;
	a[0].d = 0, a[0].k = -1;
	while (1) {
		n = wf_step(tl, ts, ql, qs, n, a);
		if (n < 0) break;
		++s;
	}
	return s;
}

#ifdef WF_ED_MAIN
#include <assert.h>
#include <zlib.h>
#include <stdlib.h>
#include "ketopt.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	gzFile fp1, fp2;
	kseq_t *ks1, *ks2;
	ketopt_t o = KETOPT_INIT;
	int c, s;
	uint8_t *mem;

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

	mem = (uint8_t*)malloc((ks1->seq.l + ks2->seq.l) * 16);
	s = wf_ed(ks1->seq.l, ks1->seq.s, ks2->seq.l, ks2->seq.s, mem);
	free(mem);
	printf("%s\t%s\t%d\n", ks1->name.s, ks2->name.s, s);

	kseq_destroy(ks1);
	kseq_destroy(ks2);
	gzclose(fp1);
	gzclose(fp2);
	return 0;
}
#endif
