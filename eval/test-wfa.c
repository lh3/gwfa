#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include "ketopt.h"
#include "kseq.h"
#include "gap_affine/affine_wavefront_align.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	gzFile fp1, fp2;
	kseq_t *ks1, *ks2;
	ketopt_t o = KETOPT_INIT;
	int c;

	while ((c = ketopt(&o, argc, argv, 1, "", 0)) >= 0) {
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: test-wfa <in1.fa> <in2.fa>\n");
		return 1;
	}

	fp1 = gzopen(argv[o.ind+0], "r");
	fp2 = gzopen(argv[o.ind+1], "r");
	assert(fp1 && fp2);
	ks1 = kseq_init(fp1);
	ks2 = kseq_init(fp2);
	kseq_read(ks1);
	kseq_read(ks2);

	mm_allocator_t* const mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
	affine_penalties_t pan = { .match = 0, .mismatch = 1, .gap_opening = 1, .gap_extension = 1 }; // Init Affine-WFA
	affine_wavefronts_t *wf =
		affine_wavefronts_new_complete(ks2->seq.l, ks1->seq.l, &pan, NULL, mm_allocator);
	affine_wavefronts_align(wf, ks2->seq.s, ks2->seq.l, ks1->seq.s, ks1->seq.l);
	const int score = edit_cigar_score_gap_affine(&wf->edit_cigar, &pan);
	printf("%s\t%s\t%d\n", ks1->name.s, ks2->name.s, score);

	kseq_destroy(ks1);
	kseq_destroy(ks2);
	gzclose(fp1);
	gzclose(fp2);
	return 0;
}
