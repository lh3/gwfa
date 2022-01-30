#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include "ketopt.h"
#include "kseq.h"
#include "edlib.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	gzFile fp1, fp2;
	kseq_t *ks1, *ks2;
	ketopt_t o = KETOPT_INIT;
	int c;
	EdlibAlignResult rst;

	while ((c = ketopt(&o, argc, argv, 1, "", 0)) >= 0) {
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: test-edlib <in1.fa> <in2.fa>\n");
		return 1;
	}

	fp1 = gzopen(argv[o.ind+0], "r");
	fp2 = gzopen(argv[o.ind+1], "r");
	assert(fp1 && fp2);
	ks1 = kseq_init(fp1);
	ks2 = kseq_init(fp2);
	kseq_read(ks1);
	kseq_read(ks2);

	rst = edlibAlign(ks2->seq.s, ks2->seq.l, ks1->seq.s, ks1->seq.l,
			edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_DISTANCE, NULL, 0));
	printf("%s\t%s\t%d\n", ks1->name.s, ks2->name.s, rst.editDistance);

	kseq_destroy(ks1);
	kseq_destroy(ks2);
	gzclose(fp1);
	gzclose(fp2);
	return 0;
}
