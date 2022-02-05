## Getting Started
```sh
git clone https://github.com/lh3/gwfa
cd gwfa && make
./gwf-test test/C4-90.gfa.gz test/C4-NA19240.1.fa.gz
./gwf-test test/C4-NA19240.1.fa.gz test/C4-NA19240.2.fa.gz
```

## Introduction

GWFA (Graph WaveFront Alignment) is an algorithm to align a sequence against a
sequence graph. It adapts the [WFA algorithm][wfa] for graphs. This repo
presents a _proof-of-concept_ implementation of GWFA that computes the edit
distance between a graph and a sequence without backtracing. The algorithm
assumes the start of the sequence to be aligned with the start of the first
segment in the graph and requires the query sequence to be fully aligned. This
behavior is similar to the SHW mode of [edlib][edlib]. It is not intended for
mapping reads against a whole-genome graph. This is not an enduser tool.

GWFA is optimized for graphs consisting of long segments. It is largely reduced
to [my implementation][mylv89] of the [Landau-Vishkin algorithm][lv89] given
two linear sequences as input. Similar to WFA, GWFA is fast when the edit
distance is small. It can align a ~120kb sequence to a ~160kb graph at 0.1%
divergence in 0.02 second, much faster than the ordinary Needleman-Wunsch
formulation.

## Evaluation

To evaluate the performance of GWFA, we constructed two small graphs with
minigraph. The first graph includes ~120kb region around the C4A and C4B genes.
The second includes ~5Mb of HLA, consisting of 980 segments and 1399 links.
Both graphs are available [via Zenodo][zenodo]. Neither has cycles.

For the C4-90 graph, [Haowen Zhang][haowen] found GWFA can report the same edit
distance in comparison to [his implementations][hz-sga] of various
sequence-to-graph algorithms.

For the HLA-57 graph, GWFA took 3m38s to align the HG002.2 haplotype (extracted
from HLA-61.agc on Zenodo) with 41932 edits. Other exact solutions did not
finish in a day. [GraphAligner][graphaligner]-1.0.13 failed to report an
alignment given the entire ~5Mb query sequence. For the first 2Mb subsequence
of HG002.2, GraphAligner reported 7075 edits in the vg mode and 6973 edits in
the dbg mode. GWFA found a smaller edit distance of 6447.

We have not tested more complex graphs with cycles. Although we think the GWFA
algorithm should be correct in theory, we are not sure if the current
implmentation is correct in all corner cases. Please use with caution.

## Contributors

[Shiqi Wu][shiqi] formulated the initial GWFA algorithm. [Haowen Zhang][haowen]
inspected the code and helped the evaluation.

[mylv89]: https://github.com/lh3/lv89
[lv89]: https://doi.org/10.1016/0196-6774(89)90010-2
[wfa]: https://github.com/smarco/WFA
[zenodo]: https://zenodo.org/record/5976063
[haowen]: https://github.com/haowenz
[hz-sga]: https://github.com/haowenz/SGA
[graphaligner]: https://github.com/maickrau/GraphAligner
[shiqi]: https://github.com/Shiqi-Wu
[edlib]: https://github.com/Martinsos/edlib
