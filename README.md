span-analysis
=============

Analysis and comparison scripts for [SPAN](https://research.jetbrains.org/groups/biolabs/tools/span-peak-analyzer).

Benchmarks:

* `mcgill` - original approach for peak caller tuning
* `SPAN Comparison Opendnase.ipynb` - benchmark vs DNAse
* `SPAN H3K4me3 H3K36me3 vs RNA-Seq benchmark.ipynb` 
* `SPAN noise experiment.ipynb` & `noise` - SPAN analysis in different noise levels

Datasets:

* `GSE74310 - scATAC-Seq SPAN MACS2 vs DNAse.ipynb` - better sensitivity vs mice DNAse hypersensitivity.
* `GSE65360 - scATACSeq SPAN MACS2 vs DNAse.ipynb` - scATAC-Seq GSE65360 vs human DNAse hypersensitivity for 3 out of 5
  cell lines - better sensitivity vs DNAse
  (K562, H1-ESC, BJ, GM12878, HL-60) 3 out of 5.
* `GSE65360 - SPAN failed models analysis.ipynb` - SPAN failed models analysis
* `GSE112622 - transcription vs H3K36me3.ipynb` - ULI H3K36me3 vs RNA-Seq data with replicates - better consistency with
  K36me3.
* `ULI ABF Monocytes aging dataset` - better consistency with different signal-to-noise ratio.
  http://artyomovlab.wustl.edu/aging/
* `PRJN392905 Immgen - ATAC-Seq MACS vs SPAN automatic markup.ipynb`
  , `PRJN392905 Immgen - ATAC-Seq MACS2, SPAN replicates consistency.ipynb` - monocytes ATAC-Seq data with replicates -
  better consistency within replicates + automatic markup
* `GSE53643 - Replicated H3K4me2 consistency, SPAN automarkup.ipynb` - Replicated H3K4me2 consistency
* `GSE104284 - Replicated K27ac K4me1 K4me3 mice injury.ipynb`
* `GSE16256 - RoadmapEpigenomics.ipynb`
* `GSE26320 - ENCODE.ipynb`

ABF Paper reviews:

* `SPAN ABF Paper markup analysis.ipynb`
* `SPAN ABF Paper test error swarmplots.ipynb`
