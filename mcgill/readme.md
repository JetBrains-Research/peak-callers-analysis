### SPAN benchmarking

We were inspired by the article ["Optimizing ChIP-seq peak detectors using visual labels and supervised machine learning", Bioinformatics 33(4), 2017, T. D. Hocking et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5408812/). In this work, among other things, the authors manually labeled several small regions of some ChIP-seq tracks using a visual inspection in a genome browser. They then determined how well various peak callers' output corresponded with those labels.

We've decided to apply this approach to SPAN. For the sake of reproducibility, we've created a set of scripts that generate the benchmarking results. We've also included the end results that are normally obtained on the steps 8 and 9.

## Usage:

1. Make sure you have about **120G+** free disk space; the required BED files are pretty big. Make sure that `R` is installed and upgraded to the latest version.

2. Download the scripts (`load.sh`, `call_peaks.sh`, `benchmark.R`, `benchmark_other.R`, `4fold.R`, `4fold_other.R`) and SPAN jar executive (e.g. `span-0.4.0.jar`) in the same folder.

3. Launch
```bash
$ load.sh
```
It will download the necessary data. The process can take **several days;** if the script fails to complete for any reason, just relaunch it; it should resume downloading right where it stopped.

The downloaded BED files will be in `chip-seq-benchmark` subdirectory.

4. Launch
```bash
$ call_peaks.sh
```
If necessary, you can supply FDR control value and GAP parameter as arguments, e.g. 
```bash
$ call_peaks.sh 1E-6 4
```
By default, FDR is controlled at `1E-4`, and `GAP` is `5`. The first launch of this script will likely take significant time; the consequent launches with different FDR control values will be much faster due to caching.

To launch the script for several FDR control and GAP values, you can use the following syntax:
```bash
$ for fdr in 1E-3 1E-4 1E-5; do for gap in 3 4 5; do call_peaks.sh $fdr $gap; done; done;
```

The script will generate multiple `peaks<sample>_<gap>_<fdr>.bed` peak files in the same folder where the corresponding BED signal files are.

5. Launch
```bash
$ Rscript benchmark.R
```
If necessary, you can supply FDR control and GAP values as arguments in comma-separated lists, e.g. 
```bash
Rscript benchmark.R 1E-4,1E-5,1E-6 3,4,5
```
The values should be exactly the same as the ones previously passed to `call_peaks.sh`. By default, they're `1E-4` and `5`.

If the script fails, complaining at the missing peak files, it's possible that the previous step failed to complete fully. Just launch `call_peaks.sh` again, it will skip the existing peak files and generate the missing ones only.

The script will generate multiple `<dataset>_<gap>_<fdr>.csv` and `summary_<gap>_<fdr>.csv` files in the working directory. They are tab-separated unquoted tables.

Dataset files:

Each row corresponds to a labelled genomic region. `chrom`, `chromStart` and `chromEnd` define the region in question, and `annotation` provides the label. `status` describes whether the label was confirmed by the peak file (takes values `correct`, `false positive` and `false negative`). `chunk` and `sample.id` contain corresponding info about the region.

Summary files:

Each row corresponds to one of the seven datasets. The columns `fpr`, `fnr` and `error.rate` are in percent; they represent an estimate of false positive rate, false negative rate and overall error percentage respectively for each data set.

6. If necessary, repeat steps 4.-5. with different FDR control and GAP values. (Alternatively, use the multiple parameter solutions given in 4. and 5.)

7. To compare SPAN to other peak callers investigated in the original paper, launch
```bash
Rscript benchmark_other.R
```
from the console.

The script will generate a `summary_other.csv` file in the working directory. It's a tab-separated unquoted table. Each row corresponds to a unique caller-dataset-parameter combination, described in the first three columns. The other columns (`fpr`, `fnr` and `error.rate`) have the same meaning as in 5.

8. To conduct a 4-fold cross-validation experiment similar to the one described in the article, launch
```bash
$ Rscript 4fold.R <FDRs> <GAPs>
```
where `FDRs` and `GAPs` are comma-separated lists of corresponding parameters' values, e.g.
```bash
$ Rscript 4fold.R 1E-4,1E-5,1E-6 3,4,5
```
The values should be exactly the same as the ones previously passed to `benchmark.R`.

The script will generate a set of `<dataset>_4fold_test.csv` files in the working directory. Each file will consist of four rows corresponding to the four folds of the cross-validation experiment, denoted by the `fold.id` column. For each fold, `optimal.fdr` and `optimal.gap` will show the trained values of the parameters, and the `training.error.rate` and `test.error.rate` are self-explanatory.

The files included in the repository were produced by calling
```bash
$ Rscript 4fold.R 1E-1,5E-2,1E-2,5E-3,1E-3,5E-4,1E-4,5E-5,1E-5,5E-6,1E-6 0,1,2,3,4,5,6,7,8,9,10
```
as well as all the corresponding previous steps.

9. To compare SPAN to other peak callers investigated in the original paper, launch
```bash
Rscript 4fold_other.R
```
from the console.

The script will generate a `summary_4fold_others.csv` file in the working directory. It's a tab-separated unquoted table. Each row corresponds to a unique caller-dataset-fold number combination, described in the first three columns. In then gives `training.error`, `optimal.param` and `test.error` values. Note that for the `.default` callers no training was done, thus `optimal.param` is always the same.

We've included a copy of the `summary_4fold_others.csv` file in the repository. Feel free to remove and recalculate it.

## Conclusion

Use the `<dataset>_*.csv` and `summary_*.csv` files to compare SPAN with other peak callers. In particular, compare the mean test error in the 4-fold cross-validation experiment.
