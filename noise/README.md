This pipeline is based on https://github.com/JetBrains-Research/chipseq-smk-pipeline

1. BAM samples are in <work_dir>/<histone_modifications>

2. Uncomment the line `rules.all_snr_results.input` in `Snakefile` to compute all signal to noise ratios.

3. Pick best quality tracks for each modification.

4. Run `bash run.bash` to launch peak callers with different settings.

5. Analyze results in the notebook `SPAN noise experiment`.

