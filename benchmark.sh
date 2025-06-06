# Keep several datasets for benchmarking under folders bam1, bam2, ...
cd ~/data/2025_benchmark

for BAMSF in bams1 bams2 bams3 bams4 bams5; do

	echo $BAMSF; rm -rf $(pwd)/bams; ln -sf $(pwd)/$BAMSF bams;

	# Default
	time snakemake -p -s ~/work/chipseq-smk-pipeline/Snakefile all --use-conda --cores all  --directory $(pwd) --config fastq_dir=$(pwd) start_with_bams=True genome=hg38 macs2=True --rerun-incomplete --rerun-triggers mtime;
	time snakemake -p -s ~/work/chipseq-smk-pipeline/Snakefile all --use-conda --cores all  --directory $(pwd) --config fastq_dir=$(pwd) start_with_bams=True genome=hg38 macs2=True macs2_mode=broad macs2_params="--broad --broad-cutoff 0.1" macs2_suffix=broad0.1 --rerun-incomplete --rerun-triggers mtime;
	time snakemake -p -s ~/work/chipseq-smk-pipeline/Snakefile all --use-conda --cores all  --directory $(pwd) --config fastq_dir=$(pwd) start_with_bams=True genome=hg38 sicer=True --rerun-incomplete --rerun-triggers mtime;
	time snakemake -p -s ~/work/chipseq-smk-pipeline/Snakefile all --use-conda --cores all  --directory $(pwd) --config fastq_dir=$(pwd) start_with_bams=True genome=hg38 span=True  --rerun-incomplete --rerun-triggers mtime;

	# MACS3 / SICER2
	time snakemake -p -s ~/work/chipseq-smk-pipeline/Snakefile all --use-conda --cores all  --directory $(pwd) --config fastq_dir=$(pwd) start_with_bams=True genome=hg38 macs3=True --rerun-incomplete --rerun-triggers mtime;
	time snakemake -p -s ~/work/chipseq-smk-pipeline/Snakefile all --use-conda --cores all  --directory $(pwd) --config fastq_dir=$(pwd) start_with_bams=True genome=hg38 macs3=True macs3_mode=broad macs3_params="--broad --broad-cutoff 0.1" macs3_suffix=broad0.1 --rerun-incomplete --rerun-triggers mtime;
	time snakemake -p -s ~/work/chipseq-smk-pipeline/Snakefile all --use-conda --cores all  --directory $(pwd) --config fastq_dir=$(pwd) start_with_bams=True genome=hg38 sicer2=True --rerun-incomplete --rerun-triggers mtime;

	# Ok
	time snakemake -p -s ~/work/chipseq-smk-pipeline/Snakefile all --use-conda --cores all  --directory $(pwd) --config fastq_dir=$(pwd) start_with_bams=True genome=hg38 homer=True --rerun-incomplete --rerun-triggers mtime;
	time snakemake -p -s ~/work/chipseq-smk-pipeline/Snakefile all --use-conda --cores all  --directory $(pwd) --config fastq_dir=$(pwd) start_with_bams=True genome=hg38 hotspot=True hotspot_executable=~/data/2025_peak_callers/hotspot-4.1.1/hotspot-distr/hotspot-deploy/bin/hotspot --rerun-incomplete --rerun-triggers mtime;
	time snakemake -p -s ~/work/chipseq-smk-pipeline/Snakefile all --use-conda --cores all  --directory $(pwd) --config fastq_dir=$(pwd) start_with_bams=True genome=hg38 peakseq=True peakseq_executable=~/data/2025_peak_callers/PeakSeq/bin/PeakSeq --rerun-incomplete --rerun-triggers mtime;

	# Slow
	time snakemake -p -s ~/work/chipseq-smk-pipeline/Snakefile all --use-conda --cores all  --directory $(pwd) --config fastq_dir=$(pwd) start_with_bams=True genome=hg38 fseq2=True --rerun-incomplete --rerun-triggers mtime;

	# Extremely slow, more than 1 hour per sample
	# time snakemake -p -s ~/work/chipseq-smk-pipeline/Snakefile all --use-conda --cores all  --directory $(pwd) --config fastq_dir=$(pwd) start_with_bams=True genome=hg38 gps=True --rerun-incomplete --rerun-triggers mtime;
	# time snakemake -p -s ~/work/chipseq-smk-pipeline/Snakefile all --use-conda --cores all  --directory $(pwd) --config fastq_dir=$(pwd) start_with_bams=True genome=hg38 bayespeak=True bayespeak_rscript_executable=~/miniconda3/envs/bayespeak/bin/Rscript --rerun-incomplete --rerun-triggers mtime;

	stop;
done;
