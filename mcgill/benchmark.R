# This R script processes the output of "call_peaks.sh" and generates several error tables.
# Usage: Rscript benchmark.R [<FDRs> [<GAPs>]]
# <FDRs> and <GAPs> are comma-separated lists of FDR control values and GAP parameters, respectively.
# The default values are 1E-4 and 5.
# It's expected that "call_peaks.sh" had already been called for every combination of these parameters.

# Install packages if necessary.

if (!require(PeakError)) {
	if(!require(devtools)) { install.packages("devtools", repos = "http://cran.us.r-project.org") }
	devtools::install_github("tdhock/PeakError")
}

library(PeakError)

# Parse command line arguments

args <- commandArgs(TRUE)

if (length(args) > 0) {
	fdrs <- strsplit(args[1], ",")[[1]]	
} else {
	fdrs <- "1E-4"
}

if (length(args) > 1) {
	gaps <- strsplit(args[2], ",")[[1]]	
} else {
	gaps <- "5"
}

# Load the table of annotated regions.

pre <- "http://members.cbio.mines-paristech.fr/~thocking/chip-seq-chunk-db/"

regions.table <- read.table(paste0(pre, "RData.annotated.regions.txt"), header = TRUE)

# Extract dataset, chunk and ChIP-seq target information from the table.

regions.table$dataset <- factor(sub("/.*", "", regions.table$file))

regions.table$chunk <- as.numeric(sub(".*/", "", sub("/regions.RData", "", regions.table$file)))

regions.table$target <- sub("_.*", "", regions.table$dataset)

# Compute the error for each dataset, each chunk of the dataset, each sample in the chunk, and each FDR and GAP provided to the script.

summary.errors <- data.frame(dataset = c(), fpr = c(), fnr = c(), error.rate = c(), fdr = c(), gap = c())

for (dataset in levels(regions.table$dataset)) {
	message(sprintf("Processing dataset %s...", dataset))
	dataset.errors <- data.frame()
	for (i in which(regions.table$dataset == dataset)) {
		chunk <- regions.table$chunk[[i]]
		cat(sprintf("\tchunk %d", chunk))
		region.url <- url(paste0(pre, regions.table$file[[i]]))
		load(region.url)
		close(region.url)
		sample.ids <- levels(factor(regions$sample.id))
		for (sample.id in sample.ids) {
			sample.no <- sub("McGill", "", sample.id)
			for (fdr in fdrs) {
				for (gap in gaps) {
					peak.file <- paste0("./chip-seq-benchmark/", regions.table$target[[i]], "/peaks", sample.no, "_", gap, "_", fdr, ".bed")
					if (file.size(peak.file) > 0) {
						peaks <- read.table(peak.file, header = FALSE)
						colnames(peaks)[1:3] <- c("chrom", "chromStart", "chromEnd")
					} else {
						# For simplicity, we create a non-empty mock dataframe
						peaks <- data.frame(chrom = "chr1", chromStart = 1, chromEnd = 2)
					}
					dataset.errors <- rbind(dataset.errors, cbind(PeakError(peaks, regions[regions$sample.id == sample.id,]), chunk = chunk, sample.id = sample.id, fdr = fdr, gap = gap))
				}
			}
			cat(".")
		}
		cat("\n")

	}
	# Output a detailed dataset-FDR-GAP error table and add a line to the summary table.
	for (fdr in fdrs) {
		for (gap in gaps) {
			current.dataset.errors <- dataset.errors[dataset.errors$fdr == fdr & dataset.errors$gap == gap,]
			write.table(current.dataset.errors, sprintf("./%s_%s_%s.csv", dataset, gap, fdr), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
			fp <- sum(current.dataset.errors$fp)
			fn <- sum(current.dataset.errors$fn)
			possible.fp <- sum(current.dataset.errors$possible.fp)
			possible.fn <- sum(current.dataset.errors$possible.tp)
			fpr <- fp * 100 / possible.fp
			fnr <- fn * 100 / possible.fn
			error.rate <- (fp + fn) * 100 / nrow(current.dataset.errors)
			summary.errors <- rbind(summary.errors, data.frame(dataset = dataset, fpr = fpr, fnr = fnr, error.rate = error.rate, fdr = fdr, gap = gap))
		}
	}
}

# Output the summary table for each FDR and GAP provided to the script.

for (fdr in fdrs) {
	for (gap in gaps) {
		write.table(summary.errors[summary.errors$fdr == fdr & summary.errors$gap == gap,], paste0("./summary_", gap, "_", fdr, ".csv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
	}
}

