# This R script summarizes the error rates of the peak callers tested in the original experiment.
# Usage: Rscript benchmark_other.R

# Load the error data frame.

pre <- "http://members.cbio.mines-paristech.fr/~thocking/chip-seq-chunk-db/"

errors.table <- read.table(paste0(pre, "RData.model.error.txt"), header = TRUE)

# Extract dataset, peak caller, chunk and ChIP-seq target information from the table.

errors.table$dataset <- factor(sub("/.*", "", errors.table$file))

errors.table$caller <- factor(sub(".*/", "", sub("\\.RData", "", errors.table$file)))

errors.table$chunk <- as.numeric(sub("/.*", "", sub("[^/]*/", "", errors.table$file)))

errors.table$target <- sub("_.*", "", errors.table$dataset)

# Summarize errors for each dataset and each caller.

summary.errors <- data.frame(caller = c(), dataset = c(), fpr = c(), fnr = c(), error.rate = c())

for (dataset in levels(errors.table$dataset)) {
	message(sprintf("Processing dataset %s...", dataset))
	for (caller in levels(errors.table$caller)) {
		message(sprintf("\tProcessing caller %s...", caller))
		dataset.errors <- data.frame()
		for (i in which(errors.table$dataset == dataset & errors.table$caller == caller)) {
			message(sprintf("\t\tchunk %d", errors.table$chunk[[i]]))
			errors.url <- url(paste0(pre, errors.table$file[[i]]))
			load(errors.url)
			close(errors.url)
			dataset.errors <- rbind(dataset.errors, error)
		}
		params <- levels(dataset.errors$param.name)
		for (param in params) {
			current.errors <- dataset.errors[dataset.errors$param.name == param,]
			fp <- sum(current.errors$fp)
			fn <- sum(current.errors$fn)
			possible.fp <- sum(current.errors$possible.fp)
			possible.fn <- sum(current.errors$possible.tp)
			fpr <- fp * 100 / possible.fp
			fnr <- fn * 100 / possible.fn
			error.rate = (fp + fn) * 100 / nrow(current.errors)
			summary.errors <- rbind(summary.errors, data.frame(caller = caller, dataset = dataset, param.name = param, fpr = fpr, fnr = fnr, error.rate = error.rate))
		}
	}
}

# Output the summary error rate table.

write.table(summary.errors, "./summary_other.csv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

