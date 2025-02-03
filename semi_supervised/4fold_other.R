# This R script calculates the test error rates of the peak callers tested in the original experiment.
# Usage: Rscript 4fold_other.R

# Load the error data frame.

pre <- "http://members.cbio.mines-paristech.fr/~thocking/chip-seq-chunk-db/"

errors.table <- read.table(paste0(pre, "RData.model.error.txt"), header = TRUE)

errors.table$dataset <- factor(sub("/.*", "", errors.table$file))
errors.table$caller <- factor(sub(".*/", "", sub("\\.RData", "", errors.table$file)))
errors.table$chunk <- factor(sub("/.*", "", sub("[^/]*/", "", errors.table$file)))
errors.table$target <- sub("_.*", "", errors.table$dataset)

test.folds <- read.table("./4foldcv-test-folds.csv", sep = ",", header = TRUE)
test.folds$dataset <- factor(sub("/.*", "", test.folds$test.chunk))
test.folds$chunk <- factor(sub(".*/", "", test.folds$test.chunk))

# Summarize errors for each dataset, each caller, each test fold.

summary.errors <- data.frame()

for (caller in unique(errors.table$caller)) {
	message(sprintf("Porcessing caller %s", caller))
	for (dataset in unique(errors.table$dataset)) {
		cat(sprintf("\tProcessing dataset %s", dataset))
		test.folds.dataset <- test.folds[test.folds$dataset == dataset,]		
		for (fold in unique(test.folds.dataset$fold.id)) {
			cat(".")
			test.chunks <- test.folds.dataset[test.folds.dataset$fold.id == fold,]$chunk
			training.chunks <- test.folds.dataset[test.folds.dataset$fold.id != fold,]$chunk
			dataset.training.errors <- data.frame()
			dataset.test.errors <- data.frame()
			for (i in which(errors.table$dataset == dataset & errors.table$caller == caller)) {
				errors.url <- url(paste0(pre, errors.table$file[[i]]))
				load(errors.url)
				close(errors.url)
				if (errors.table$chunk[[i]] %in% test.chunks) {
					dataset.test.errors <- rbind(dataset.test.errors, error)
				} else {
					dataset.training.errors <- rbind(dataset.training.errors, error)
				}
			}
			dataset.training.errors.param <- split(dataset.training.errors, dataset.training.errors$param.name)
			error.rates <- as.numeric(lapply(dataset.training.errors.param, function(df) sum(df$fp + df$fn) * 100 /nrow(df)))
			# Select optimal parameter configuration for the training set.
			optimal.training.error.rate <- min(error.rates)
			optimal.param <- dataset.training.errors.param[[which.min(error.rates)]]$param.name[[1]]
			test.error.rate <- with(dataset.test.errors[dataset.test.errors$param.name == optimal.param,], sum(fp + fn) * 100 / length(fp))
			summary.errors <- rbind(summary.errors, data.frame(dataset = dataset, caller = caller, fold.id = fold,
					training.error = optimal.training.error.rate, optimal.param = optimal.param, test.error = test.error.rate))
		}
		cat("\n")
	}
}

write.table(summary.errors, "./summary_4fold_others.csv", sep = "\t", quote = FALSE)
