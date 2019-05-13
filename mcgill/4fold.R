# Parse command-line arguments

args <- commandArgs(TRUE)
if (length(args) != 2) {
	stop("Usage:\n\tRscript 4fold.R <fdrs> <gaps>\nwhere <fdrs> and <gaps> are comma-separated lists of corresponding parameters.")
}
fdrs <- strsplit(args[1], ",")[[1]]	
gaps <- strsplit(args[2], ",")[[1]]

# Read the test fold description

test.folds <- read.table("./4foldcv-test-folds.csv", sep = ",", header = TRUE)
test.folds$dataset <- factor(sub("/.*", "", test.folds$test.chunk))
test.folds$chunk <- factor(sub(".*/", "", test.folds$test.chunk))

# Run cross-validation for each dataset

for (dataset in levels(test.folds$dataset)) {
	message(sprintf("Processing dataset %s...", dataset))
	
	# Create a data frame for all parameter combinations
	dataset.errors <- data.frame()
	roc.table <- data.frame()
	for (fdr in fdrs) {
		for (gap in gaps) {
			table.path <- sprintf("./%s_%s_%s.csv", dataset, gap, fdr)
			current.errors <- read.table(table.path, header = TRUE, sep = "\t", colClasses = c(gap = "factor", fdr = "factor"))
			dataset.errors <- rbind(dataset.errors, current.errors)
			fpr <- with(current.errors, sum(fp) / sum(possible.fp))
			tpr <- with(current.errors, sum(tp) / sum(possible.tp))
			roc.table <- rbind(roc.table, data.frame(fdr = fdr, gap = gap, fpr = fpr, tpr = tpr))
		}
	}

	# Plot ROC curves
	png(sprintf("./roc_%s.png", dataset))
	with(roc.table[roc.table$fdr == fdrs[[1]],], plot(tpr ~ fpr, type = "b", col = fdr, xlim = c(min(roc.table$fpr), max(roc.table$fpr)*1.2), ylim = c(min(roc.table$tpr), max(roc.table$tpr)), main = dataset))
	for (fdr in fdrs[2:length(fdrs)]) {
		with(roc.table[roc.table$fdr == fdr,], lines(tpr ~ fpr, type = "b", col = fdr))
	}
	legend("bottomright", legend = unique(roc.table$fdr), fill = unique(roc.table$fdr))
	dev.off()

	test.errors <- data.frame()
	for (fold in unique(test.folds[test.folds$dataset == dataset,]$fold.id)) {
		message(sprintf("\tfold id %s", fold))
		# Define training and test sets
		test.chunks <- test.folds[test.folds$dataset == dataset & test.folds$fold.id == fold,]$chunk
		training.chunks <- test.folds[test.folds$dataset == dataset & test.folds$fold.id != fold,]$chunk
		training.set <- dataset.errors[dataset.errors$chunk %in% training.chunks,]
		test.set <- dataset.errors[dataset.errors$chunk %in% test.chunks,]

		# Calculate aggregated errors for training set
		training.errors <- data.frame()
		for (fdr in fdrs) {
			for (gap in gaps) {
				current.training.set <- training.set[training.set$fdr == fdr & training.set$gap == gap,]
				fp <- sum(current.training.set$fp)
				fn <- sum(current.training.set$fn)
				error.rate <- (fp + fn) * 100 / nrow(current.training.set)
				training.errors <- rbind(training.errors, data.frame(fdr = fdr, gap = gap, error.rate = error.rate))
			}
		}

		# Determine the optimal parameter combination and calculate the test error
		optimal.index <- which.min(training.errors$error.rate)
		training.error.rate <- training.errors$error.rate[[optimal.index]]
		optimal.fdr <- training.errors$fdr[[optimal.index]]
		optimal.gap <- training.errors$gap[[optimal.index]]
		optimal.test.set <- test.set[test.set$fdr == optimal.fdr & test.set$gap == optimal.gap,]
		fp <- sum(optimal.test.set$fp)
		fn <- sum(optimal.test.set$fn)
		test.error.rate <- (fp + fn) * 100 / nrow(optimal.test.set)
		test.errors <- rbind(test.errors, data.frame(fold.id = fold, optimal.fdr = optimal.fdr, optimal.gap = optimal.gap, training.error.rate = training.error.rate, test.error.rate = test.error.rate))		
	}
	write.table(test.errors, sprintf("./%s_4fold_test.csv", dataset), sep = "\t", quote = FALSE, row.names = FALSE)
	message(sprintf("\tMean test error: %g", mean(test.errors$test.error.rate)))
}
