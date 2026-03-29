#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
opt <- list(
    base     = args[1],
    gamma    = args[2],
    gamma_pf = args[3]
)

# --- Helpers ---

parse_nll <- function(results_dir, label) {
    files <- list.files(results_dir, pattern="^[^_]+_results\\.txt$", full.names=TRUE, recursive=TRUE)
    if (length(files) == 0) {
        warning(paste("No results file found in", results_dir, "- returning NA"))
        return(NA)
    }
    lines <- readLines(files[1])
    score_line <- grep("Final Likelihood", lines, value=TRUE)
    if (length(score_line) == 0) {
        warning(paste("No 'Final Likelihood' line found in", files[1]))
        return(NA)
    }
    val <- as.numeric(gsub(".*Final Likelihood \\(-lnL\\):\\s*", "", score_line[1]))
    cat(label, "NLL:", val, "\n")
    val
}

aic <- function(nll, k) 2 * k + 2 * nll   # nll is already positive

lrt  <- function(nll_null, nll_alt, df) {
    stat <- 2 * (nll_null - nll_alt)       # both positive, null >= alt
    list(
        statistic = stat,
        p_value   = pchisq(stat, df=df, lower.tail=FALSE)
    )
}

# --- Parse scores ---

nll <- list(
    base            = parse_nll(opt$base,     "Base"),
    gamma           = parse_nll(opt$gamma,    "Gamma"),
    gamma_per_family= parse_nll(opt$gamma_pf, "Gamma_pf")
)

# Number of free parameters per model
k <- list(
    base             = 1,   # lambda only
    gamma            = 2,   # lambda + alpha shape
    gamma_per_family = 3    # lambda + alpha + Poisson root
)

# --- AIC table ---

results <- data.frame(
    model     = names(nll),
    NLL       = unlist(nll),
    params    = unlist(k),
    AIC       = mapply(aic, unlist(nll), unlist(k)),
    row.names = NULL
)
results$delta_AIC <- results$AIC - min(results$AIC)
results$weight    <- {
    w <- exp(-0.5 * results$delta_AIC)
    round(w / sum(w), 4)
}
results$best <- results$delta_AIC == 0

# --- LRT: base vs gamma (nested, df=1) ---

lrt_base_vs_gamma <- lrt(nll$base, nll$gamma, df=1)

results$LRT_vs_base_stat <- NA
results$LRT_vs_base_p    <- NA
results$LRT_vs_base_stat[results$model == "gamma"] <- lrt_base_vs_gamma$statistic
results$LRT_vs_base_p[results$model    == "gamma"] <- lrt_base_vs_gamma$p_value

# --- Write outputs ---

write.table(results, "cafe_model_comparison.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

best_model <- results$model[results$best]
writeLines(best_model, "best_model.txt")

cat("\n=== CAFE Model Comparison ===\n")
print(results[, c("model","NLL","params","AIC","delta_AIC","weight","best")])
cat("\nLRT base vs gamma: statistic =", round(lrt_base_vs_gamma$statistic, 3),
    ", p =", signif(lrt_base_vs_gamma$p_value, 3), "\n")
cat("Best model:", best_model, "\n")
