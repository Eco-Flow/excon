#!/usr/bin/Rscript

library(pheatmap)

# Helper function to safely plot pheatmap
safe_pheatmap <- function(mat, outfile, width=9, height=5, my_palette) {
    pdf(outfile, width=width, height=height)
    if (nrow(mat) == 0 || ncol(mat) == 0) {
        plot.new()
        text(0.5, 0.5, "No data available", cex=1.5)
    } else if (length(unique(as.vector(mat))) <= 1) {
        plot.new()
        text(0.5, 0.5, "No significant GO terms found", cex=1.5)
    } else {
        tryCatch({
            pheatmap::pheatmap(log10(mat), col=my_palette, cluster_rows=F,
                               treeheight_row=0, treeheight_col=0, legend=T)
        }, error = function(e) {
            plot.new()
            text(0.5, 0.5, paste("Plot failed:", e$message), cex=1.2)
        })
    }
    dev.off()
}

my_palette <- colorRampPalette(c("purple","red", "orange", "yellow", "white"))(n=20)

# --- Negative GO summary ---
erefd <- read.table("Go_summary_neg.tsv", h=T, sep="\t")
newdata <- erefd[order(erefd$Count_significant, decreasing=T),]
rownames(newdata) <- paste(newdata$GO_ID, newdata$GO_term)
df <- subset(newdata, select=-c(GO_term, GO_ID, Count_significant))
df[is.na(df)] <- 1
df2 <- as.matrix(df)

safe_pheatmap(head(df2, n=30), "Go_summary_neg.pdf", my_palette=my_palette)

mat_clean <- df2[, !grepl("Node", colnames(df2)), drop=FALSE]
safe_pheatmap(head(mat_clean, n=30), "Go_summary_neg_noNode.pdf", my_palette=my_palette)

# --- Positive GO summary ---
erefd <- read.table("Go_summary_pos.tsv", h=T, sep="\t")
newdata <- erefd[order(erefd$Count_significant, decreasing=T),]
rownames(newdata) <- paste(newdata$GO_ID, newdata$GO_term)
df <- subset(newdata, select=-c(GO_term, GO_ID, Count_significant))
df[is.na(df)] <- 1
df2 <- as.matrix(df)

safe_pheatmap(head(df2, n=30), "Go_summary_pos.pdf", my_palette=my_palette)

mat_clean2 <- df2[, !grepl("Node", colnames(df2)), drop=FALSE]
safe_pheatmap(head(mat_clean2, n=30), "Go_summary_pos_noNode.pdf", my_palette=my_palette)
