# Statistical analyses reported in "A new exceptionally preserved biota from
# South China expands Silurian biodiversity".
#
# Run from the repository root:  Rscript analysis.R
# Inputs are Data S2 (census) and Data S3 (major groups), as published.
# Writes cluster membership and Z-scores (Fig. 4A) to output/.

census <- readxl::read_excel("data/dataS2_census.xlsx", sheet = "Data",
                             na = c("", "NA"))
groups <- readxl::read_excel("data/dataS3_major_groups.xlsx", sheet = "Data",
                             range = "A1:AS29")
dir.create("output", showWarnings = FALSE)

# ---- Diversity -------------------------------------------------------------
counts <- census$Count
N <- sum(counts)
cat("Taxa:", length(counts), "  Specimens:", N, "\n")

# Zahl's unbiased Simpson index
zahl <- 1 - sum(counts * (counts - 1)) / (N * (N - 1))
cat("Zahl's unbiased Simpson index:", round(zahl, 2), "\n")

# ACE (Chao & Lee 1992), with Chao's (1987) log-normal 95% confidence interval,
# which is bounded below by the observed richness
rich <- vegan::estimateR(counts)
unseen <- rich[["S.ACE"]] - rich[["S.obs"]]
C <- exp(1.96 * sqrt(log(1 + rich[["se.ACE"]]^2 / unseen^2)))
cat(sprintf("ACE: %.1f (95%% CI %.0f-%.0f)\n", rich[["S.ACE"]],
            rich[["S.obs"]] + unseen / C, rich[["S.obs"]] + unseen * C))

# ---- Trophic composition (Fig. 4, C and D) ---------------------------------
Proportions <- function(category) {
  tab <- tapply(census$Count, category, sum)
  round(100 * tab / sum(tab))
}
cat("\nFeeding habit (% of specimens):\n")
print(Proportions(census$`Feeding habit`))
cat("\nSubstrate niche (% of specimens):\n")
print(Proportions(census$`Substrate niche`))

# ---- Faunal similarity (Fig. 4A, fig. S8) ----------------------------------
presence <- as.data.frame(groups[, -(1:2)])
rownames(presence) <- groups$Short_title

set.seed(0)
rcDist <- vegan::raupcrick(presence, nsimul = 2^14)

# Linkage with the highest cophenetic correlation; median and centroid
# linkage can produce non-monotonic dendrograms, so are not considered
linkages <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty")
trees <- lapply(linkages, function(method) hclust(rcDist, method = method))
cophCor <- vapply(trees, function(tr) cor(rcDist, cophenetic(tr)), numeric(1))
names(cophCor) <- linkages
cat("\nCophenetic correlation:\n")
print(round(cophCor, 3))
hc <- trees[[which.max(cophCor)]]

# Number of clusters maximizing mean silhouette width; k-means++ for comparison
kRange <- 2:20
silHC <- silKPP <- setNames(numeric(length(kRange)), kRange)
for (k in kRange) {
  silHC[[as.character(k)]] <- summary(
    cluster::silhouette(cutree(hc, k), rcDist))$avg.width
  silKPP[[as.character(k)]] <- summary(
    cluster::silhouette(TreeDist::KMeansPP(rcDist, k)$cluster, rcDist))$avg.width
}
kBest <- kRange[which.max(silHC)]
cat(sprintf("\n%s linkage: maximum silhouette %.2f with %d clusters\n",
            linkages[which.max(cophCor)], max(silHC), kBest))
cat(sprintf("k-means++: maximum silhouette %.3f\n", max(silKPP)))

cluster <- cutree(hc, kBest)
write.csv(data.frame(Biota = groups$Biota, Cluster = cluster),
          "output/clusters.csv", row.names = FALSE)
cat("\nCluster membership:\n")
print(split(groups$Short_title, cluster))

# Z-score of each group's occurrence frequency across clusters (Fig. 4A)
freq <- t(aggregate(presence, by = list(cluster), FUN = mean)[, -1])
colnames(freq) <- paste("Cluster", seq_len(kBest))
z <- t(scale(t(freq)))
write.csv(round(z, 3), "output/zscores.csv")

sessionInfo()
