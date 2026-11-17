# Niulan biota: statistical analyses

Code to reproduce the statistical results reported in

> J. Yang *et al.*, A new exceptionally preserved biota from South China
> expands Silurian biodiversity.

## Data

| File | Manuscript | Content |
|------|------------|---------|
| `data/dataS2_census.xlsx` | Data S2 | Specimen counts, feeding habit and substrate niche of each metazoan taxon |
| `data/dataS3_major_groups.xlsx` | Data S3 | Presence of major groups in Silurian Lagerstätten |

## Running

```r
install.packages(c("readxl", "vegan", "cluster", "TreeDist"))
```

```sh
Rscript analysis.R
```

`analysis.R` prints each statistic and writes cluster membership
(`output/clusters.csv`) and the Z-scores shaded in Fig. 4A
(`output/zscores.csv`).

## Expected results

| Statistic | Value |
|-----------|-------|
| Taxa; specimens | 70; 1567 |
| Zahl's unbiased Simpson index | 0.96 |
| ACE estimate of species richness | 72.5 (95% CI 70–94) |
| Feeding habits: suspension feeders, predators, detritivores | 66%, 17%, 17% |
| Substrate niches: epifaunal suspension feeders, infaunal suspension feeders, vagrant detritivores | 55%, 23%, 22% |
| Linkage with highest cophenetic correlation | UPGMA (average) |
| Optimal clustering | 12 clusters; silhouette 0.56 |
| k-means++ | silhouette < 0.42 |

Raup–Crick dissimilarities are estimated by simulation (`nsimul = 2^14`,
`set.seed(0)`); other seeds give slightly different values.
Results were generated with R (development version, 2026-06-29), readxl 1.5.0, vegan 2.7.5, cluster 2.1.8 and
TreeDist 2.15.

## Licence

GPL-3; see `LICENSE`.
