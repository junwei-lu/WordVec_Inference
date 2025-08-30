library(readr)
library(dplyr)
library(Matrix)
p <- 100
load("SVD_SPPMI.Rdata")
cosine_similarity <- ans$fit$u[,1:p] %*% diag(sqrt(ans$fit$d[1:p]))
rownames(cosine_similarity) <- rownames(ans$fit$u)
rm(ans)
cosine_similarity <- cosine_similarity/apply(cosine_similarity, 1, norm, '2')
cosine_similarity <- cosine_similarity %*% t(cosine_similarity)
coscov <- read_csv("cosineTest.csv")
coscov <- as.matrix(coscov)
pvalue <- cosine_similarity/sqrt(coscov)
rm(coscov)
pvalue <- pnorm(-pvalue, log.p = TRUE)
pvalue[lower.tri(pvalue, diag = TRUE)] <- 0
pvalue[pvalue>log(0.05)] <- 0
pvalue <- as(pvalue, "sparseMatrix")
namelist <- rownames(pvalue)
myn <- length(namelist)
myn <- myn*(myn-1)/2
pvalue <- summary(pvalue) %>%
  arrange(x) 
pvalue$thre <- log((1:nrow(pvalue))/myn*0.05)
idx <- max(which(pvalue$x <= pvalue$thre))
pvalue <- pvalue[1:idx,] %>%
  select(i, j)
pvalue$cos <- cosine_similarity[cbind(pvalue$i, pvalue$j)]
pvalue <- pvalue %>%
  mutate(i = namelist[i], j = namelist[j])
print(nrow(pvalue))
print(nrow(pvalue)/myn)
write.csv(pvalue, file = "KG_pvalue.csv", row.names = FALSE)