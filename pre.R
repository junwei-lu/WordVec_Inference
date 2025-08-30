library(readr)
library(dplyr)
library(rsvd)
library(Matrix)
get_cooc = function(cooc, dict){
  cooc$i <- dict$feature_id[match(cooc$i, dict$WordIndex)]
  cooc$j <- dict$feature_id[match(cooc$j, dict$WordIndex)]
  return(cooc)
}
get_SPPMI = function(colist, thre_co = 100, thre_rs = 300){
  mapdict <- data.frame(code = unique(c(unique(colist$i), unique(colist$j))))
  mapdict$WordIndex <- 1:nrow(mapdict)
  colist$i <- match(colist$i, mapdict$code)
  colist$j <- match(colist$j, mapdict$code)
  colist <- sparseMatrix(i = colist$i, j = colist$j, x = colist$count, 
                         dims = c(nrow(mapdict), nrow(mapdict)),
                         dimnames = list(mapdict$code, mapdict$code))
  colist <- colist + t(colist)
  diag(colist) <- diag(colist)/2
  colist[which(colist<thre_co)] <- 0
  rs <- rowSums(colist)
  names(rs) <- mapdict$code
  rm(mapdict)
  idx <- which(rs >= thre_rs)
  colist <- colist[idx,idx]
  rs <- rowSums(colist)
  s <- sum(rs)
  PMI <- t((colist*s)/rs)/rs
  rm(colist, idx)
  PMI <- log(PMI)
  PMI[PMI<0] <- 0
  idx <- which(rowSums(PMI!=0)>0)
  PMI <- PMI[idx,idx]
  PMI <- as(PMI, "sparseMatrix")
  set.seed(214)
  fit <- rsvd(PMI, k = min(5000,nrow(PMI)))
  idx <- which(sign(fit$u[1,])==sign(fit$v[1,]))
  fit <- list(u = fit$u[,idx], d = fit$d[idx])
  rownames(fit$u) <- rownames(PMI)
  return(list(fit = fit, rs = rs))
}
dict <- read_csv("worddict.csv")
cooc <- read_csv("cooc.csv")
cooc <- get_cooc(cooc, dict)
ans <- get_SPPMI(cooc, thre_co = 200, thre_rs = 500)
save(ans, file = "SVD_SPPMI.Rdata")


## compute the pre_info for calculating the variance of cosine similarity
library(rhdf5)
library(h5)
p <- 100
load("SVD_SPPMI.Rdata")
## input the following information based on pre.py
n # number of patient
len # days/patient
dw # number of words/day/patient
## calculate the following information
q = 30*dw # window size = 30 days
Tw = dw*len # T: length of word sequence 
T1 = Tw*q-q*(q+1)/2 
p_w = ans$rs/(n*(2*Tw*q-q^2-q))
write.csv(p_w, file = "p_w.csv", row.names = FALSE)
nlist <- c(nlist, n)
T1list <- c(T1list, T1)
U = ans$fit$u[,1:p]
cosine_similarity <- ans$fit$u[,1:p] %*% diag(sqrt(ans$fit$d[1:p]))
PMI_lr <- cosine_similarity %*% t(cosine_similarity)
cosine_similarity <- cosine_similarity/apply(cosine_similarity, 1, norm, '2')
cosine_similarity <- cosine_similarity %*% t(cosine_similarity)
myA = cosine_similarity/PMI_lr
write.csv(myA, file = "myA.csv", row.names = FALSE)
P = U%*%t(U)
write.csv(P, file = "P.csv", row.names = FALSE)
PDP = P%*%diag(1/p_w)%*%P
write.csv(PDP, file = "PDP.csv", row.names = FALSE)