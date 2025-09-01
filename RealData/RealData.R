library(knitr)
library(dplyr)

source("Evaluate_cos.R")

####################################################################################
# Table2: The AUC of detecting known relation pairs with different methods.
load("cosine_similarities.Rdata") # read in cosine similarities from different methods
AUClist <- lapply(cosine, function(x){
  Evaluate(cosU = x, LPrm = TRUE)
})
names(AUClist) <- names(cosine)
AUC <- cbind(do.call("cbind", lapply(AUClist, function(x) x[,2])), AUClist[[1]][,1])
# print(AUC)
kable(AUC, "latex", 3)

####################################################################################
# Table3: Power of detecting known relation pairs
load("cosine_similarities.Rdata") # read in cosine similarities from different methods
load("log_p_matrix.Rdata") # read in the log of p values for identifying whether each element is zero
load("AllRelationPairsWillNull.Rdata") # read in known relationship pairs with negative samples

powerlist <- lapply(cosine, function(x){
  return(Get_Power(log_p_matrix = log_p_matrix, cosine_similarity = x, 
                   AllRelationPairs = AllRelationPairs, target_FDR = c(0.05)))
})
names(powerlist) <- names(cosine)
power <- cbind(powerlist[[1]][,1], do.call("cbind", lapply(powerlist, function(x) x[,2])))
# print(power)
kable(power, "latex", 3)

####################################################################################
# Table4: Spearman rank correlation test between the score and GPT ratings
load("cosine_similarities.Rdata") # read in cosine similarities from different methods
load("log_p_matrix.Rdata") # read in the log of p values for identifying whether each element is zero
load("GPT_scores.Rdata") # read in the GPT scores for 100 clinical features with AD
N <- 1000
ref <- c('chatgpt','gpt4')
set.seed(214)
corrtest <- lapply(c(list('KNIT' = -log_p_matrix), cosine), function(x){
  mycosine <- x[match(score$code, rownames(x)),which(colnames(x)=="PheCode:290.11")]
  subcorr <- c()
  for(currref in ref){
    r <- cor(mycosine, score[[currref]], method = "spearman")
    r_permute <- lapply(1:N, function(o){
      s <- score[[currref]][sample(1:length(mycosine), length(mycosine))]
      return(cor(mycosine, s, method = "spearman"))
    })
    subcorr = c(cubcorr, r, 1-sum(r>r_permute)/(1+N))
  }
  return(subcorr)
})
corrtest <- do.call("cbind", corrtest)
colnames(corrtest) <- c('KNIT', names(cosine))
# print(corrtest)
kable(corrtest, "latex", 3)

####################################################################################
# Figure 6: The estimated low-rank PMI with the smallest p-values when quantifying their relationships with AD
load("cosine_similarities.Rdata") # read in cosine similarities from different methods
load("log_p_matrix.Rdata") # read in the log of p values for identifying whether each element is zero
cosine <- cosine[['SVDPPMI']]
idphe <- which(rownames(cosine)=="PheCode:290.11")
p.log <- log_p_matrix[idphe,-idphe]
names(p.log) <- rownames(cosine)
p.log <- sort(p.log)
df <- data.frame("myrank" = 1:length(p.log), 
                 "cosine" = cosine[idphe, match(names(p.log), colnames(cosine))],
                 "p_log" = p.log, 
                 "newp" = log(-p.log))
minrank <- max(which(df$p_log < (log(1:nrow(df))-log(length(p.log))-log(log(length(p.log))+1))+log(0.05)))
df$myrank <- 1:nrow(df)

library(ggplot2) 
ggp <- ggplot(df[1:5000,])  +  
  geom_bar(aes(x=myrank, y=cosine),stat="identity", fill="#00BFC4",colour="#00BFC4") + 
  geom_line(aes(x=myrank, y=newp/5),stat="identity",color="#F8766D", linewidth = 1.5, alpha = 0.8) + 
  labs(title= "", 
       x="Rank",y="low-rank PMI estimator")+ 
  scale_y_continuous(sec.axis=sec_axis(~.*0.2,name="log(-log(p))")) +
  geom_vline(xintercept=minrank, linewidth = 1.5, color = "black") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
ggp 

####################################################################################
# Figure 7: t-SNE
library(Rtsne)
library(ggplot2)
library(dplyr)
load("embedding.Rdata") # read in the embedding from SVD-PPMI
tsne_out <- Rtsne(embed, 3)
tsne_plot <- data.frame(x = tsne_out$Y[,1], 
                        y = tsne_out$Y[,2])
group <- sapply(rownames(embed), function(x) strsplit(x,":")[[1]][1])
group[which(group=="Other lab")] <- "Local Lab"
group[which(group=="ShortName")] <- "Local Lab"
tsne_plot$group <- group
tsne_plot$code <- rownames(embed)
tsne_plot$desc <- dict$desc[match(tsne_plot$code,dict$code)]
tsne_plot$group <- factor(tsne_plot$group, levels = c("PheCode", "Local Lab", "LOINC", "RXNORM", "CCS"))
tsne_plot <- tsne_plot %>% arrange(group)

p1 <- ggplot2::ggplot(tsne_plot,label=group) + 
  geom_point(aes(x=x,y=y,color=group),size=0.3) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(color = NULL) +
  theme_void() +
  theme(legend.title=element_text(size=20), 
        legend.text=element_text(size=20))

tsne_plot <- data.frame(x = tsne_out$Y[,2], 
                        y = tsne_out$Y[,3])
group <- sapply(rownames(embed), function(x) strsplit(x,":")[[1]][1])
group[which(group=="Other lab")] <- "Local Lab"
group[which(group=="ShortName")] <- "Local Lab"
tsne_plot$group <- group
tsne_plot$code <- rownames(embed)
tsne_plot$desc <- dict$desc[match(tsne_plot$code,dict$code)]
tsne_plot$group <- factor(tsne_plot$group, levels = c("PheCode", "Local Lab", "LOINC", "RXNORM", "CCS"))
tsne_plot <- tsne_plot %>% arrange(group)

p2 <- ggplot2::ggplot(tsne_plot,label=group) + 
  geom_point(aes(x=x,y=y,color=group),size=0.3) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(color = NULL) +
  theme_void() +
  theme(legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
