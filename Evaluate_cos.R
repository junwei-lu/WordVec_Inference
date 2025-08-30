library(dplyr)
library(pROC)

load("pairs.Rdata") # read in the known relationship pairs

## functions for computing AUC for detecting known relationship pairs ####
Cos_boostrap<-function(Cos,dict_keep,dicts,g1,g2){
  dicts$id1 = match(dicts$feature_id1,dict_keep$feature_id)
  dicts$id2 = match(dicts$feature_id2,dict_keep$feature_id)
  dicts = dicts %>% na.omit()
  Cos_h = Cos[cbind(dicts$id1,dicts$id2)]
  
  num = nrow(dicts)
  #print(paste('Number of Pairs:',num))
  
  temp_dict = data.frame(id1=c(dicts$id1,dicts$id2),
                         id2=c(dicts$id2,dicts$id1))
  
  set.seed(1)
  IdA = sample(which(dict_keep$group%in%g1),size=2*num,replace=T)
  IdB = sample(which(dict_keep$group%in%g2),size=2*num,replace=T)
  
  temp_ID = data.frame(id1=IdA,id2=IdB)
  temp_ID = setdiff(temp_ID,temp_dict)
  temp_ID = subset(temp_ID,id1!=id2)
  print(paste('Number of Non-overlap:',nrow(temp_ID)))
  temp_ID = temp_ID[1:num,]
  
  Cos_null = Cos[cbind(temp_ID$id1,temp_ID$id2)]
  
  return(list('num'=num,'Cos_h'=Cos_h,'Cos_null'=Cos_null,p=c(Cos_h,Cos_null),y=c(rep(1,length(Cos_h)),rep(0,length(Cos_null)) )))
}

ROC<-function(Cos,dict_keep,dicts,g1,g2){
  Boot = Cos_boostrap(Cos,dict_keep,dicts,g1,g2)
  roc0 = roc(Boot$y, Boot$p)
  cut0 = as.numeric(quantile(Boot$Cos_null,0.99))
  cut1 = as.numeric(quantile(Boot$Cos_null,0.95))
  cut2 = as.numeric(quantile(Boot$Cos_null,0.9))
  
  tpr0 = roc0$sensitivities[which(roc0$thresholds>cut0)[1]]
  tpr1 = roc0$sensitivities[which(roc0$thresholds>cut1)[1]]
  tpr2 = roc0$sensitivities[which(roc0$thresholds>cut2)[1]]
  c('#pairs'=Boot$num,'auc'=roc0$auc,'cut/0.01'=cut0[1],'cut/0.05'=cut1[1],'cut/0.1'=cut2[1],
    'TPR/0.01'=tpr0,'TPR/0.05'=tpr1,'TPR/0.1'=tpr2)
}

Lab_PheCode<-function(Cos,dict_keep,dicts,g1,g2){
  dicts = subset(dicts,(feature_id1%in%dict_keep$feature_id)&(feature_id2%in%dict_keep$feature_id))
  dicts$id1 = match(dicts$feature_id1, rownames(Cos))
  dicts$id2 = match(dicts$feature_id2, rownames(Cos))
  dicts$cos = Cos[cbind(dicts$id1,dicts$id2)]
  Boot = Cos_boostrap(Cos,dict_keep,dicts,g1,g2)
  dicts$cos = Boot$Cos_h 
  temp = dicts %>% group_by(combination) %>% summarise(cos=mean(cos))
  
  num = length(temp$cos)
  y = c(rep(1,num),rep(0,num))
  p = c(temp$cos,Boot$Cos_null[1:num])
  
  roc0 = roc(y,p)
  cut0 = as.numeric(quantile(Boot$Cos_null[1:num],0.99))
  cut1 = as.numeric(quantile(Boot$Cos_null[1:num],0.95))
  cut2 = as.numeric(quantile(Boot$Cos_null[1:num],0.9))
  
  tpr0 = roc0$sensitivities[which(roc0$thresholds>cut0)[1]]
  tpr1 = roc0$sensitivities[which(roc0$thresholds>cut1)[1]]
  tpr2 = roc0$sensitivities[which(roc0$thresholds>cut2)[1]]
  c('#pairs'=num,'auc'=roc0$auc,'cut/0.01'=cut0[1],'cut/0.05'=cut1[1],'cut/0.1'=cut2[1],
    'TPR/0.01'=tpr0,'TPR/0.05'=tpr1,'TPR/0.1'=tpr2)
}

SNR<-function(Cos,dict_keep,dicts,g1='PheCode',g2='PheCode'){
  dicts = subset(dicts,(feature_id1%in%dict_keep$feature_id)&(feature_id2%in%dict_keep$feature_id))
  dicts$id1 = match(dicts$feature_id1,dict_keep$feature_id)
  dicts$id2 = match(dicts$feature_id2,dict_keep$feature_id)
  dicts = dicts %>% na.omit()
  
  within = Cos[cbind(dicts$id1,dicts$id2)]
  Cos[cbind(dicts$id1,dicts$id2)] = 0
  Cos[cbind(dicts$id2,dicts$id1)] = 0
  diag(Cos) = 0
  
  id1 = which(dict_keep$group==g1)
  id2 = which(dict_keep$group==g2)
  Cos = Cos[id1,id1]
  between = Cos[Cos!=0]
  return(mean(within)/mean(between))
}

Evaluate = function(cosU, LPrm = TRUE){
  mydict = data.frame(feature_id = rownames(cosU))
  mydict$group = sapply(mydict$feature_id,function(s) return(strsplit(s,":")[[1]][1]))
  if(LPrm){
    mydict$group[grep("LOINC:LP",mydict$feature_id)]="LOINC:LP"
  }
  
  lab.pairs = data.frame(feature_id1 = c(PR.pair$feature_id1,
                                         VA.pair$feature_id1,
                                         lab.pairs.LP$feature_id1),
                         feature_id2 = c(PR.pair$feature_id2,
                                         VA.pair$feature_id2,
                                         lab.pairs.LP$feature_id2))
  ## changing now ######################################
  res1 = ROC(cosU,mydict,rxnorm_rxnorm,g1='RXNORM',g2='RXNORM')
  res2 = ROC(cosU,mydict,PheCode_PheCode_int,g1='PheCode',g2='PheCode')
  res3 = ROC(cosU,mydict,lab.pairs,g1='LOINC',g2='LOINC')
  res4 = ROC(cosU,mydict,dict.wiki.keep,g1='PheCode',g2='PheCode')
  res5 = ROC(cosU,mydict,rbind(dict.online,Phecode_Rxnorm_Medrt_Snomed[,1:2]),g1='RXNORM',g2='PheCode')
  res6 = ROC(cosU,mydict,CCStoPhecode,g1='PheCode',g2='CCS')
  res7 = Lab_PheCode(cosU,mydict,rbind(PheCode.Lab.PR,PheCode.Lab.VA),g1='PheCode',g2=c('LOINC','ShortName'))
  
  Res = rbind('Drug-Drug'=res1,'Diease-Diease'=res2,'Lab-Lab'=res3,
              'Disease-Diease'=res4,'Disease-Drug'=res5,'Disease-Procedure'=res6,
              'Diease-Lab'=res7)
  
  Res = rbind(Res, 'similar' = colSums(Res[1:3,1]*Res[1:3,])/sum(Res[1:3,1]),
              'related' = colSums(Res[4:7,1]*Res[4:7,])/sum(Res[4:7,1]))
  Res[8,1] = sum(Res[1:3,1]); Res[9,1] = sum(Res[4:7,1])
  return(Res)
}

## functions for computing power for detecting known relationship pairs ####

get_type = function(AllRelationPairs){
  pairs = lapply(AllRelationPairs, function(x){ 
    x = x[,c('type',"pair","RELA","id")]
    x = x[!duplicated(x),]
  })
  pairs = do.call("rbind",pairs)
  pairs = pairs[!duplicated(pairs),]
  tn = lapply(1:28, function(i){
    if(sum(pairs$id==i)==0) return(c(type = "NA", name = "NA"))
    if(i<=9){
      name = unique(pairs$pair[which(pairs$id==i)])
      type = ifelse(i<9,unique(pairs$type[which(pairs$id==i)]),"rm")
      name = paste(name,"(",ifelse(i<=3,"sim",ifelse(i<=8,"rela","")),")",sep="")
    }else{
      name = ifelse(i<=27, pairs$RELA[which(pairs$id==i)], "CUI-code")
      type = unique(pairs$type[which(pairs$id==i)])
      name = paste(name,"(",ifelse(type=="similarity","sim","rela"),")",sep="")
    }
    stopifnot(length(type)==1&length(name)==1)
    return(c(type = type, name = name))
  })
  tn = do.call("rbind", tn)
  return(list(type = tn[,1], name = tn[,2]))
}


get_power_pmitest_sub = function(log_p_matrix, cosine_similarity, pairs, target_FDR){
  target_FDR = sort(unique(target_FDR), decreasing = TRUE)
  if(nrow(pairs)==0) return(list(ans = numeric(length(target_FDR)*3+2), 
                                 ptable = NULL))
  grouplist = data.frame(g1 = dict$group[match(pairs$code1, dict$code)],
                         g2 = dict$group[match(pairs$code2, dict$code)])
  group0 = grouplist[!duplicated(grouplist),]
  group = data.frame(g1 = apply(group0, 1, min),
                     g2 = apply(group0, 1, max))
  rm(group0)
  group = group[!duplicated(group),]
  codegroup = dict$group[match(rownames(log_p_matrix),dict$code)]
  pvalue_list = lapply(1:nrow(group), function(i){
    if(group[i,1]==group[i,2]){
      idx = which(codegroup==group[i,1])
      idlist = as.data.frame(t(combn(idx, 2)))
      colnames(idlist) = c("id1", "id2")
      idgrouplist = which(grouplist$g1==group[i,1] & grouplist$g2==group[i,2])
    }else{
      idx1 = which(codegroup==group[i,1])
      idx2 = which(codegroup==group[i,2])
      idlist = data.frame(id1 = rep(idx1, each = length(idx2)),
                          id2 = rep(idx2, length(idx1)))
      idgrouplist = which(grouplist$g1==group[i,1] & grouplist$g2==group[i,2])
      idgrouplist = unique(c(idgrouplist,
                             which(grouplist$g1==group[i,2] & grouplist$g2==group[i,1])))
    }
    grouplist_i = data.frame(id1 = match(pairs$code1[idgrouplist], rownames(log_p_matrix)),
                             id2 = match(pairs$code2[idgrouplist], rownames(log_p_matrix)))
    grouplist_i = data.frame(id1 = c(grouplist_i$id1, grouplist_i$id2),
                             id2 = c(grouplist_i$id2, grouplist_i$id1))
    idlist0 = setdiff(idlist, grouplist_i)
    idlist0$relation = 0
    idlist1 = intersect(idlist, grouplist_i)
    idlist1$relation = 1
    idlist = rbind(idlist1, idlist0)
    rm(idlist0, idlist1)
    idlist$logp = log_p_matrix[cbind(idlist$id1, idlist$id2)]
    idlist$cos = cosine_similarity[cbind(idlist$id1, idlist$id2)]
    return(idlist)
  })
  pvalue_list = do.call("rbind", pvalue_list)
  pvalue_list = pvalue_list[order(pvalue_list$logp),]
  pvalue_list$fdr = 1
  nmax = n = nrow(pvalue_list)
  ans = numeric()
  NRelation = sum(pvalue_list$relation)
  pvalue_list2 <- pvalue_list %>% arrange(desc(cos))
  for(fdr in target_FDR){
    nmax = max(c(which(pvalue_list$logp[1:nmax]<= log(fdr) + log(c(1:nmax) / n / (log(n)+1)))),0)
    pvalue_list$fdr[1:nmax] = fdr
    ans = c(ans, 
            sum(pvalue_list$relation[1:nmax])/NRelation)
    tmpid <- min(c(which(pvalue_list2$cos<=0), nrow(pvalue_list)))
    if(tmpid >= nmax){
      ans = c(ans, sum(pvalue_list2$relation[1:nmax])/NRelation)
    }else{
      tmpcount <- sum(pvalue_list2$relation[1:tmpid])
      tmpcount <- tmpcount + (NRelation - tmpcount)/(nrow(pvalue_list2)-tmpid)*(nmax-tmpid)
      ans = c(ans, tmpcount/NRelation)
    }
    ans = c(ans,nmax)
  }
  ans = c(ans, NRelation, n)
  names(ans) = c(paste(rep(c("power_test", "power_cos", "NSel"), length(target_FDR)),
                       rep(target_FDR, each = 2)),
                 "NPairs", "NHypo")
  return(list(ans = ans, ptable = pvalue_list))
}

Get_Power_PMItest = function(log_p_matrix, cosine_similarity, AllRelationPairs, target_FDR, echo = TRUE){
  tn = get_type(AllRelationPairs)
  fulltable = lapply(3:1, function(k){
    anslist = lapply(1:28, function(i){
      if(echo) cat("pair",k,"type",i,"\n")
      pairs = AllRelationPairs[[k]]
      pairs = pairs[which(pairs$id==i),]
      ans = get_power_pmitest_sub(log_p_matrix, cosine_similarity, pairs, target_FDR)
      if(echo) cat(ans$ans,"\n")
      return(ans$ans)
    })
    ans = do.call("rbind", anslist)
    ans = as.data.frame(ans)
    rownames(ans) = tn$name
    pairtype = tn$type
    id = which(pairtype=="similarity")
    ans = rbind(ans, apply(ans[id,],2,weighted.mean,ans[id,ncol(ans)-1]))
    ans[nrow(ans),ncol(ans)-1] = sum(ans[id,ncol(ans)-1])
    ans[nrow(ans),ncol(ans)] = sum(ans[id,ncol(ans)])
    rownames(ans)[nrow(ans)] = "weighted.sim"
    id = which(pairtype=="related")
    ans = rbind(ans, apply(ans[id,],2,weighted.mean,ans[id,ncol(ans)-1]))
    ans[nrow(ans),ncol(ans)-1] = sum(ans[id,ncol(ans)-1])
    ans[nrow(ans),ncol(ans)] = sum(ans[id,ncol(ans)])
    rownames(ans)[nrow(ans)] = "weighted.rela"
    ans = ans[c(which(pairtype=="similarity"),
                which(pairtype=="related"),
                which(rownames(ans)=="weighted.sim"),
                which(rownames(ans)=="weighted.rela")),]
    kable(ans,"latex",3)
    return(ans)
  })
  names(fulltable) = c("codi-codi","CUI-codi","CUI-CUI")
  return(fulltable)
}

Get_Power = function(log_p_matrix, cosine_similarity, AllRelationPairs, target_FDR){
  k = 3
  anslist = lapply(c(1,2,3,4,5,6,8), function(i){
    cat("pair",k,"type",i,"\n")
    pairs = AllRelationPairs[[k]]
    pairs = pairs[which(pairs$id==i),]
    ans = get_power_pmitest_sub(log_p_matrix, cosine_similarity, pairs, target_FDR)
    return(ans$ans)
  })
  return(do.call("rbind", anslist))
}