#load the libraries
library(tidyverse)
#library(glmnet)
#library(DescTools)
library(ggplot2)
library(ggrepel)
#library(siggitRausti)#from github
#library(decoupleR)
#library(rstatix)
#library(pROC)
library(ggvenn)
#load the data
#load("hilleDATAall.Rdata")
library(flextable)

#Sort out the data
#hille <- select(hille,-c(Heigth,Weight,gcs,paoind,pafind,fioind,Vasopressoryesno,TotalfluidhoursbeforeICUarrival,Bilirubinmicromoll,CreatininemicromolL,PlateletsL))

#####Remove variables which are missing in 10% or more patients
#missing_variables <- colnames(hille)[which(colSums(is.na(hille)) > (nrow(hille)/10))]
#hille <- select(hille,which(colSums(is.na(hille)) <= (nrow(hille)/10)))


####Then remove any patient missing 10% of columns
#missing_patients <- hille$MSOmicsID[(rowSums(is.na(hille)) > (ncol(hille)/10))]
#hille <- filter(hille, (rowSums(is.na(hille)) <= (ncol(hille)/10)))

##select just mets then gap fill
#mets <- select(hille,Alanine:ncol(hille))
#mets <- select(mets,!Hypoxanthine)

#mets <- missRanger::missRanger(mets,formula = .  ~ ., pmm.k = 3, verbose = FALSE)

##PAreto scale the metabolites

#scaleP <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
#metsS <- mets %>% mutate_all(log10) %>% mutate_all(scaleP)

####set up the function roc_with
##library(pROC)
##library(caret)

#load("~/SepsisPaper/Nov22mets/hilleDATAall.Rdata")
#fac <- select(hille,c(dead))

#cmb <- cbind(fac,metsS)
#cmb <- mutate_if(cmb,is.factor,~as.numeric(as.character(.x)))
## Now start the analysis (~ 60 seconds)
## Set the correct parameters:
#repeat_no <- 500 # How often to bootstrap resample
#alpha_val = 1 # Choose between ridge or lasso regression (0 = ridge, 1 = lasso, 0.5 = elastic net)
#fraction_orig <- 1 # The fraction of original sample to work with. 
#type_response = 'binomial'  # Take the value of 'gaussian' or 'binomial' depending on the analysis required. 
##Setup your data ( I call it houston regression. You can call it that too, so you dont need to change the 
## name in the rest of the code):

#houston_regression <- data.frame(cmb)
#library(siggitRausti)
## create empty vectors for variables (preallocating would be better here but eeeghhhhh)
#lambda_vec_total <- c()
#rmse_vec_total <- c()
#rsq_vec_total <- c()
#auc_vec_total <- c()
## Change the sampling matrix so that it is consistently sampled at 90/10 ratio alive/dead:
#id_alive <- which(houston_regression$dead != 1)
#id_dead <- which(houston_regression$dead == 1)
#sample_matrix1 <- sample_vector(id_alive,nrow(houston_regression)*0.63,repeat_no,replicate = T)
#sample_matrix2 <- sample_vector(id_dead,nrow(houston_regression)*0.37,repeat_no,replicate = T)
#sample_matrix <- rbind(sample_matrix1,sample_matrix2)


#for (n in 1:repeat_no){ # 
#  houston_regression2 <- houston_regression[sample_matrix[,n],]
#  if (type_response == 'binomial'){
#    x <- as.matrix(houston_regression2[,-1])
#    y <- houston_regression2$dead
#    fraction_0<-sum(y==0)/length(y)
#    fraction_1<-sum(y==1)/length(y)
#    # assign that value to a "weights" vector
#    weights<-rep(NA,length(y))
#   weights[which(y == 0)] <- fraction_0
#    weights[which(y == 1)] <- fraction_1
#    penfac <- rep(1,ncol(x))
#    penfac[c(1,2)] <- 1
#    cv.out <- tryCatch({cv.glmnet(x,y,alpha =alpha_val,type.measure = "class",family = 'binomial',
#                                  penalty.factor = penfac,weights = weights,nfolds = 5)}, error=function(e) e)
#    if(inherits(cv.out, "error")) next
#    lambda_lse <- cv.out$lambda.min
#    lambda_vec_total <- c(lambda_vec_total,lambda_lse)
#  } else if (type_response == 'gaussian') {
#    #Use the test and train data partitions however you desire...
#    x <- as.matrix(houston_regression2[,-1])
#    y <- houston_regression2$dead
#    cv.out <- cv.glmnet(x,y,alpha=alpha_val,type.measure = "mse",penalty.factor = penfac,nfolds = 5)
#    lambda_lse <- cv.out$lambda.min
#    lambda_vec_total <- c(lambda_vec_total,lambda_lse)
#  }
#}
#CHOSEN_LAMBDA <- 10^median(log(lambda_vec_total,10)) # mean value of log-transformed distribution
#if (type_response == 'binomial'){
#  houston_regression_3 <- houston_regression[complete.cases(houston_regression),]
#  x <- as.matrix(houston_regression_3[,-1])
#  y <- houston_regression_3$dead
#  fraction_0<-sum(y==0)/length(y)
#  fraction_1<-sum(y==1)/length(y)
#  # assign that value to a "weights" vector
#  weights<-rep(NA,length(y))
#  weights[which(y == 0)] <- fraction_0
#  weights[which(y == 1)] <- fraction_1
#  cv.fin <- glmnet(x,y,alpha=alpha_val,type.measure = "class",lambda = CHOSEN_LAMBDA,penalty.factor = penfac, weights = weights)
#  heyhey  <- coef(cv.fin,s=CHOSEN_LAMBDA)
#  inds<-which(abs(heyhey[,1])>0)
#  variables<-row.names(heyhey)[inds]
#  variables<-setdiff(variables,c('(Intercept)'))
#  print(variables)
#  modules_of_interest <- setdiff(variables,c('(Intercept)'))
#  res <- predict(cv.fin,x,s = CHOSEN_LAMBDA,type = 'class') # Lambda 1se or lambda min. Choices choices choices...
#  roccurve <- roc(houston_regression_3[,1] ~ as.numeric(res[,1]),quiet = TRUE,ci = T)
#  auc_val <- pROC::auc(roccurve)
#  houston_regression_4 <- houston_regression_3[,c(1,which(colnames(houston_regression_3) %in% variables))]
#  houston_regression_4$dead <- as.factor(houston_regression_4$dead)
#  #ROC_maker(houston_regression_4)
#} else if (type_response == 'gaussian'){
#  houston_regression_3 <- houston_regression[complete.cases(houston_regression),]
#  x <- as.matrix(houston_regression_3[,-1])
#  y <- houston_regression_3$dead
#  cv.fin <- glmnet(x,y,alpha=alpha_val,type.measure = "mse",penalty.factor = penfac,lambda = median(lambda_vec_total))
#  heyhey  <- coef(cv.fin,s=CHOSEN_LAMBDA)
#  inds<-which(abs(heyhey[,1])>0)
#  variables<-row.names(heyhey)[inds]
#  variables<-setdiff(variables,c('(Intercept)'))
#  print(variables)
#  modules_of_interest <- setdiff(variables,c('(Intercept)'))
#  res <- predict(cv.fin,x,s = CHOSEN_LAMBDA) 
#  rmse_error <- sqrt(mean((houston_regression_3$dead - res)^2))
#  rss <- sum((res - houston_regression_3$dead) ^ 2)  ## residual sum of squares
#  tss <- sum((houston_regression_3$dead - mean(houston_regression_3$dead)) ^ 2)  ## total sum of squares
#  rsq <- 1 - rss/tss
#}

#tadaaDead <- tibble(HilleName = c(names(heyhey[2:length(heyhey[ ,1]),1])),t = c(heyhey[2:length(heyhey[ ,1]),1]))

#fac <- select(hille,c(sTMngmL))

#cmb <- cbind(fac,metsS)
#cmb <- mutate_if(cmb,is.factor,~as.numeric(as.character(.x)))
## Now start the analysis (~ 60 seconds)


## Now start the analysis (~ 60 seconds)
## Set the correct parameters:
#repeat_no <- 500 # How often to bootstrap resample
#alpha_val = 1 # Choose between ridge or lasso regression (0 = ridge, 1 = lasso, 0.5 = elastic net)
#fraction_orig <- 1 # The fraction of original sample to work with. 
#type_response = 'gaussian'  # Take the value of 'gaussian' or 'binomial' depending on the analysis required. 
##Setup your data ( I call it houston regression. You can call it that too, so you dont need to change the 
## name in the rest of the code):

#houston_regression <- data.frame(cmb)
#library(siggitRausti)
## create empty vectors for variables (preallocating would be better here but eeeghhhhh)
#lambda_vec_total <- c()
#rmse_vec_total <- c()
#rsq_vec_total <- c()
#auc_vec_total <- c()
## Change the sampling matrix so that it is consistently sampled at 90/10 ratio alive/dead:
#id_alive <- which(houston_regression$sTMngmL <= median(cmb$sTMngmL))
#id_dead <- which(houston_regression$sTMngmL > median(cmb$sTMngmL))
#sample_matrix1 <- sample_vector(id_alive,nrow(houston_regression)*0.63,repeat_no,replicate = T)
#sample_matrix2 <- sample_vector(id_dead,nrow(houston_regression)*0.37,repeat_no,replicate = T)
#sample_matrix <- rbind(sample_matrix1,sample_matrix2)


#for (n in 1:repeat_no){ # 
#  houston_regression2 <- houston_regression[sample_matrix[,n],]
#  if (type_response == 'binomial'){
#    x <- as.matrix(houston_regression2[,-1])
#    y <- houston_regression2$sTMngmL
#    fraction_0<-sum(y==0)/length(y)
#    fraction_1<-sum(y==1)/length(y)
#    # assign that value to a "weights" vector
#    weights<-rep(NA,length(y))
#    weights[which(y == 0)] <- fraction_0
#    weights[which(y == 1)] <- fraction_1
#    cv.out <- tryCatch({cv.glmnet(x,y,alpha =alpha_val,type.measure = "class",family = 'binomial',
#                                  weights = weights,nfolds = 5)}, error=function(e) e)
#    if(inherits(cv.out, "error")) next
#    lambda_lse <- cv.out$lambda.min
#    lambda_vec_total <- c(lambda_vec_total,lambda_lse)
#  } else if (type_response == 'gaussian') {
#    #Use the test and train data partitions however you desire...
#    x <- as.matrix(houston_regression2[,-1])
#    y <- houston_regression2$sTMngmL
#    cv.out <- cv.glmnet(x,y,alpha=alpha_val,type.measure = "mse",nfolds = 5)
#    lambda_lse <- cv.out$lambda.min
#    lambda_vec_total <- c(lambda_vec_total,lambda_lse)
#  }
#}
#CHOSEN_LAMBDA <- 10^median(log(lambda_vec_total,10)) # mean value of log-transformed distribution
#if (type_response == 'binomial'){
#  houston_regression_3 <- houston_regression[complete.cases(houston_regression),]
#  x <- as.matrix(houston_regression_3[,-1])
#  y <- houston_regression_3$sTMngmL
#  fraction_0<-sum(y==0)/length(y)
#  fraction_1<-sum(y==1)/length(y)
#  # assign that value to a "weights" vector
#  weights<-rep(NA,length(y))
#  weights[which(y == 0)] <- fraction_0
#  weights[which(y == 1)] <- fraction_1
#  cv.fin <- glmnet(x,y,alpha=alpha_val,type.measure = "class",lambda = CHOSEN_LAMBDA, weights = weights)
#  heyhey  <- coef(cv.fin,s=CHOSEN_LAMBDA)
#  inds<-which(abs(heyhey[,1])>0)
#  variables<-row.names(heyhey)[inds]
#  variables<-setdiff(variables,c('(Intercept)'))
#  print(variables)
#  modules_of_interest <- setdiff(variables,c('(Intercept)'))
#  res <- predict(cv.fin,x,s = CHOSEN_LAMBDA,type = 'class') # Lambda 1se or lambda min. Choices choices choices...
#  roccurve <- roc(houston_regression_3[,1] ~ as.numeric(res[,1]),quiet = TRUE,ci = T)
#  auc_val <- pROC::auc(roccurve)
#  houston_regression_4 <- houston_regression_3[,c(1,which(colnames(houston_regression_3) %in% variables))]
#  houston_regression_4$sTMngmL <- as.factor(houston_regression_4$sTMngmL)
#  #ROC_maker(houston_regression_4)
#} else if (type_response == 'gaussian'){
#  houston_regression_3 <- houston_regression[complete.cases(houston_regression),]
#  x <- as.matrix(houston_regression_3[,-1])
#  y <- houston_regression_3$sTMngmL
#  cv.fin <- glmnet(x,y,alpha=alpha_val,type.measure = "mse",lambda = median(lambda_vec_total))
#  heyhey  <- coef(cv.fin,s=CHOSEN_LAMBDA)
#  inds<-which(abs(heyhey[,1])>0)
#  variables<-row.names(heyhey)[inds]
#  variables<-setdiff(variables,c('(Intercept)'))
#  print(variables)
#  modules_of_interest <- setdiff(variables,c('(Intercept)'))
#  res <- predict(cv.fin,x,s = CHOSEN_LAMBDA) 
#  rmse_error <- sqrt(mean((houston_regression_3$sTMngmL - res)^2))
#  rss <- sum((res - houston_regression_3$sTMngmL) ^ 2)  ## residual sum of squares
#  tss <- sum((houston_regression_3$sTMngmL - mean(houston_regression_3$sTMngmL)) ^ 2)  ## total sum of squares
#  rsq <- 1 - rss/tss
#}

#tadaaSTM <- tibble(HilleName = c(names(heyhey[2:length(heyhey[ ,1]),1])),t = c(heyhey[2:length(heyhey[ ,1]),1]))

#fac <- select(hille,c(PeCAMngmL))

#cmb <- cbind(fac,metsS)
#cmb <- mutate_if(cmb,is.factor,~as.numeric(as.character(.x)))

## Now start the analysis (~ 60 seconds)
## Set the correct parameters:
#repeat_no <- 500 # How often to bootstrap resample
#alpha_val = 1 # Choose between ridge or lasso regression (0 = ridge, 1 = lasso, 0.5 = elastic net)
#fraction_orig <- 1 # The fraction of original sample to work with. 
#type_response = 'gaussian'  # Take the value of 'gaussian' or 'binomial' depending on the analysis required. 
##Setup your data ( I call it houston regression. You can call it that too, so you dont need to change the 
## name in the rest of the code):

#houston_regression <- data.frame(cmb)
#library(siggitRausti)
## create empty vectors for variables (preallocating would be better here but eeeghhhhh)
#lambda_vec_total <- c()
#rmse_vec_total <- c()
#rsq_vec_total <- c()
#auc_vec_total <- c()
## Change the sampling matrix so that it is consistently sampled at 90/10 ratio alive/dead:
#id_alive <- which(houston_regression$PeCAMngmL <= median(cmb$PeCAMngmL))
#id_dead <- which(houston_regression$PeCAMngmL > median(cmb$PeCAMngmL))
#sample_matrix1 <- sample_vector(id_alive,nrow(houston_regression)*0.63,repeat_no,replicate = T)
#sample_matrix2 <- sample_vector(id_dead,nrow(houston_regression)*0.37,repeat_no,replicate = T)
#sample_matrix <- rbind(sample_matrix1,sample_matrix2)


#for (n in 1:repeat_no){ # 
#  houston_regression2 <- houston_regression[sample_matrix[,n],]
#  if (type_response == 'binomial'){
#    x <- as.matrix(houston_regression2[,-1])
#    y <- houston_regression2$PeCAMngmL
#    fraction_0<-sum(y==0)/length(y)
#    fraction_1<-sum(y==1)/length(y)
#    # assign that value to a "weights" vector
#    weights<-rep(NA,length(y))
#    weights[which(y == 0)] <- fraction_0
#    weights[which(y == 1)] <- fraction_1
#    cv.out <- tryCatch({cv.glmnet(x,y,alpha =alpha_val,type.measure = "class",family = 'binomial',
#                                  weights = weights,nfolds = 5)}, error=function(e) e)
#    if(inherits(cv.out, "error")) next
#    lambda_lse <- cv.out$lambda.min
#    lambda_vec_total <- c(lambda_vec_total,lambda_lse)
#  } else if (type_response == 'gaussian') {
#    #Use the test and train data partitions however you desire...
#    x <- as.matrix(houston_regression2[,-1])
#    y <- houston_regression2$PeCAMngmL
#    cv.out <- cv.glmnet(x,y,alpha=alpha_val,type.measure = "mse",nfolds = 5)
#    lambda_lse <- cv.out$lambda.min
#    lambda_vec_total <- c(lambda_vec_total,lambda_lse)
#  }
#}
#CHOSEN_LAMBDA <- 10^median(log(lambda_vec_total,10)) # mean value of log-transformed distribution
#if (type_response == 'binomial'){
#  houston_regression_3 <- houston_regression[complete.cases(houston_regression),]
#  x <- as.matrix(houston_regression_3[,-1])
#  y <- houston_regression_3$PeCAMngmL
#  fraction_0<-sum(y==0)/length(y)
#  fraction_1<-sum(y==1)/length(y)
#  # assign that value to a "weights" vector
#  weights<-rep(NA,length(y))
#  weights[which(y == 0)] <- fraction_0
#  weights[which(y == 1)] <- fraction_1
# cv.fin <- glmnet(x,y,alpha=alpha_val,type.measure = "class",lambda = CHOSEN_LAMBDA, weights = weights)
#  heyhey  <- coef(cv.fin,s=CHOSEN_LAMBDA)
#  inds<-which(abs(heyhey[,1])>0)
#  variables<-row.names(heyhey)[inds]
#  variables<-setdiff(variables,c('(Intercept)'))
#  print(variables)
#  modules_of_interest <- setdiff(variables,c('(Intercept)'))
#  res <- predict(cv.fin,x,s = CHOSEN_LAMBDA,type = 'class') # Lambda 1se or lambda min. Choices choices choices...
#  roccurve <- roc(houston_regression_3[,1] ~ as.numeric(res[,1]),quiet = TRUE,ci = T)
#  auc_val <- pROC::auc(roccurve)
#  houston_regression_4 <- houston_regression_3[,c(1,which(colnames(houston_regression_3) %in% variables))]
#  houston_regression_4$PeCAMngmL <- as.factor(houston_regression_4$dead)
#  #ROC_maker(houston_regression_4)
#} else if (type_response == 'gaussian'){
#  houston_regression_3 <- houston_regression[complete.cases(houston_regression),]
#  x <- as.matrix(houston_regression_3[,-1])
#  y <- houston_regression_3$PeCAMngmL
#  cv.fin <- glmnet(x,y,alpha=alpha_val,type.measure = "mse",lambda = median(lambda_vec_total))
#  heyhey  <- coef(cv.fin,s=CHOSEN_LAMBDA)
#  inds<-which(abs(heyhey[,1])>0)
#  variables<-row.names(heyhey)[inds]
#  variables<-setdiff(variables,c('(Intercept)'))
#  print(variables)
#  modules_of_interest <- setdiff(variables,c('(Intercept)'))
#  res <- predict(cv.fin,x,s = CHOSEN_LAMBDA) 
#  rmse_error <- sqrt(mean((houston_regression_3$PeCAMngmL - res)^2))
#  rss <- sum((res - houston_regression_3$PeCAMngmL) ^ 2)  ## residual sum of squares
#  tss <- sum((houston_regression_3$PeCAMngmL - mean(houston_regression_3$PeCAMngmL)) ^ 2)  ## total sum of squares
#  rsq <- 1 - rss/tss
#}

#tadaaPECAM <- tibble(HilleName = c(names(heyhey[2:length(heyhey[ ,1]),1])),t = c(heyhey[2:length(heyhey[ ,1]),1]))


#name_map <- read.csv("name_map.csv")

#toptableDead <- left_join(tadaaDead,name_map,by = "HilleName") %>% mutate(t = abs(t))
#toptableSTM <- left_join(tadaaSTM,name_map,by = "HilleName") %>% mutate(t = abs(t))
#toptablePECAM <- left_join(tadaaPECAM,name_map,by = "HilleName") %>% mutate(t = abs(t))

#pathways <- read.delim("2022-03-10_CPDB_pathways_metabolites.tab",header = TRUE,sep = "\t")

#pathways_L <- separate_rows(pathways,metabolites,sep = ",")%>%
#  separate(metabolites,c("type","ID"),sep=":")%>%
#  filter(type == "kegg")%>%
#  filter(source == c("KEGG"))%>%
#  select(pathway,ID)%>%
#  distinct()

#degDead <- toptableDead %>%
#  select(KEGG, t) %>% 
#  filter(!is.na(t)) %>% 
#  column_to_rownames(var = "KEGG") %>%
#  as.matrix()

#degSTM <- toptableSTM %>%
#  select(KEGG, t) %>% 
#  filter(!is.na(t)) %>% 
#  column_to_rownames(var = "KEGG") %>%
#  as.matrix()

#degPECAM <- toptablePECAM %>%
#  select(KEGG, t) %>% 
#  filter(!is.na(t)) %>% 
#  column_to_rownames(var = "KEGG") %>%
#  as.matrix()

#res_gseaDead <- run_ora(mat=degDead, net=pathways_L, .source='pathway', .target='ID',minsize = 3,
#                      seed = 1)

#res_gseaSTM <- run_ora(mat=degSTM, net=pathways_L, .source='pathway', .target='ID',minsize = 3,
#                        seed = 1)

#res_gseaPECAM <- run_ora(mat=degPECAM, net=pathways_L, .source='pathway', .target='ID',minsize = 3,
#                          seed = 1)


#f_contrast_actsDead <- res_gseaDead %>%
#  filter(statistic == 'ora')%>%
#  adjust_pvalue(method = "fdr",p.col = "p_value",output.col = "adjustedpvalue")%>%
#  mutate(tag = ifelse(adjustedpvalue <= 0.05,source,NA)) %>% mutate(logP = -log10(adjustedpvalue)) %>% mutate(Comparision = "Death")

#f_contrast_actsSTM <- res_gseaSTM %>%
#  filter(statistic == 'ora')%>%
#  adjust_pvalue(method = "fdr",p.col = "p_value",output.col = "adjustedpvalue")%>%
#  mutate(tag = ifelse(adjustedpvalue <= 0.05,source,NA)) %>% mutate(logP = -log10(adjustedpvalue)) %>% mutate(Comparision = "sTM")

#f_contrast_actsPECAM <- res_gseaPECAM %>%
#  filter(statistic == 'ora')%>%
#  adjust_pvalue(method = "fdr",p.col = "p_value",output.col = "adjustedpvalue")%>%
#  mutate(tag = ifelse(adjustedpvalue <= 0.05,source,NA)) %>% mutate(logP = -log10(adjustedpvalue)) %>% mutate(Comparision = "PECAM")


#ORAall <- f_contrast_actsDead %>% full_join(f_contrast_actsSTM) %>%
#  full_join(f_contrast_actsPECAM) %>% filter(!is.na(tag)) %>% rename(KEGG = source)

#save(ORAall, file = "ORAall.Rdata")
#save(f_contrast_actsDead, file = "f_contrast_actsDead.Rdata")
#save(f_contrast_actsSTM, file = "f_contrast_actsSTM.Rdata")
#save(f_contrast_actsPECAM, file = "f_contrast_actsPECAM.Rdata")
load("ORAall.Rdata")
load("f_contrast_actsDead.Rdata")
load("f_contrast_actsSTM.Rdata")
load("f_contrast_actsPECAM.Rdata")

ORAall$Comparision <- str_replace_all(ORAall$Comparision,"Death","Non-Survival")

# Plot
pathspoints <- ggplot(ORAall, aes(x=KEGG, y=logP, color = Comparision,size = score)) +
  geom_point( stat="identity")+ theme(axis.text.x = element_text(size = 6,angle=60, hjust = 1),axis.text.y = element_text(size = 6, hjust = 1),legend.position="top",legend.title = element_text(size=6))+
  scale_color_manual(values=c("#0073C2FF", "#EFC000FF", "#CD534CFF"))+
  coord_flip()


PathsDead <- unique(f_contrast_actsDead$tag)
PathsDead <- PathsDead[which(PathsDead != "NA")]
Pathsstm <- unique(f_contrast_actsSTM$tag)
Pathsstm <- Pathsstm[which(Pathsstm != "NA")]
PathsPECAM <- unique(f_contrast_actsPECAM$tag)
PathsPECAM <- PathsPECAM[which(PathsPECAM != "NA")]

paths <- list(NonSurvival = c(PathsDead),sTM = c(Pathsstm),PECAM = c(PathsPECAM))
pathsvenn <- ggvenn(
  paths, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

allpaths <- intersect(intersect(PathsDead,Pathsstm),PathsPECAM)
endopaths <- intersect(Pathsstm,PathsPECAM)
stmdead <- intersect(Pathsstm,PathsDead)
pecamdead <- intersect(PathsDead,PathsPECAM)
onlysTM <- setdiff(Pathsstm,PathsDead)
onlystM <- setdiff(onlysTM,PathsPECAM)
onlyDead <- setdiff(PathsDead,Pathsstm)
onlyDead <- setdiff(onlyDead,PathsPECAM)
onlyPECAM <- setdiff(PathsPECAM,PathsDead)
onlyPECAM <- setdiff(onlyPECAM,Pathsstm)

ORAtbl <- tibble(ORAall) %>% select(-c("statistic","condition","p_value","tag","logP")) %>%
  relocate(Comparision, .before = KEGG) %>% mutate(score = round(score,3),adjustedpvalue = round(adjustedpvalue,3))

###### make a pretty flex table
myftpath <- flextable(ORAtbl)
myftpath <- theme_vanilla(myftpath)
myftpath <- merge_v(myftpath, j = "Comparision")
myftpath <- italic(myftpath,j = "Comparision", italic = TRUE)
myftpath <- labelizor(
  x = myftpath, 
  part = "header", 
  labels = c("Comparision" = "Comparison","score" = "ORA Score","adjustedpvalue" = "Adjusted P-value"))
myftpath <- width(myftpath, j = c("Comparision","KEGG","score","adjustedpvalue"), c(1.25,2,0.75,0.85))
