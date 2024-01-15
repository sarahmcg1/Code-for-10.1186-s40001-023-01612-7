#load the libraries
library(tidyverse)
library(glmnet)
library(DescTools)
library(ggplot2)
library(ggrepel)
library(siggitRausti)#from github
library(pROC)
library(readxl)
library(patchwork)


###set up the function roc_with
roc_with_ci <- function(obj){
  ciobj <- ci.se(obj, specificities = seq(0, 1, l = 25))
  dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                       lower = ciobj[, 1],
                       upper = ciobj[, 3])
  
  ggroc(obj,size=1) +
    theme_minimal() +
    geom_abline(
      slope = 1,
      intercept = 1,
      linetype = "dashed",
      alpha = 0.7,
      color = "grey10",
      size = 1) + 
    coord_equal() +
    geom_ribbon(
      data = dat.ci,
      aes(x = x, ymin = lower, ymax = upper),
      fill = "steelblue",
      alpha = 0.2) + 
    #ggtitle(capture.output(obj$ci)) + 
    ylab('Sensitivity') + xlab('Specificity') + 
    theme(plot.title   = element_text(size=20, hjust= .5, vjust = 2, face = "bold"),
          axis.text.x = element_text(face="bold", color="grey10", 
                                     size=12),
          axis.text.y = element_text(face="bold", color="grey10", 
                                     size=12),
          axis.title.x = element_text(face="bold", color="grey10", 
                                      size=18),
          axis.title.y = element_text(face="bold", color="grey10", 
                                      size=18), 
          #legend.title = element_blank(),
          #legend.spacing.y = unit(0, "mm"),
          panel.border = element_rect(colour = "white", fill=NA),
          aspect.ratio = 0.9, axis.text = element_text(colour = 1, size = 12),
          #legend.background = element_blank(),
          #legend.box.background = element_rect(colour = "black"),
          axis.line = element_line(colour ="black", size = 0.8),
          axis.ticks = element_line(colour ="black", size = 0.8),
          legend.text  = element_text(size =14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"))
} 

#load the data
load("hilleDATAall.Rdata")

#Sort out the data
hille <- select(hille,-c(Heigth,Weight,gcs,paoind,pafind,fioind,Vasopressoryesno,TotalfluidhoursbeforeICUarrival,Bilirubinmicromoll,CreatininemicromolL,PlateletsL))

#####Remove variables which are missing in 10% or more patients
missing_variables <- colnames(hille)[which(colSums(is.na(hille)) > (nrow(hille)/10))]
hille <- select(hille,which(colSums(is.na(hille)) <= (nrow(hille)/10)))


####Then remove any patient missing 10% of columns
missing_patients <- hille$MSOmicsID[(rowSums(is.na(hille)) > (ncol(hille)/10))]
hille <- filter(hille, (rowSums(is.na(hille)) <= (ncol(hille)/10)))

#select just mets then gap fill
mets <- select(hille,Alanine:ncol(hille))
mets <- select(mets,!Hypoxanthine)

mets <- missRanger::missRanger(mets,formula = .  ~ ., pmm.k = 3, verbose = FALSE)

#PAreto scale the metabolites

scaleP <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
metsS <- mets %>% mutate_all(log10) %>% mutate_all(scaleP)

load("hilleDATAall.Rdata")
fac <- select(hille,c(dead))
metsS <- rename(metsS, Glutamate = Glutamicacid,Tryptophan = Tryptophane, Aspartate = Asparticacid,Pyruvate = Pyruvicacid, Fumarate = Fumaricacid,Succinate = Succinicacid,Lactate = Lacticacid,Malate = Malicacid,
                Arachidonate = Arachidonicacid,DihomoLinolenicate = DihomoLinolenicacid,Docosahexanoate = Docosahexanoicacid,Docosapentanoate = Docosapentanoicacid,Eicosapentanoate = Eicosapentanoicacid,
                gammaLinolenicateA = gammaLinolenicacidalpha,Linoleate = Linoleicacid,Linolenate = Linolenicacid,Oleate = Oleicacid,Palmitate = Palmitoleicacid,Urate = Uricacid)

cmb <- cbind(fac,metsS)
cmb <- mutate_if(cmb,is.factor,~as.numeric(as.character(.x)))
# Now start the analysis (~ 60 seconds)
# Set the correct parameters:
repeat_no <- 500 # How often to bootstrap resample
alpha_val = 1 # Choose between ridge or lasso regression (0 = ridge, 1 = lasso, 0.5 = elastic net)
fraction_orig <- 1 # The fraction of original sample to work with. 
type_response = 'binomial'  # Take the value of 'gaussian' or 'binomial' depending on the analysis required. 
#Setup your data ( I call it houston regression. You can call it that too, so you dont need to change the 
# name in the rest of the code):

houston_regression <- data.frame(cmb)
library(siggitRausti)
# create empty vectors for variables (preallocating would be better here but eeeghhhhh)
lambda_vec_total <- c()
rmse_vec_total <- c()
rsq_vec_total <- c()
auc_vec_total <- c()
# Change the sampling matrix so that it is consistently sampled at 90/10 ratio alive/dead:
id_alive <- which(houston_regression$dead != 1)
id_dead <- which(houston_regression$dead == 1)
sample_matrix1 <- sample_vector(id_alive,nrow(houston_regression)*0.63,repeat_no,replicate = T)
sample_matrix2 <- sample_vector(id_dead,nrow(houston_regression)*0.37,repeat_no,replicate = T)
sample_matrix <- rbind(sample_matrix1,sample_matrix2)


for (n in 1:repeat_no){ # 
  houston_regression2 <- houston_regression[sample_matrix[,n],]
  if (type_response == 'binomial'){
    x <- as.matrix(houston_regression2[,-1])
    y <- houston_regression2$dead
    fraction_0<-sum(y==0)/length(y)
    fraction_1<-sum(y==1)/length(y)
    # assign that value to a "weights" vector
    weights<-rep(NA,length(y))
    weights[which(y == 0)] <- fraction_0
    weights[which(y == 1)] <- fraction_1
    penfac <- rep(1,ncol(x))
    penfac[c(1,2)] <- 1
    cv.out <- tryCatch({cv.glmnet(x,y,alpha =alpha_val,type.measure = "class",family = 'binomial',
                                  penalty.factor = penfac,weights = weights,nfolds = 5)}, error=function(e) e)
    if(inherits(cv.out, "error")) next
    lambda_lse <- cv.out$lambda.min
    lambda_vec_total <- c(lambda_vec_total,lambda_lse)
  } else if (type_response == 'gaussian') {
    #Use the test and train data partitions however you desire...
    x <- as.matrix(houston_regression2[,-1])
    y <- houston_regression2$dead
    cv.out <- cv.glmnet(x,y,alpha=alpha_val,type.measure = "mse",penalty.factor = penfac,nfolds = 5)
    lambda_lse <- cv.out$lambda.min
    lambda_vec_total <- c(lambda_vec_total,lambda_lse)
  }
}
CHOSEN_LAMBDA <- 10^median(log(lambda_vec_total,10)) # mean value of log-transformed distribution
if (type_response == 'binomial'){
  houston_regression_3 <- houston_regression[complete.cases(houston_regression),]
  x <- as.matrix(houston_regression_3[,-1])
  y <- houston_regression_3$dead
  fraction_0<-sum(y==0)/length(y)
  fraction_1<-sum(y==1)/length(y)
  # assign that value to a "weights" vector
  weights<-rep(NA,length(y))
  weights[which(y == 0)] <- fraction_0
  weights[which(y == 1)] <- fraction_1
  cv.fin <- glmnet(x,y,alpha=alpha_val,type.measure = "class",lambda = CHOSEN_LAMBDA,penalty.factor = penfac, weights = weights)
  heyhey  <- coef(cv.fin,s=CHOSEN_LAMBDA)
  inds<-which(abs(heyhey[,1])>0)
  variables<-row.names(heyhey)[inds]
  variables<-setdiff(variables,c('(Intercept)'))
  print(variables)
  modules_of_interest <- setdiff(variables,c('(Intercept)'))
  res <- predict(cv.fin,x,s = CHOSEN_LAMBDA,type = 'class') # Lambda 1se or lambda min. Choices choices choices...
  roccurve <- roc(houston_regression_3[,1] ~ as.numeric(res[,1]),quiet = TRUE,ci = T)
  auc_val <- pROC::auc(roccurve)
  houston_regression_4 <- houston_regression_3[,c(1,which(colnames(houston_regression_3) %in% variables))]
  houston_regression_4$dead <- as.factor(houston_regression_4$dead)
  #ROC_maker(houston_regression_4)
} else if (type_response == 'gaussian'){
  houston_regression_3 <- houston_regression[complete.cases(houston_regression),]
  x <- as.matrix(houston_regression_3[,-1])
  y <- houston_regression_3$dead
  cv.fin <- glmnet(x,y,alpha=alpha_val,type.measure = "mse",penalty.factor = penfac,lambda = median(lambda_vec_total))
  heyhey  <- coef(cv.fin,s=CHOSEN_LAMBDA)
  inds<-which(abs(heyhey[,1])>0)
  variables<-row.names(heyhey)[inds]
  variables<-setdiff(variables,c('(Intercept)'))
  print(variables)
  modules_of_interest <- setdiff(variables,c('(Intercept)'))
  res <- predict(cv.fin,x,s = CHOSEN_LAMBDA) 
  rmse_error <- sqrt(mean((houston_regression_3$dead - res)^2))
  rss <- sum((res - houston_regression_3$dead) ^ 2)  ## residual sum of squares
  tss <- sum((houston_regression_3$dead - mean(houston_regression_3$dead)) ^ 2)  ## total sum of squares
  rsq <- 1 - rss/tss
}
tadaaDead <- tibble(HilleName = c(names(heyhey[2:length(heyhey[ ,1]),1])),t = c(heyhey[2:length(heyhey[ ,1]),1]))
lambda_vec_totaldead <- lambda_vec_total
# Save variables that have positive and negative association with the response variable:
Death_pos <- names(heyhey[,1])[which(heyhey > 0)]
Death_neg <- names(heyhey[,1])[which(heyhey < 0)]
Death_pos <- Death_pos[which(Death_pos != "(Intercept)")]
Death_neg <- Death_neg[which(Death_neg != "(Intercept)")]
D <- roc_with_ci(roccurve)
De <- D + annotate(geom = 'text',x= 0.9, y= 0.9, label= print(paste0('AUC: ',round(roccurve$auc,2))))

De


####Redo for stM
fac <- select(hille,c(sTMngmL))

cmb <- cbind(fac,metsS)
cmb <- mutate_if(cmb,is.factor,~as.numeric(as.character(.x)))
# Now start the analysis (~ 60 seconds)


# Now start the analysis (~ 60 seconds)
# Set the correct parameters:
repeat_no <- 500 # How often to bootstrap resample
alpha_val = 1 # Choose between ridge or lasso regression (0 = ridge, 1 = lasso, 0.5 = elastic net)
fraction_orig <- 1 # The fraction of original sample to work with. 
type_response = 'gaussian'  # Take the value of 'gaussian' or 'binomial' depending on the analysis required. 
#Setup your data ( I call it houston regression. You can call it that too, so you dont need to change the 
# name in the rest of the code):

houston_regression <- data.frame(cmb)

# create empty vectors for variables (preallocating would be better here but eeeghhhhh)
lambda_vec_total <- c()
rmse_vec_total <- c()
rsq_vec_total <- c()
auc_vec_total <- c()
# Change the sampling matrix so that it is consistently sampled at 90/10 ratio alive/dead:
id_alive <- which(houston_regression$sTMngmL <= median(cmb$sTMngmL))
id_dead <- which(houston_regression$sTMngmL > median(cmb$sTMngmL))
sample_matrix1 <- sample_vector(id_alive,nrow(houston_regression)*0.63,repeat_no,replicate = T)
sample_matrix2 <- sample_vector(id_dead,nrow(houston_regression)*0.37,repeat_no,replicate = T)
sample_matrix <- rbind(sample_matrix1,sample_matrix2)


for (n in 1:repeat_no){ # 
  houston_regression2 <- houston_regression[sample_matrix[,n],]
  if (type_response == 'binomial'){
    x <- as.matrix(houston_regression2[,-1])
    y <- houston_regression2$sTMngmL
    fraction_0<-sum(y==0)/length(y)
    fraction_1<-sum(y==1)/length(y)
    # assign that value to a "weights" vector
    weights<-rep(NA,length(y))
    weights[which(y == 0)] <- fraction_0
    weights[which(y == 1)] <- fraction_1
    cv.out <- tryCatch({cv.glmnet(x,y,alpha =alpha_val,type.measure = "class",family = 'binomial',
                                  weights = weights,nfolds = 5)}, error=function(e) e)
    if(inherits(cv.out, "error")) next
    lambda_lse <- cv.out$lambda.min
    lambda_vec_total <- c(lambda_vec_total,lambda_lse)
  } else if (type_response == 'gaussian') {
    #Use the test and train data partitions however you desire...
    x <- as.matrix(houston_regression2[,-1])
    y <- houston_regression2$sTMngmL
    cv.out <- cv.glmnet(x,y,alpha=alpha_val,type.measure = "mse",nfolds = 5)
    lambda_lse <- cv.out$lambda.min
    lambda_vec_total <- c(lambda_vec_total,lambda_lse)
  }
}
CHOSEN_LAMBDA <- 10^median(log(lambda_vec_total,10)) # mean value of log-transformed distribution
if (type_response == 'binomial'){
  houston_regression_3 <- houston_regression[complete.cases(houston_regression),]
  x <- as.matrix(houston_regression_3[,-1])
  y <- houston_regression_3$sTMngmL
  fraction_0<-sum(y==0)/length(y)
  fraction_1<-sum(y==1)/length(y)
  # assign that value to a "weights" vector
  weights<-rep(NA,length(y))
  weights[which(y == 0)] <- fraction_0
  weights[which(y == 1)] <- fraction_1
  cv.fin <- glmnet(x,y,alpha=alpha_val,type.measure = "class",lambda = CHOSEN_LAMBDA, weights = weights)
  heyhey  <- coef(cv.fin,s=CHOSEN_LAMBDA)
  inds<-which(abs(heyhey[,1])>0)
  variables<-row.names(heyhey)[inds]
  variables<-setdiff(variables,c('(Intercept)'))
  print(variables)
  modules_of_interest <- setdiff(variables,c('(Intercept)'))
  res <- predict(cv.fin,x,s = CHOSEN_LAMBDA,type = 'class') # Lambda 1se or lambda min. Choices choices choices...
  roccurve <- roc(houston_regression_3[,1] ~ as.numeric(res[,1]),quiet = TRUE,ci = T)
  auc_val <- pROC::auc(roccurve)
  houston_regression_4 <- houston_regression_3[,c(1,which(colnames(houston_regression_3) %in% variables))]
  houston_regression_4$sTMngmL <- as.factor(houston_regression_4$sTMngmL)
  #ROC_maker(houston_regression_4)
} else if (type_response == 'gaussian'){
  houston_regression_3 <- houston_regression[complete.cases(houston_regression),]
  x <- as.matrix(houston_regression_3[,-1])
  y <- houston_regression_3$sTMngmL
  cv.fin <- glmnet(x,y,alpha=alpha_val,type.measure = "mse",lambda = median(lambda_vec_total))
  heyhey  <- coef(cv.fin,s=CHOSEN_LAMBDA)
  inds<-which(abs(heyhey[,1])>0)
  variables<-row.names(heyhey)[inds]
  variables<-setdiff(variables,c('(Intercept)'))
  print(variables)
  modules_of_interest <- setdiff(variables,c('(Intercept)'))
  res <- predict(cv.fin,x,s = CHOSEN_LAMBDA) 
  rmse_error <- sqrt(mean((houston_regression_3$sTMngmL - res)^2))
  rss <- sum((res - houston_regression_3$sTMngmL) ^ 2)  ## residual sum of squares
  tss <- sum((houston_regression_3$sTMngmL - mean(houston_regression_3$sTMngmL)) ^ 2)  ## total sum of squares
  rsq <- 1 - rss/tss
}
tadaasTM <- tibble(HilleName = c(names(heyhey[2:length(heyhey[ ,1]),1])),t = c(heyhey[2:length(heyhey[ ,1]),1]))
# Save variables that have positive and negative association with the response variable:
sTMngmL_pos <- names(heyhey[,1])[which(heyhey > 0)]
sTMngmL_neg <- names(heyhey[,1])[which(heyhey < 0)]
sTMngmL_pos <- sTMngmL_pos[which(sTMngmL_pos != "(Intercept)")]
sTMngmL_neg <- sTMngmL_neg[which(sTMngmL_neg != "(Intercept)")]

lambda_vec_totalstm <- lambda_vec_total

##Finally for PECAM
fac <- select(hille,c(PeCAMngmL))

cmb <- cbind(fac,metsS)
cmb <- mutate_if(cmb,is.factor,~as.numeric(as.character(.x)))

# Now start the analysis (~ 60 seconds)
# Set the correct parameters:
repeat_no <- 500 # How often to bootstrap resample
alpha_val = 1 # Choose between ridge or lasso regression (0 = ridge, 1 = lasso, 0.5 = elastic net)
fraction_orig <- 1 # The fraction of original sample to work with. 
type_response = 'gaussian'  # Take the value of 'gaussian' or 'binomial' depending on the analysis required. 
#Setup your data ( I call it houston regression. You can call it that too, so you dont need to change the 
# name in the rest of the code):

houston_regression <- data.frame(cmb)

# create empty vectors for variables (preallocating would be better here but eeeghhhhh)
lambda_vec_total <- c()
rmse_vec_total <- c()
rsq_vec_total <- c()
auc_vec_total <- c()
# Change the sampling matrix so that it is consistently sampled at 90/10 ratio alive/dead:
id_alive <- which(houston_regression$PeCAMngmL <= median(cmb$PeCAMngmL))
id_dead <- which(houston_regression$PeCAMngmL > median(cmb$PeCAMngmL))
sample_matrix1 <- sample_vector(id_alive,nrow(houston_regression)*0.63,repeat_no,replicate = T)
sample_matrix2 <- sample_vector(id_dead,nrow(houston_regression)*0.37,repeat_no,replicate = T)
sample_matrix <- rbind(sample_matrix1,sample_matrix2)


for (n in 1:repeat_no){ # 
  houston_regression2 <- houston_regression[sample_matrix[,n],]
  if (type_response == 'binomial'){
    x <- as.matrix(houston_regression2[,-1])
    y <- houston_regression2$PeCAMngmL
    fraction_0<-sum(y==0)/length(y)
    fraction_1<-sum(y==1)/length(y)
    # assign that value to a "weights" vector
    weights<-rep(NA,length(y))
    weights[which(y == 0)] <- fraction_0
    weights[which(y == 1)] <- fraction_1
    cv.out <- tryCatch({cv.glmnet(x,y,alpha =alpha_val,type.measure = "class",family = 'binomial',
                                  weights = weights,nfolds = 5)}, error=function(e) e)
    if(inherits(cv.out, "error")) next
    lambda_lse <- cv.out$lambda.min
    lambda_vec_total <- c(lambda_vec_total,lambda_lse)
  } else if (type_response == 'gaussian') {
    #Use the test and train data partitions however you desire...
    x <- as.matrix(houston_regression2[,-1])
    y <- houston_regression2$PeCAMngmL
    cv.out <- cv.glmnet(x,y,alpha=alpha_val,type.measure = "mse",nfolds = 5)
    lambda_lse <- cv.out$lambda.min
    lambda_vec_total <- c(lambda_vec_total,lambda_lse)
  }
}
CHOSEN_LAMBDA <- 10^median(log(lambda_vec_total,10)) # mean value of log-transformed distribution
if (type_response == 'binomial'){
  houston_regression_3 <- houston_regression[complete.cases(houston_regression),]
  x <- as.matrix(houston_regression_3[,-1])
  y <- houston_regression_3$PeCAMngmL
  fraction_0<-sum(y==0)/length(y)
  fraction_1<-sum(y==1)/length(y)
  # assign that value to a "weights" vector
  weights<-rep(NA,length(y))
  weights[which(y == 0)] <- fraction_0
  weights[which(y == 1)] <- fraction_1
  cv.fin <- glmnet(x,y,alpha=alpha_val,type.measure = "class",lambda = CHOSEN_LAMBDA, weights = weights)
  heyhey  <- coef(cv.fin,s=CHOSEN_LAMBDA)
  inds<-which(abs(heyhey[,1])>0)
  variables<-row.names(heyhey)[inds]
  variables<-setdiff(variables,c('(Intercept)'))
  print(variables)
  modules_of_interest <- setdiff(variables,c('(Intercept)'))
  res <- predict(cv.fin,x,s = CHOSEN_LAMBDA,type = 'class') # Lambda 1se or lambda min. Choices choices choices...
  roccurve <- roc(houston_regression_3[,1] ~ as.numeric(res[,1]),quiet = TRUE,ci = T)
  auc_val <- pROC::auc(roccurve)
  houston_regression_4 <- houston_regression_3[,c(1,which(colnames(houston_regression_3) %in% variables))]
  houston_regression_4$PeCAMngmL <- as.factor(houston_regression_4$dead)
  #ROC_maker(houston_regression_4)
} else if (type_response == 'gaussian'){
  houston_regression_3 <- houston_regression[complete.cases(houston_regression),]
  x <- as.matrix(houston_regression_3[,-1])
  y <- houston_regression_3$PeCAMngmL
  cv.fin <- glmnet(x,y,alpha=alpha_val,type.measure = "mse",lambda = median(lambda_vec_total))
  heyhey  <- coef(cv.fin,s=CHOSEN_LAMBDA)
  inds<-which(abs(heyhey[,1])>0)
  variables<-row.names(heyhey)[inds]
  variables<-setdiff(variables,c('(Intercept)'))
  print(variables)
  modules_of_interest <- setdiff(variables,c('(Intercept)'))
  res <- predict(cv.fin,x,s = CHOSEN_LAMBDA) 
  rmse_error <- sqrt(mean((houston_regression_3$PeCAMngmL - res)^2))
  rss <- sum((res - houston_regression_3$PeCAMngmL) ^ 2)  ## residual sum of squares
  tss <- sum((houston_regression_3$PeCAMngmL - mean(houston_regression_3$PeCAMngmL)) ^ 2)  ## total sum of squares
  rsq <- 1 - rss/tss
}
tadaaPECAM <- tibble(HilleName = c(names(heyhey[2:length(heyhey[ ,1]),1])),t = c(heyhey[2:length(heyhey[ ,1]),1]))
# Save variables that have positive and negative association with the response variable:
PeCAMngmL_pos <- names(heyhey[,1])[which(heyhey > 0)]
PeCAMngmL_neg <- names(heyhey[,1])[which(heyhey < 0)]

PeCAMngmL_pos <- PeCAMngmL_pos[which(PeCAMngmL_pos != "(Intercept)")]
PeCAMngmL_neg <- PeCAMngmL_neg[which(PeCAMngmL_neg != "(Intercept)")]

lambda_vec_totalpecam <- lambda_vec_total

up <- list(NonSurvival = c(Death_pos),sTM = c(sTMngmL_pos),PECAM = c(PeCAMngmL_pos))
down <- list(NonSurvival = c(Death_neg),sTM = c(sTMngmL_neg),PECAM = c(PeCAMngmL_neg))
