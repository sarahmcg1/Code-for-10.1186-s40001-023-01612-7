#load the libraries
library(tidyverse)

library(ggplot2)
library(ggrepel)

library(ggvenn)

library(flextable)


#From LASSO

load("tadaaPECAM.Rdata")
load("tadaasTM.Rdata")
load("tadaaDead.Rdata")

name_map <- read.csv("name_map.csv")

toptableDead <- left_join(tadaaDead,name_map,by = "HilleName") %>% mutate(t = abs(t))
toptableSTM <- left_join(tadaaSTM,name_map,by = "HilleName") %>% mutate(t = abs(t))
toptablePECAM <- left_join(tadaaPECAM,name_map,by = "HilleName") %>% mutate(t = abs(t))

pathways <- read.delim("2022-03-10_CPDB_pathways_metabolites.tab",header = TRUE,sep = "\t")

pathways_L <- separate_rows(pathways,metabolites,sep = ",")%>%
  separate(metabolites,c("type","ID"),sep=":")%>%
  filter(type == "kegg")%>%
  filter(source == c("KEGG"))%>%
  select(pathway,ID)%>%
  distinct()

degDead <- toptableDead %>%
  select(KEGG, t) %>% 
  filter(!is.na(t)) %>% 
  column_to_rownames(var = "KEGG") %>%
  as.matrix()

degSTM <- toptableSTM %>%
  select(KEGG, t) %>% 
  filter(!is.na(t)) %>% 
  column_to_rownames(var = "KEGG") %>%
  as.matrix()

degPECAM <- toptablePECAM %>%
  select(KEGG, t) %>% 
  filter(!is.na(t)) %>% 
  column_to_rownames(var = "KEGG") %>%
 as.matrix()

res_gseaDead <- run_ora(mat=degDead, net=pathways_L, .source='pathway', .target='ID',minsize = 3,
                      seed = 1)

res_gseaSTM <- run_ora(mat=degSTM, net=pathways_L, .source='pathway', .target='ID',minsize = 3,
                        seed = 1)
res_gseaPECAM <- run_ora(mat=degPECAM, net=pathways_L, .source='pathway', .target='ID',minsize = 3,
                          seed = 1)


f_contrast_actsDead <- res_gseaDead %>%
  filter(statistic == 'ora')%>%
  adjust_pvalue(method = "fdr",p.col = "p_value",output.col = "adjustedpvalue")%>%
  mutate(tag = ifelse(adjustedpvalue <= 0.05,source,NA)) %>% mutate(logP = -log10(adjustedpvalue)) %>% mutate(Comparision = "Death")

f_contrast_actsSTM <- res_gseaSTM %>%
  filter(statistic == 'ora')%>%
  adjust_pvalue(method = "fdr",p.col = "p_value",output.col = "adjustedpvalue")%>%
  mutate(tag = ifelse(adjustedpvalue <= 0.05,source,NA)) %>% mutate(logP = -log10(adjustedpvalue)) %>% mutate(Comparision = "sTM")

f_contrast_actsPECAM <- res_gseaPECAM %>%
  filter(statistic == 'ora')%>%
  adjust_pvalue(method = "fdr",p.col = "p_value",output.col = "adjustedpvalue")%>%
  mutate(tag = ifelse(adjustedpvalue <= 0.05,source,NA)) %>% mutate(logP = -log10(adjustedpvalue)) %>% mutate(Comparision = "PECAM")


ORAall <- f_contrast_actsDead %>% full_join(f_contrast_actsSTM) %>%
 full_join(f_contrast_actsPECAM) %>% filter(!is.na(tag)) %>% rename(KEGG = source)

#save(ORAall, file = "ORAall.Rdata")
#save(f_contrast_actsDead, file = "f_contrast_actsDead.Rdata")
#save(f_contrast_actsSTM, file = "f_contrast_actsSTM.Rdata")
#save(f_contrast_actsPECAM, file = "f_contrast_actsPECAM.Rdata")
#load("ORAall.Rdata")
#load("f_contrast_actsDead.Rdata")
#load("f_contrast_actsSTM.Rdata")
#load("f_contrast_actsPECAM.Rdata")

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
