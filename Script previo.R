#library(readxl)
library(dplyr)
library(tidyverse)
library(EnvStats)


## Reading data

read.csv('Input files/resultados Placa TEG 29_09_23_dXPRESS.csv')->expr_data

expr_data<-expr_data%>%dplyr::filter(!(Name %in% c("PC 8","PC 9","PC 14",   
                                           "PC 22","PC 30","PC 45",
                                           "PC 63","775", "3014"))) #Excluding cases with neoadjuvancy and CECs
as.numeric(expr_data$Value)->expr_data$Value

#Subsetting housekeeping genes and showing expression levels

expr_data<-expr_data%>%dplyr::select(-ID)

# According with in-parallel analysis ran in NormFinder (REF), genes A and E were the most relevant to normalize data
# In this plot, their combination is shown as 'norm'

p<-expr_data%>%
  subset(Gene %in% c("A","B","C","D","E","F","G"))%>%
  pivot_wider(names_from=Gene,values_from=Value)%>%
  mutate(norm=(A+E)/2)%>% #Selected genes using normfinder
  pivot_longer(!c(Name),names_to="Gene",values_to="Value")%>%
  ggplot(aes(x=Gene,y=Value))+geom_boxplot(outlier.shape = NA)+geom_jitter(size=1,width = .15)+
  theme_bw()+ stat_mean_sd_text(y.pos = 25,size=3.5) 

png("Results/Housekeeping Profile.png",width = 6500,height = 2700,res=600)
p
dev.off()


dados_norm<-expr_data%>%
  pivot_wider(names_from=Gene,values_from=Value)%>%
  group_by(Name)%>%
  mutate(norm=mean(c(A,E)))%>%
  pivot_longer(!c(Name,norm),names_to="Gene",values_to="Value")%>%
  group_by(Name)%>%
  mutate(Value=norm-Value)%>%
  select(-norm)%>%
  pivot_wider(names_from=Gene,values_from=Value)%>%
  select(-A,-B,-C,-D,-E,-F,-G)
data.frame(dados_norm)->dados_norm
dados_norm$Name->rownames(dados_norm)
png("Resultados/heatmap_data_norm.png",width=4000,height = 4000,res=600)
heatmap(as.matrix(dados_norm[,-1]))
dev.off()
colnames(dados_norm)

names(tail(sort(colSums(is.na(dados_norm))/nrow(dados_norm)),10))->exclude_20

dados_norm<-dados_norm%>%
  select(-any_of(exclude_20))

dados_norm$Name->rownames(dados_norm)
names(tail(sort(colSums(is.na(t(dados_norm)))/ncol(dados_norm)),3))->exclude_20

dados_norm<-dados_norm%>%
  filter(!(Name %in% exclude_20))


png("Resultados/heatmap_data_norm_after.png",width=4000,height = 4000,res=600)
heatmap(as.matrix(dados_norm[,-1]))
dev.off()

names(tail(sort(colSums(is.na(dados_norm))/nrow(dados_norm)),5))->exclude_10

dados_norm<-dados_norm%>%
  select(-any_of(exclude_10))

dados_norm$Name->rownames(dados_norm)
colSums(is.na(t(dados_norm)))/ncol(dados_norm)->exclude_10

names(exclude_10[which(exclude_10!=0)])->exclude_10

dados_norm<-dados_norm%>%
  filter(!(Name %in% exclude_10))

