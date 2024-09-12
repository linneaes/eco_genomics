setwd("/gpfs1/cl/pbio3990/Intro_to_R")
library(tidyverse)
Met<-read_csv("ThermalStressMetabolic.csv")
dim(Met)
head(Met)
tail(Met)
str(Met)
colnames(Met)
Met %>%
  rename(SMR = `SMR (g/hr)`)
Met_trans<-Met %>%
  rename(SMR = `SMR (g/hr)`,
         MMR = `MMR (g/hr)`)%>%
  mutate(ThermalScope= MMR-SMR)%>%
  select(Temp, Line, SMR, MMR)
Met_trans
Met_trans %>%
  filter(Temp == 34)
Met_trans %>%
  filter(Temp == 34 & Line =="LS1")
Met_trans %>%
  filter(Temp == 34 | Temp ==28)
Met_trans %>%
  filter(Temp %in% c(28, 34))
Met_trans %>%
  arrange(SMR)
Met_trans %>%
  arrange(desc(SMR))
Met_trans %>%
  distinct()
Met_trans %>%
  group_by(Temp, Line)
Met_trans %>%
  group_by(Temp) %>%
  summarise(
    avg_SMR = mean(SMR)
  )
Met_trans %>%
  group_by(Temp) %>%
  summarise(
    avg_SMR = mean(SMR, na.rm = T)
  )
Met_trans %>%
  group_by(Temp) %>%
  drop_na() %>%
  summarise(
    mean=mean(SMR),
    sd=sd(SMR),
    min=min(SMR),
    max=max(SMR),
    n=n()
  )
library(ggplot2)
Met_trans$Temp<-as.factor(Met_trans$Temp)
ggplot(Met_trans, aes(x=Temp, y=SMR))+
  geom_boxplot()
ggplot(Met_trans, aes(x=Temp, y=SMR, colour=Temp))+
  geom_boxplot()
ggplot(Met_trans, aes(x=Temp, y=SMR, fill=Temp))+
  geom_boxplot()
ggplot(Met_trans, aes(x=Temp, y=SMR, fill=Temp))+
  geom_boxplot()+
  scale_fill_manual(values=c("blue", "orange", "red"))
ggplot(Met_trans, aes(x=Temp, y=SMR, fill=Temp))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Dark2")
ggplot(Met_trans, aes(x=MMR, y=SMR, colour=Temp))+
  geom_point()+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~Line)+
  theme_bw()
