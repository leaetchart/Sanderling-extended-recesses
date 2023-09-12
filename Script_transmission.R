##Packages##

library(stringr)
library(tidyr)
library(lubridate)
library(tidyverse)
library(readxl)
library(car)
library(faraway)
library(lme4)
library(frair)
library(effsize)
library(MuMIn)
options(na.action = "na.fail")
library(boot)
library(lme4)
library(lmtest)
library(glmmTMB)
library(fitdistrplus)
library(scales)
library(DHARMa)
library(glmmTMB)
library(lmtest)
library(corrplot)
library(cowplot)
library(sjPlot)
library(ggeffects)
library(AICcmodavg)
library(ggResidpanel)
library(effects)
library(jtools)

####Recesses####
d<-read.table("Script & Datasets/Recesses.txt",sep="\t",h=T)
d<-d[d$recess_size_max>2,]

##1. Occurrence of extended recesses##
#Time window choice
df_MS <-d %>% drop_na(previous_temperature24)
g24 <- glmer(recesslong ~ previous_temperature24 + (1|nest_id), family = binomial, data=df_MS, nAGQ = 0)
g12 <- glmer(recesslong ~ previous_temperature12 + (1|nest_id), family = binomial, data=df_MS, nAGQ = 0)
g6 <- glmer(recesslong ~ previous_temperature6 + (1|nest_id), family = binomial, data=df_MS, nAGQ = 0)
g2 <- glmer(recesslong ~ previous_temperature2 + (1|nest_id), family = binomial, data=df_MS, nAGQ = 0)
g1 <- glmer(recesslong ~ previous_temperature1 + (1|nest_id), family = binomial, data=df_MS, nAGQ = 0)

models <- list(g1, g2, g6, g12, g24)
mod.names <- c('1h', '2h', '6h', '12h', '24h')
aictab(cand.set = models, modnames = mod.names, second.ord=FALSE)

#Analysis#
df_log <-d %>% drop_na(corporal_condition)
df_log <-df_log %>% drop_na(sexe)
df_log <-df_log %>% drop_na(nest_age)
df_log <-df_log %>% drop_na(previous_temperature12)
df_log <-df_log %>% drop_na(strategy)

#Model selection#
g_int <- glmer(recesslong ~ scale(corporal_condition) + sexe + scale(nest_age) + scale(previous_temperature12) + strategy + scale(corporal_condition):scale(previous_temperature12) + day + (1|nest_id), family = binomial, data=df_log, nAGQ = 0)
g_simple <- glmer(recesslong ~ scale(corporal_condition) + sexe + scale(nest_age) + scale(previous_temperature12) + strategy + day + (1|nest_id), family = binomial, data=df_log, nAGQ = 0)

r.squaredGLMM(g_int)
r.squaredGLMM(g_simple)

#Diagnostic plots#
simulationOutput <- simulateResiduals(fittedModel = g_simple)
plot(simulationOutput)
plotResiduals(simulationOutput, df_log$corporal_condition, quantreg = T)
plotResiduals(simulationOutput, df_log$sexe, quantreg = T)
plotResiduals(simulationOutput, df_log$nest_age, quantreg = T)
plotResiduals(simulationOutput, df_log$previous_temperature12, quantreg = T)
testDispersion(simulationOutput)
testResiduals(simulationOutput)

r_int<- ranef(g_simple)$nest_id$`(Intercept)`
qqnorm(r_int)
qqline(r_int)
shapiro.test(r_int)

#Plot#
predict<-ggpredict(g_simple, terms="previous_temperature12 [all]")
ggplot() +
  geom_line(data = predict, aes(x=x, y=predicted), colour="black")+
  geom_ribbon(data = predict, aes(y=predicted, ymin=conf.low, ymax=conf.high, x=x), alpha=0.4, fill="grey", colour="grey") +
  xlab("Ground-temperature (12h)") +
  ylab("Probability of performing an extended recess")+
  ylim(0,0.15)+
  xlim(0,27)+
  theme_classic()+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold"),
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")))+
  annotate("text", x=24, y=0.15, label= "R²m = 0.21", size=5)+
  annotate("text", x=24, y=0.14, label= "R²c = 0.43", size=5)


##2. Duration of extended recesses##
#Time window choice
df_MS2<- df_MS[df_MS$recesslong==1,]

m24<-lmerTest::lmer(recess_size_max ~previous_temperature24 + (1|nest_id), REML = FALSE, data=df_MS2)
m12<-lmerTest::lmer(recess_size_max ~ previous_temperature12 + (1|nest_id), REML = FALSE, data=df_MS2)
m6<-lmerTest::lmer(recess_size_max ~ previous_temperature6 + (1|nest_id), REML = FALSE, data=df_MS2)
m2<-lmerTest::lmer(recess_size_max ~ previous_temperature2 + (1|nest_id), REML = FALSE, data=df_MS2)
m1<-lmerTest::lmer(recess_size_max ~ previous_temperature1 + (1|nest_id), REML = FALSE, data=df_MS2)

models <- list(m24, m12, m6, m2, m1)
mod.names <- c('24', '12', '6', '2', '1')
aictab(cand.set = models, modnames = mod.names, second.ord=FALSE)

#Analysis#
df_RL<-d[d$recesslong==1,]
df_RL %>% group_by(nest_id) %>%tally()
df_RL <- df_RL %>% drop_na(previous_temperature6)
df_RL <- df_RL %>% drop_na(corporal_condition)
df_RL <- df_RL %>% drop_na(sexe)

#Model selection#
m_int<- lmerTest::lmer(recess_size_max ~ strategy + scale(nest_age) + sexe + scale(corporal_condition) + scale(previous_temperature6) + scale(corporal_condition):scale(previous_temperature6) + (1|nest_id), REML = TRUE, data=df_RL)

m_simple <- lmerTest::lmer(recess_size_max ~ strategy + scale(nest_age) + sexe + scale(corporal_condition) + scale(previous_temperature6) + (1|nest_id), REML = TRUE, data=df_RL)

r.squaredGLMM(m_int)
r.squaredGLMM(m_simple)
summary(m_int)
confint(m_int)

#Diagnostic plots#
res <- residuals(m_int)
pred <- predict(m_int)
resid_auxpanel(res, pred, plots = "default", bins = 30,
               smoother = FALSE, qqline = TRUE, qqbands = FALSE, scale = 1,
               theme = "bw", axis.text.size = 10, title.text.size = 12,
               title.opt = TRUE, nrow = NULL)

#Interaction#
interaction<- lmerTest::lmer(recess_size_max ~ previous_temperature6*corporal_condition + nest_age + sexe + strategy + (1|nest_id), REML = TRUE, data=df_RL)
interactions::sim_slopes(interaction, pred = previous_temperature6, modx = corporal_condition, modxvals = c(41.3, 57.8, 68.2), johnson_neyman = TRUE,
                         control.fdr = TRUE)
#Plot
predict_int<-ggpredict(m_int, terms=c("previous_temperature6 [all]", "corporal_condition [41.3, 57.8, 68.24]"))
predict_int <- predict_int %>% 
  mutate(Type = ifelse(group=="68.24", "ns", "p < 0.05"))
ggplot(predict_int, aes(x=x, y=predicted, group=group, linetype=Type))+
  scale_linetype_manual(values=c("dashed", "solid"))+
  geom_line(size=0.7)+
  xlab("Ground-temperature (6h)") +
  ylab("Duration of extended recesses (minutes)") +
  labs(col="Body condition")+
  geom_ribbon(data=predict_int, aes(y = predicted, ymin = conf.low, ymax = conf.high, fill= group), colour = "lightgrey", alpha=0.3)+
  theme_classic()+
  theme(axis.text=element_text(size=22),
        axis.title.y=element_text(size=22,face="bold"),
        axis.title.x=element_text(size=22,face="bold", vjust = -1),
        legend.title = element_text(color = "black", size = 16),
        legend.text = element_text(face = "plain", size = 14), 
        legend.spacing.y = unit(.5, 'cm'))+
  scale_fill_manual(values=c("hotpink2", "tan1", "cadetblue"), name="Body condition", labels=c('41.3', '57.8', '68.2'))+
  annotate("text", x=19, y=750, label= "R²m = 0.15", size=7)+
  annotate("text", x=19, y=700, label= "R²c = 0.23", size=7)+ 
  guides(linetype="none")

##3. Duration of short recesses##
df_RC<-d[d$recesslong !=1,]
df_RC<- df_RC %>% drop_na(corporal_condition)
df_RC<- df_RC %>% drop_na(sexe)
df_RC6<- df_RC %>% drop_na(previous_temperature6)

#Model selection#
m_int<-lmerTest::lmer(recess_size_max ~ strategy + scale(nest_age) + sexe + scale(corporal_condition) + scale(previous_temperature6) + scale(corporal_condition):scale(previous_temperature6) + (1|nest_id), REML = TRUE, data=df_RC6)
m_simple<-lmerTest::lmer(recess_size_max ~ strategy + scale(nest_age) + sexe + scale(corporal_condition) + scale(previous_temperature6) + (1|nest_id), REML = TRUE, data=df_RC6)
r.squaredGLMM(m_int)
r.squaredGLMM(m_simple)

summary(m_simple)
confint(m_simple)

#Diagnostic plots#
res <- residuals(m_int)
pred <- predict(m_int)
resid_auxpanel(res, pred, plots = "default", bins = 30,
               smoother = FALSE, qqline = TRUE, qqbands = FALSE, scale = 1,
               theme = "bw", axis.text.size = 10, title.text.size = 12,
               title.opt = TRUE, nrow = NULL)

#Plot#
predict_court<-ggpredict(m_simple, terms="previous_temperature6 [all]")
fitdata_6 = as.data.frame(Effect(c("previous_temperature6"),m_simple, xlevels=list(previous_temperature6=seq(0.5594935, 29.06393, 0.1))))
ggplot(data = df_RC6, aes(y=recess_size_max,x=previous_temperature6)) +
  geom_point(position = "jitter", alpha = 0.3)+
  geom_line(data = fitdata_6, aes(y=fit), colour="white", size=1)+
  geom_ribbon(data = fitdata_6, aes(y=fit, ymin=lower, ymax=upper), alpha=0.4, colour = "grey", fill="grey")+
  xlab("Ground-temperature (6h)") +
  ylab("Duration of short recesses (minutes)")+
  theme_classic()+
  theme(axis.text=element_text(size=22),
        axis.title.y=element_text(size=22,face="bold"),
        axis.title.x=element_text(size=22,face="bold", vjust = -1),
        legend.title = element_text(color = "black", size = 16),
        legend.text = element_text(face = "plain", size = 14), 
        legend.spacing.y = unit(.5, 'cm'))+
  annotate("text", x=25, y=150, label= "R²m = 0.002", size=7)+
  annotate("text", x=25, y=140, label= "R²c = 0.063", size=7)

####TDR####
rm(list=ls())
d2<-read.table("Script & Datasets/TDR.txt",sep="\t",h=TRUE)

##1. TDR for days without extended recesses##
df_courts<-d2[d2$TDR_long==0,]
df_courts <-df_courts %>% drop_na(temperature_during_24hperiod)
df_courts <- df_courts %>% drop_na(corporal_condition)
df_courts <- df_courts %>% drop_na(sexe)

#Model selection#
court_int<-lmerTest::lmer(TDR_begin ~ strategy + scale(nest_age) + sexe + scale(corporal_condition) + scale(temperature_during_24hperiod) + scale(corporal_condition):scale(temperature_during_24hperiod) + (1|nest_id), REML = TRUE, data=df_courts)
court_simple<-lmerTest::lmer(TDR_begin ~ strategy + scale(nest_age) + sexe + scale(corporal_condition) + scale(temperature_during_24hperiod) + (1|nest_id), REML = TRUE, data=df_courts)

r.squaredGLMM(court_int)
confint(court_int)
r.squaredGLMM(court_simple)
confint(court_simple)

#Diagnostic plots#
res <- residuals(court_simple)
pred <- predict(court_simple)
resid_auxpanel(res, pred, plots = "default", bins = 30,
               smoother = FALSE, qqline = TRUE, qqbands = FALSE, scale = 1,
               theme = "bw", axis.text.size = 10, title.text.size = 12,
               title.opt = TRUE, nrow = NULL)

#Plot#
fitdata_court = as.data.frame(Effect(c("temperature_during_24hperiod"),court_simple, xlevels=list(temperature_during_24hperiod=seq(2.653183, 19.24305, 0.1))))
ggplot(data = df_courts, aes(y=TDR_begin,x=temperature_during_24hperiod)) +
  geom_point(position = "jitter", alpha = 0.5)+
  geom_line(data = fitdata_court, aes(y=fit), colour="black", size=0.5)+
  geom_ribbon(data = fitdata_court, aes(y=fit, ymin=lower, ymax=upper), alpha=0.4, colour = "white", fill="grey")+
  xlab("Mean daily ground temperature (°C)") +
  ylab("TDR (minutes)
for days with
only short recesses")+
  ylim(0, 650)+
  xlim(2, 20)+
  theme_classic()+
  theme(axis.text=element_text(size=22),
        axis.title.y=element_text(size=22,face="bold"),
        axis.title.x=element_text(size=22,face="bold", vjust = -1),
        legend.title = element_text(color = "black", size = 16),
        legend.text = element_text(face = "plain", size = 14), 
        legend.spacing.y = unit(.5, 'cm'))+
  annotate("text", x=17, y=630, label= "R²m = 0.14", size=7)+
  annotate("text", x=17, y=580, label= "R²c = 0.57", size=7)

##2. TDR for days with extended recesses##
df_longs<-d2[d2$TDR_long>0,]
df_longs <-df_longs %>% drop_na(temperature_during_24hperiod)
df_longs <- df_longs %>% drop_na(corporal_condition)
df_longs <- df_longs %>% drop_na(sexe)

##Model selection##
long_int <- lmerTest::lmer(TDR_begin ~ strategy + scale(nest_age) + sexe + scale(corporal_condition) + scale(temperature_during_24hperiod)  + scale(corporal_condition):scale(temperature_during_24hperiod) + (1|nest_id), REML = TRUE, data=df_longs)
long_simple <- lmerTest::lmer(TDR_begin ~ strategy + scale(nest_age) + sexe + scale(corporal_condition) + scale(temperature_during_24hperiod) + (1|nest_id), REML = TRUE, data=df_longs)

r.squaredGLMM(long_int)
confint(long_int)
r.squaredGLMM(long_simple)
summary(long_int)

#Diagnostic plots#
res <- residuals(long_int)
pred <- predict(long_int)
resid_auxpanel(res, pred, plots = "default", bins = 30,
               smoother = FALSE, qqline = TRUE, qqbands = FALSE, scale = 1,
               theme = "bw", axis.text.size = 10, title.text.size = 12,
               title.opt = TRUE, nrow = NULL)

#Interaction#
int <- lmerTest::lmer(TDR_begin ~ strategy + nest_age + sexe + corporal_condition*temperature_during_24hperiod + (1|nest_id), REML = TRUE, data=df_longs)
interactions::sim_slopes(int, pred = temperature_during_24hperiod, modx = corporal_condition, modxvals = c(41.3, 57.8, 68.2),
                        johnson_neyman = TRUE, control.fdr = TRUE)

#Plot#
predict_TDRlong<-ggpredict(long_int, terms=c("temperature_during_24hperiod [all]", "corporal_condition [41.3, 57.8, 68.24]"))
predict_TDRlong <- predict_TDRlong %>% rename_at('group', ~'Body condition')
predict_TDRlong <- predict_TDRlong %>% 
  mutate(Type = ifelse(`Body condition`=="68.24", "ns", "p < 0.05"))

ggplot(predict_TDRlong, aes(x, predicted, fill=`Body condition`, linetype=Type)) +
  scale_linetype_manual(values=c("dashed", "solid"))+
  geom_line(size=0.7) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .25)+
  labs(x = "Mean daily ground temperature (°C)",
       y = "TDR (minutes)
for days with short
and extended recesses")+
  scale_y_continuous(breaks=seq(0,1100,200), limits=c(0,1100))+
  xlim(2,17)+
  theme_classic()+
  theme(axis.text=element_text(size=22),
        axis.title.x=element_text(size=22,face="bold"),
        axis.title.y = element_text( vjust=0.5, hjust=0.5, size=22,face="bold"),
        legend.title = element_text(color = "black", size = 16),
        legend.text = element_text(face = "plain", size = 14), 
        legend.spacing.y = unit(.5, 'cm'))+
  scale_fill_manual(values=c("hotpink2", "tan1", "cadetblue"), name="Body condition", labels=c('41.3', '57.8', '68.2'))+
  annotate("text", x=14, y=1100, label= "R²m = 0.23", col="black", size=7)+
  annotate("text", x=14, y=1000, label= "R²c = 0.47", col="black", size=7)+ guides(linetype="none")
