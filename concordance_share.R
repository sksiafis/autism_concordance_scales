library(readxl)
library(metafor)
library(tidyverse)
library(meta)
library(ggthemr)
library(RColorBrewer)
library(irr)
library(VGAM)
#R version 4.0.3, readxl_1.3.1, metafor_3.0-2, tidyverse_1.3.0, meta_5.1-1 , ggthemr_1.1.0, RColorBrewer_1.1-2, irr_0.84.1 

dust_theme <- ggthemr('dust', set_theme = FALSE)

#Dataset import
repetitive_pair<-read_xlsx("repetitive_pair_share.xlsx") #the dataset contains the SMDs for clinician and caregiver rating scales of restricted-repetitive behaviors in autism as measured in placebo-controlled pharmacological RCTs

#Creating datasets for the different pairs of clinician-caregiver scales
repetitive_pair<-repetitive_pair %>% filter(!is.na(clinic_te) & !is.na(careg_te)) #primary pairs based on the hierachy preferring YBOCS-versions for the clinicians and ABC-S for the caregivers
repetitive_pair_ybocs_abc<-repetitive_pair %>% filter(!is.na(abc_pair_te) & !is.na(cybocs_pair_te))
repetitive_pair_ybocs_rbs<-repetitive_pair %>% filter(!is.na(rbs_pair_te) & !is.na(cybocs_pair_te))
repetitive_pair_ados_rbs<-repetitive_pair %>% filter(!is.na(rbs_pair_te) & !is.na(ados_pair_te))

#1.Primary analysis

#1.a. intraclass correlation coefficients (ICC) using a two-way model and setting the unit of analysis as the mean of ratings 
icc(data.frame(caregiver=repetitive_pair$careg_te, clinician=repetitive_pair$clinic_te), 
    model="twoway", type="agreement", unit="average")


#1.b. meta-regression analysis using SMDclinician as the dependent and SMDcaregiver as the independent variable

reg_clinic_careg<-metareg(metagen(data=repetitive_pair, 
                                  TE=clinic_te, seTE=sqrt(clinic_vi)), ~careg_te)

reg_clinic_careg

ggplot(data=repetitive_pair, #meta-regression plot (with intercept)
       aes(x=careg_te, xmin=careg_te-1.96*sqrt(careg_vi), 
           xmax=careg_te+1.96*sqrt(careg_vi),
           y=clinic_te, ymax=clinic_te+1.96*sqrt(clinic_vi), 
           ymin=clinic_te-1.96*sqrt(clinic_vi)))+
  geom_point(colour="royalblue4", size=1/(sqrt(repetitive_pair$clinic_vi+repetitive_pair$careg_vi/2)))+
  geom_abline(intercept = reg_clinic_careg$beta[1,1], slope=reg_clinic_careg$beta[2,1], colour="royalblue4")+
  coord_cartesian(xlim=c(-0.2,1), ylim = c(-0.2,1))+
  geom_abline(intercept = 0, alpha=0.5, linetype="dashed")+
  geom_hline(yintercept = 0, alpha=1)+
  geom_vline(xintercept = 0, alpha=1)+
  scale_x_continuous(name="SMDcaregiver")+
  scale_y_continuous(name="SMDclinician")+
  dust_theme$theme


#1.3. meta-analysis of differences between SMDclinician and SMDcaregiver
repetitive_pair<-repetitive_pair %>% mutate(ri=0.3) %>% mutate(dg=clinic_te-careg_te, #calculate the difference assumming r=0.3
                                                               dg_se=sqrt(clinic_vi+careg_vi-2*((ri*(1/n1new+1/n2new))+((ri*ri*clinic_te*careg_te)/(2*(sample-2))))), #Cov formula from Gleser 2009
                                                               mean_dg=(clinic_te+careg_te)/2)

meta_pair<-metagen(studlab=plot_name, data=repetitive_pair, #generic-invariance model
                   TE=dg, seTE=dg_se, fixed = FALSE,prediction = TRUE) 

forest(meta_pair, sortvar = TE, prediction = TRUE, smlab="Dg",  layout="Revman5", #manual edits were done in the version of the manuscript 
       leftlabs = c("Study", "Dg", "seDg", "Weight", "Dg [95% CI]"))


#1.4. Bland-Altman 
s_prim<-shapiro.test(repetitive_pair$dg) #examining the normal distribution of the difference with a shapiro-wilk test

ggplot(data=repetitive_pair, aes(x=dg))+ #examining the normal distribituion with histogram and density plots
  geom_histogram(alpha=0.5)+
  geom_density(fill="royalblue4", alpha=0.5)+
  xlab("SMDclinician-SMDcaregiver")+
  ggtitle(paste0("Shapiro-Wilk test W=", round(s_prim$statistic, digits=2),", p-value=", round(s_prim$p.value, digits=2)))+
  dust_theme$theme


ggplot(data=repetitive_pair, aes(sample=dg))+ #examining the normal distribution with a qqplot
  stat_qq() + 
  stat_qq_line()+
  xlab("theoretical")+
  ylab("observed")+
  dust_theme$theme

metareg_dg_mean<-metareg(meta_pair,  ~mean_dg) #meta-regression to investigate the relationship between difference and mean of SMDclinician and SMDcaregiver


ggplot(data=repetitive_pair, aes(x=mean_dg, y=dg))+ #bland-altman plot
  geom_point(size=1/repetitive_pair$dg_se, colour="royalblue4", alpha=0.8)+
  ylab("SMDclinician-SMDcaregiver")+
  xlab("(SMDclinician+SMDcaregiver)/2")+
  geom_abline(slope=metareg_dg_mean$b[2,1], intercept = metareg_dg_mean$b[1,1], 
              size=1.2,colour="firebrick1")+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = meta_pair$TE.random, colour="royalblue4", size=1.2)+
  geom_hline(yintercept = meta_pair$TE.random - 1.96*meta_pair$seTE.random, linetype="dashed", colour="royalblue4", size=1.2)+
  geom_hline(yintercept = meta_pair$TE.random + 1.96*meta_pair$seTE.random, linetype="dashed", colour="royalblue4", size=1.2)+
  geom_hline(yintercept = meta_pair$lower.predict, linetype="dashed", colour="firebrick", size=1.2)+
  geom_hline(yintercept = meta_pair$upper.predict, linetype="dashed", colour="firebrick", size=1.2)+
  dust_theme$theme


#1.5 Funnel plot and Egger's test
funnel(meta_pair)

metabias(meta_pair)


#2. Sensitivity analyses using different assumptions of r
#2.1. r=0
repetitive_pair_0<-repetitive_pair %>% mutate(ri=0) %>% mutate(dg=clinic_te-careg_te,
                                                               dg_se=sqrt(clinic_vi+careg_vi-2*((ri*(1/n1new+1/n2new))+((ri*ri*clinic_te*careg_te)/(2*(sample-2))))), #Cov formula from Gleser 2009
                                                               mean_dg=(clinic_te+careg_te)/2)

metagen(studlab=plot_name, data=repetitive_pair_0, 
        TE=dg, seTE=dg_se, fixed = FALSE, sm="MD", prediction = TRUE)

#2.2. r=0.2
repetitive_pair_02<-repetitive_pair %>% mutate(ri=0.2) %>% mutate(dg=clinic_te-careg_te,
                                                                  dg_se=sqrt(clinic_vi+careg_vi-2*((ri*(1/n1new+1/n2new))+((ri*ri*clinic_te*careg_te)/(2*(sample-2))))), #Cov formula from Gleser 2009
                                                                  mean_dg=(clinic_te+careg_te)/2)

metagen(studlab=plot_name, data=repetitive_pair_02, 
         TE=dg, seTE=dg_se, fixed = FALSE, sm="MD", prediction = TRUE)

#2.3. r=0.5
repetitive_pair_05<-repetitive_pair %>% mutate(ri=0.5) %>% mutate(dg=clinic_te-careg_te,
                                                                  dg_se=sqrt(clinic_vi+careg_vi-2*((ri*(1/n1new+1/n2new))+((ri*ri*clinic_te*careg_te)/(2*(sample-2))))), #Cov formula from Gleser 2009
                                                                  mean_dg=(clinic_te+careg_te)/2)

metagen(studlab=plot_name, data=repetitive_pair_05, 
        TE=dg, seTE=dg_se, fixed = FALSE, sm="MD", prediction = TRUE)

#2.4. r=0.8
repetitive_pair_08<-repetitive_pair %>% mutate(ri=0.8) %>% mutate(dg=clinic_te-careg_te,
                                                                  dg_se=sqrt(clinic_vi+careg_vi-2*((ri*(1/n1new+1/n2new))+((ri*ri*clinic_te*careg_te)/(2*(sample-2))))), #Cov formula from Gleser 2009
                                                                  mean_dg=(clinic_te+careg_te)/2)

metagen(studlab=plot_name, data=repetitive_pair_08, 
                      TE=dg, seTE=dg_se, fixed = FALSE, sm="MD", prediction = TRUE)


#3. Subgroup analyses
#3.1. Pairs of rating scales: A study could have more than one pairs, so the effect sizes of each pair were calculated separately, then combined in a dataset, and subgroup analysis was conducted as below
repetitive_pair_ybocs_abc<-repetitive_pair_ybocs_abc %>% mutate(ri=0.3) %>% mutate(dg=cybocs_pair_te-abc_pair_te,
                                                                                   dg_se=sqrt(cybocs_pair_vi+abc_pair_vi-2*((ri*(1/n1new+1/n2new))+((ri*ri*cybocs_pair_te*abc_pair_te)/(2*(sample-2))))), #Cov formula from Gleser 2009
                                                                                   mean_dg=(cybocs_pair_te+abc_pair_te)/2)
meta_pair_ybocs_abc<-metagen(studlab=plot_name, data=repetitive_pair_ybocs_abc, 
                             TE=dg, seTE=dg_se, fixed = FALSE, sm="MD", prediction = TRUE)

repetitive_pair_ybocs_rbs<-repetitive_pair_ybocs_rbs %>% mutate(ri=0.3) %>%mutate(dg=cybocs_pair_te-rbs_pair_te,
                                                                                  dg_se=sqrt(cybocs_pair_vi+rbs_pair_vi-2*((ri*(1/n1new+1/n2new))+((ri*ri*cybocs_pair_te*rbs_pair_te)/(2*(sample-2))))), #Cov formula from Gleser 2009
                                                                                  mean_dg=(cybocs_pair_te+rbs_pair_te)/2)

meta_pair_ybocs_rbs<-metagen(studlab=plot_name, data=repetitive_pair_ybocs_rbs, 
                             TE=dg, seTE=dg_se, fixed = FALSE, sm="MD", prediction = TRUE)

repetitive_pair_ados_rbs<-repetitive_pair_ados_rbs %>% mutate(ri=0.3) %>% mutate(dg=ados_pair_te-rbs_pair_te,
                                                                                 dg_se=sqrt(ados_pair_vi+rbs_pair_vi-2*((ri*(1/n1new+1/n2new))+((ri*ri*ados_pair_te*rbs_pair_te)/(2*(sample-2))))), #Cov formula from Gleser 2009
                                                                                 mean_dg=(ados_pair_te+rbs_pair_te)/2)


meta_pair_ados_rbs<-metagen(studlab=plot_name, data=repetitive_pair_ados_rbs, 
                            TE=dg, seTE=dg_se, fixed = FALSE, sm="MD", prediction = TRUE)

repetitive_pair_ybocs_abc_2<-repetitive_pair_ybocs_abc %>% mutate(comparison="YBOCS vs. ABC-S") %>%
  select(study_name, dg, dg_se, comparison)
repetitive_pair_ybocs_rbs_2<-repetitive_pair_ybocs_rbs %>% mutate(comparison="YBOCS vs. RBS")%>%
  select(study_name, dg, dg_se, comparison)
repetitive_pair_ados_rbs_2<-repetitive_pair_ados_rbs %>% mutate(comparison="ADOS vs. RBS")%>%
  select(study_name, dg, dg_se, comparison)

repetitive_pair_all<-rbind(repetitive_pair_ybocs_abc_2, repetitive_pair_ybocs_rbs_2,
                           repetitive_pair_ados_rbs_2)

metagen(studlab=study_name, data=repetitive_pair_all, #subgroup analysis for the different scales (a few studies had more than one pairs)
                  TE=dg, seTE=dg_se, fixed = FALSE, sm="MD", prediction = TRUE,
                  subgroup = comparison)

#3.2. Risk of bias low (RobCategory==1) vs. unclear/high (RobCategory==2 or 3)
metagen(studlab=plot_name, data=repetitive_pair, 
        TE=dg, seTE=dg_se, fixed = FALSE, sm="MD", prediction = TRUE, subgroup = RoBCategory==1)

#3.3. Children/adolscents vs. adults
metagen(studlab=plot_name, data=repetitive_pair, 
                        TE=dg, seTE=dg_se, fixed = FALSE, sm="MD", prediction = TRUE,
                        subgroup = age_group)

#3.4. Drugs (antipsychotics==risperidone, aripiprazole or lurasidone, antidepressants/anxiolytics==citalopram, fluoxetine or buspirone, oxytocn, others==bumetanide, probiotics, guanfacine)
metagen(studlab=plot_name, data=repetitive_pair, 
                   TE=dg, seTE=dg_se, fixed = FALSE, sm="MD", prediction = TRUE,
                   subgroup = drug_category)

#4. Absolute Dgs based on the analyse-then-transform of Morrisey 2016

qfoldnorm(0.5, mean = meta_pair$TE.random, sd = sqrt(meta_pair$tau2)) #median of folded distribution
qfoldnorm(0.25, mean = meta_pair$TE.random, sd = sqrt(meta_pair$tau2)) #25% percentile
qfoldnorm(0.95, mean = meta_pair$TE.random, sd = sqrt(meta_pair$tau2)) #75% percentile



