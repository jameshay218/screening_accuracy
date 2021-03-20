#############################################################
## Author: James Hay
## Date: 19th March 2021
## Contact: jhay@hsph.harvard.edu
## Description: an agent-based model simulation of SARS-CoV-2 sensitivity, specificity and positive predictive value
##              over a short period of epidemic decling with decreasing incidence and prevalence.

library(tidyverse)
library(patchwork)
setwd("C:/Users/James//Google Drive/nCoV/testing_stats_sim/")

## Plot settings
lft_col <- "#E37222FF"
pcr_col <-"#2A6EBBFF"

linelist_color_key <- c("False negative"="#7D5CC6", 
                        "Not detected"="#CD202CFF", 
                        "True positive"="#69BE28FF", 
                        "Post isolation"="#7FA8D6FF", 
                        "Recovered"="#7FA8D6FF",
                        "False positive"="#F0AB00FF",
                        "Correctly\n isolating"="#00B2A9", 
                        "Incorrectly\n isolating"="black", 
                        "True negative"="#ACADAEFF",
                        "True negative/\nRecovered"="#ACADAEFF",
                        "True negative/\nPost infection"="#ACADAEFF",
                        "Not tested"="white")

status_levels <- c("True positive","False positive","True negative","False negative","Recovered",
                   "Post isolation","Correctly\n isolating","Incorrectly\n isolating","Not detected",
                   "True negative/\nPost infection","True negative/\nRecovered",
                   "Not tested")

theme_overall <- theme(axis.text.x=element_text(size=6),
                       axis.text.y=element_text(size=6),
                       axis.title.x=element_text(size=8),
                       axis.title.y=element_text(size=8),
                       legend.text=element_text(size=6),
                       legend.title=element_text(size=6),
                       plot.title=element_text(size=8,face="bold"),
                       plot.margin=unit(c(0,0,0,0),"cm")
                       )

############################################
## 1. Simulation control parameters
############################################
## Simulate individual-level line list as well as cohort
simulate_linelist <- TRUE
## =====================
## Epidemic control
## =====================
growth_rate <- -0.05
## Simulation duration
## Cumulative incidence in the focal period (ie. between burnin and max(times))
cumu_inc <- 0.01
## Population size
N <- 3000000
## For how long are people true infections?
duration_pos <- 21
## Burning for simulation
burnin <- duration_pos + 25
sim_duration <- 55
times <- 1:(burnin + sim_duration)
t_choose <- floor((max(times)+burnin)/2)

## =====================
## PCR surveillance
## =====================
## PCR sensitivity
pcr_gamma_mode <- 5
pcr_gamma_sd <- 10
pcr_gamma_height <- 0.95
pcr_gamma_max_sens <- 0.95

pcr_gamma_rate = ( pcr_gamma_mode + sqrt( pcr_gamma_mode^2 + 4*pcr_gamma_sd^2 ) ) / ( 2 * pcr_gamma_sd^2 )
pcr_gamma_shape = 1 + pcr_gamma_mode * pcr_gamma_rate
detect_prob_pcr <- dgamma(0:duration_pos,shape = pcr_gamma_shape,rate=pcr_gamma_rate)
detect_prob_pcr <- detect_prob_pcr*(pcr_gamma_height/max(detect_prob_pcr))
detect_prob_pcr <- pmin(detect_prob_pcr, pcr_gamma_max_sens)

## =====================
## LFT strategy
## =====================
## LFT sensitivity
lft_gamma_mode <- pcr_gamma_mode ## Assume peaks at same time
lft_gamma_sd <- 5
lft_gamma_height <- 1
lft_gamma_max_sens <- 0.95

lft_gamma_rate = ( lft_gamma_mode + sqrt( lft_gamma_mode^2 + 4*lft_gamma_sd^2 ) ) / ( 2 * lft_gamma_sd^2 )
lft_gamma_shape = 1 + lft_gamma_mode * lft_gamma_rate
detect_prob_lft <- dgamma(0:duration_pos,shape = lft_gamma_shape,rate=lft_gamma_rate)
detect_prob_lft <- detect_prob_lft*(lft_gamma_height/max(detect_prob_lft))
detect_prob_lft <- pmin(detect_prob_lft, lft_gamma_max_sens)

## =====================
## Specificity and test frequency
## =====================
## PCR spec
pcr_spec <- 0.999
## LFT spec
lft_spec <- 0.999

## LFT screening approach
test_frequency <- 3
isolation_period <- 14

## =====================
## Time-varying maximum sensitivity over infection ie. within-host model
## =====================
max_sens <- tibble(days_post_inf=0:duration_pos, 
                   detect_prob_pcr=detect_prob_pcr,
                   detect_prob_lft=detect_prob_lft)


print(mean(detect_prob_pcr[1:duration_pos]))
print(mean(detect_prob_lft[1:duration_pos]))

p_test_sens <- max_sens %>% 
  rename(`qPCR`=detect_prob_pcr,`LFT`=detect_prob_lft) %>%
  pivot_longer(-days_post_inf) %>%
  rename(Test=name) %>%
  ggplot() +
  geom_line(aes(x=days_post_inf,y=value,col=Test),size=1) +
  geom_vline(xintercept=5.5,linetype="dashed") +
  scale_color_manual(values=c("qPCR"=pcr_col,"LFT"=lft_col)) +
  scale_y_continuous(limits=c(0,1.05),breaks=seq(0,1,by=0.2),expand=c(0,0)) +
  scale_x_continuous(limits=c(0,duration_pos),breaks=seq(0,duration_pos,by=2)) +
  theme_classic() + 
  theme(legend.position=c(0.8,0.8),legend.background = element_blank()) +
  theme_overall +
  ylab("Test sensitivity") +
  xlab("Days post infection")

############################################
## 2. Run the simulation
############################################
## =====================
## Simulate an incidence curve
## =====================
## Exponentially growing/declining incidence over time
tmp <- exp(growth_rate*times)

## Scale so that only cumu_inc people are infected in the focal period
inc <- cumu_inc*tmp/sum(tmp[burnin:max(times)])
sum(inc[burnin:max(times)])

## How many infections per day?
inc_dat <- tibble(inf_time=times,per_cap_inc=inc,inc=floor(inc*N))
## Raw incidence at start and end of simulation
inc_dat %>% filter(inf_time %in% c(burnin, max(times))) %>% 
  mutate(per_cap_inc*100000)

p_inc <- ggplot(inc_dat) + 
  geom_bar(aes(x=inf_time,y=per_cap_inc),stat="identity") + 
  ylab("Incidence") +
  xlab("Time (days)") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(limits=c(burnin-0.5,max(times)+0.5),breaks=seq(burnin,max(times),by=10),
                     labels=seq(0,max(times)-burnin,by=10)) +
  theme_classic() +
  theme_overall




## =====================
## Simulate linelist
## =====================
sim_dat_linelist <- inc_dat %>% dplyr::select(-per_cap_inc) %>% uncount(inc) %>% mutate(i=1:n())
n_infected <- nrow(sim_dat_linelist)

## Subsample for speed
subsamp_frac <- 0.1
sim_dat_linelist <- sim_dat_linelist %>% 
  sample_frac(subsamp_frac) %>%  mutate(i=1:n())


infection_statuses_linelist <- expand_grid(i=1:nrow(sim_dat_linelist),t=times[times >= burnin]) %>% 
  left_join(sim_dat_linelist) %>%
  mutate(days_post_inf = t - inf_time) %>%
  left_join(max_sens) %>%
  mutate(detect_prob_pcr=ifelse(is.na(detect_prob_pcr),0,detect_prob_pcr)) %>%
  mutate(detect_prob_lft=ifelse(is.na(detect_prob_lft),0,detect_prob_lft)) %>%
  mutate(true_pos = ifelse(t >= inf_time & t < inf_time + duration_pos,1,0)) %>%
  mutate(true_pcr_pos = true_pos*rbinom(n(), size=1, prob=detect_prob_pcr),
         false_pcr_pos = (1-true_pos)*rbinom(n(), size=1, prob=1-pcr_spec),
         true_lft_pos = true_pos*rbinom(n(),size=1,prob=detect_prob_lft),
         false_lft_pos = (1-true_pos)*rbinom(n(), size=1, prob=1-lft_spec))

negative_statuses_linelist <- expand_grid(i = max(infection_statuses_linelist$i):floor(N*subsamp_frac), true_pos=0, t=times) %>% 
  mutate(false_pcr_pos=rbinom(n(),size=1,prob=1-pcr_spec),
         true_pcr_pos=0,
         true_lft_pos=0,
         false_lft_pos = (1-true_pos)*rbinom(n(), size=1, prob=1-lft_spec))

infection_statuses_linelist <- bind_rows(infection_statuses_linelist, negative_statuses_linelist)
infection_statuses_linelist <- infection_statuses_linelist %>% 
  mutate(pcr_pos = as.numeric(true_pcr_pos | false_pcr_pos),
         lft_pos=as.numeric(true_lft_pos | false_lft_pos))


## Now let's see what happens when we test every 3 days by LFT and remove positives
infection_statuses_linelist <- infection_statuses_linelist %>% 
  mutate(tested=ifelse(t %in% seq(burnin,max(times),by=test_frequency),1,0))

## If test positive by LFT, remove
first_detection_dat <- infection_statuses_linelist %>% group_by(i) %>% 
  filter(tested == 1 & lft_pos == 1) %>% filter(t == min(t)) %>% 
  mutate(first_detection=t) %>% dplyr::select(i, first_detection)

infection_statuses_surveillance_linelist <- infection_statuses_linelist %>% left_join(first_detection_dat)

## Find proportion positive
lft_pos_surveillance <- infection_statuses_surveillance_linelist %>% 
  filter(t <= first_detection | is.na(first_detection), tested==1) %>% 
  group_by(t) %>% 
  summarize(prop_pos=sum(lft_pos)/n())

observed_false_positives <- infection_statuses_surveillance_linelist %>% 
  filter(t <= first_detection | is.na(first_detection), tested==1) %>% 
  group_by(t) %>% 
  summarize(prop_pos=sum(true_lft_pos)/sum(lft_pos))

observed_false_negatives <- infection_statuses_surveillance_linelist %>% 
  filter(t <= first_detection | is.na(first_detection), tested==1) %>% 
  group_by(t) %>% 
  summarize(prop_neg=sum(true_pos == 0 & lft_pos == 0)/sum(lft_pos==0))


p_prev_linelist <- ggplot(lft_pos_surveillance) + 
  geom_point(data=lft_pos_surveillance %>% filter(t >= burnin),aes(x=t,y=prop_pos)) +
  geom_line(data=lft_pos_surveillance %>% filter(t >= burnin),aes(x=t,y=prop_pos),linetype="dotted") +
  ylab("Prevalence") + 
  xlab("Time (days)") + 
  theme_classic() + 
  scale_y_continuous(expand=c(0,0)) +
  theme_overall

p_ppv_linelist <- ggplot() + 
  geom_point(data=observed_false_positives %>% filter(t >= burnin),aes(x=t,y=prop_pos)) +
  geom_line(data=observed_false_positives %>% filter(t >= burnin),aes(x=t,y=prop_pos),linetype="dotted") +
  ylab("Positive predictive value") + 
  xlab("Time (days)") + 
  theme_classic() + 
  scale_y_continuous(expand=c(0,0),limits=c(0,1))+
  theme_overall

indiv_key <-  infection_statuses_linelist %>% 
  filter(!is.na(inf_time)) %>%
  filter(inf_time >= burnin) %>%
  dplyr::select(inf_time, i) %>%
  distinct() %>%
  arrange(inf_time) %>%
  mutate(id=1:n())

infection_statuses_pcr <- infection_statuses_linelist %>% 
  filter(t >= burnin) %>%
  mutate(Status=NA,
         Status=ifelse(true_pos == 0 & (true_pcr_pos == 0 & false_pcr_pos==0),"True negative", Status),
         Status=ifelse(true_pos == 0 & (true_pcr_pos == 1 | false_pcr_pos==1),"False positive", Status),
         Status=ifelse(true_pos == 1 & (true_pcr_pos == 0 & false_pcr_pos == 0),"False negative",Status),
         Status=ifelse(t >= inf_time + duration_pos, "Recovered", Status),
         Status=ifelse((true_pcr_pos == 1 | false_pcr_pos == 1) & true_pos == 1, "True positive", Status),
         Status=ifelse(false_pcr_pos == 1 & true_pos == 0, "False positive", Status))
infection_statuses_pcr$Status <- factor(infection_statuses_pcr$Status, levels=status_levels)
positive_infection_statuses_pcr <-  infection_statuses_pcr %>%
  filter(!is.na(inf_time)) %>%
  left_join(indiv_key)


p1 <- positive_infection_statuses_pcr %>% 
  filter(inf_time >= burnin) %>%
  ggplot() + 
  geom_tile(aes(x=t,y=id,fill=Status)) +
  theme_classic() +
  ggtitle("Daily PCR surveillance") +
  scale_fill_manual(values=linelist_color_key) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("Individual") +
  scale_x_continuous(limits=c(burnin-0.5,max(times)+0.5),expand=c(0,0),
                     breaks=seq(burnin,max(times),by=10),labels=seq(0,max(times)-burnin,by=10)) +
  xlab("Time (days)")+
  theme_overall

## Tally
pcr_surveillance_tally <- infection_statuses_pcr %>% 
  filter(t == t_choose) %>%
  mutate(Status=ifelse(is.na(inf_time) & pcr_pos == 0,"True negative",Status)) %>%
  group_by(Status) %>%
  tally() %>% pivot_wider(names_from=Status,values_from=n) %>%
  mutate(`True negative/\nRecovered`=`True negative` + `Recovered`) %>%
  mutate(i=1) %>%
  pivot_longer(-i) %>%
  filter(!(name %in% c("True negative","Recovered"))) %>%
  rename(Status=name) %>%
  mutate(value=value/subsamp_frac)
pcr_surveillance_tally$Status <- factor(pcr_surveillance_tally$Status,
                                        levels=status_levels) 
p_pcr_left <- ggplot(pcr_surveillance_tally %>% filter(Status != "True negative/\nRecovered")) +
  geom_bar(aes(x=Status,y=value,fill=Status),stat="identity") +
  scale_fill_manual(values=linelist_color_key) +
  ylab("Count") +
  theme_classic() +
  scale_y_continuous(expand=c(0,0),limits=c(0,12000)) +
  theme(legend.position="none")+
  theme_overall+
  theme(axis.text.x=element_text(size=5,angle=45,hjust=1),
        axis.text.y=element_text(size=5))

p_pcr_right <-  ggplot(pcr_surveillance_tally %>% filter(Status == "True negative/\nRecovered")) +
  geom_bar(aes(x=Status,y=value,fill=Status),stat="identity") +
  scale_fill_manual(values=linelist_color_key) +
  ylab("") +
  theme_classic() +
  scale_y_continuous(expand=c(0,0),limits=c(0,N)) +
  theme(legend.position="none")+
  theme_overall+
  theme(axis.text.x=element_text(size=5,angle=45,hjust=1),
        axis.text.y=element_text(size=5))

p_pcr_tallies <- p_pcr_left + p_pcr_right + plot_layout(widths=c(3,1))
  
infection_statuses_lft <- infection_statuses_linelist %>%
  filter(t >= burnin) %>%
  mutate(Status=NA,
         Status=ifelse(true_pos == 0 & (true_lft_pos == 0 & false_lft_pos==0),"True negative", Status),
         Status=ifelse(true_pos == 0 & (true_lft_pos == 1 | false_pcr_pos==1),"False positive", Status),
         Status=ifelse(true_pos == 1 & (true_lft_pos == 0 & false_lft_pos == 0),"False negative",Status),
         Status=ifelse(t >= inf_time + duration_pos, "Recovered", Status),
         Status=ifelse((true_lft_pos == 1 | false_lft_pos == 1) & true_pos == 1, "True positive", Status),
         Status=ifelse(false_lft_pos == 1 & true_pos == 0, "False positive", Status))


positive_infection_statuses_lft <-  infection_statuses_lft %>%
  filter(!is.na(inf_time)) %>%
  left_join(indiv_key)


## Tally
lft_surveillance_tally <- infection_statuses_lft %>% 
  filter(t == t_choose) %>%
  mutate(Status=ifelse(is.na(inf_time) & lft_pos == 0,"True negative",Status)) %>%
  group_by(Status) %>%
  tally() %>% pivot_wider(names_from=Status,values_from=n) %>%
  mutate(`True negative/\nRecovered`=`True negative` + `Recovered`) %>%
  mutate(i=1) %>%
  pivot_longer(-i) %>%
  filter(!(name %in% c("True negative","Recovered"))) %>%
  rename(Status=name) %>%
  mutate(value=value/subsamp_frac)
lft_surveillance_tally$Status <- factor(lft_surveillance_tally$Status,
                                        levels=status_levels) 
p_lft_left <- ggplot(lft_surveillance_tally %>% filter(Status != "True negative/\nRecovered")) +
  geom_bar(aes(x=Status,y=value,fill=Status),stat="identity") +
  scale_fill_manual(values=linelist_color_key) +
  ylab("Count") +
  theme_classic() +
  scale_y_continuous(expand=c(0,0),limits=c(0,12000)) +
  theme(legend.position="none")+
  theme_overall+
  theme(axis.text.x=element_text(size=5,angle=45,hjust=1),
                                 axis.text.y=element_text(size=5))

p_lft_right <-  ggplot(lft_surveillance_tally %>% filter(Status == "True negative/\nRecovered")) +
  geom_bar(aes(x=Status,y=value,fill=Status),stat="identity") +
  scale_fill_manual(values=linelist_color_key) +
  ylab("") +
  theme_classic() +
  scale_y_continuous(expand=c(0,0),limits=c(0,N)) +
  theme(legend.position="none")+
  theme_overall+
  theme(axis.text.x=element_text(size=5,angle=45,hjust=1),
        axis.text.y=element_text(size=5))

p_lft_tallies <- p_lft_left + p_lft_right + plot_layout(widths=c(3,1))

p2 <- positive_infection_statuses_lft %>% 
  filter(inf_time >= burnin) %>%
  ggplot() + geom_tile(aes(x=t,y=id,fill=Status)) +
  theme_classic() +
  ggtitle("Daily LFT surveillance") +
  scale_fill_manual(values=linelist_color_key) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("Individual") +
  scale_x_continuous(limits=c(burnin-0.5,max(times)+0.5),expand=c(0,0),
                     breaks=seq(burnin,max(times),by=10),labels=seq(0,max(times)-burnin,by=10)) +
  xlab("Time (days)")+
  theme_overall


indiv_key_surveillance <-  infection_statuses_surveillance_linelist %>% 
  filter(!is.na(inf_time)) %>%
  filter(inf_time >= burnin) %>%
  dplyr::select(inf_time, i) %>%
  distinct() %>%
  arrange(inf_time) %>%
  mutate(id=1:n())


infection_statuses_lft_surveillance <- infection_statuses_surveillance_linelist %>%
  filter(t>= burnin) %>%
  mutate(Status=NA,
         Status=ifelse(tested ==0, "Not tested",Status),
         Status=ifelse(tested==1 & true_pos == 0 & (true_lft_pos == 0 & false_lft_pos==0),"True negative", Status),
         Status=ifelse(tested==1 & true_pos == 0 & (true_lft_pos == 1 | false_lft_pos==1),"False positive", Status),
         Status=ifelse(tested == 1 & (is.na(first_detection) | t <= first_detection) & true_pos == 1 & (true_lft_pos == 0 & false_lft_pos == 0),"False negative",Status),
         Status=ifelse(t >= first_detection + isolation_period & !is.na(first_detection), "Post isolation", Status),
         Status=ifelse(t >= inf_time + duration_pos & is.na(first_detection), "Not detected", Status),
         Status=ifelse(tested==1 & t <= first_detection & (true_lft_pos == 1 | false_lft_pos == 1) & true_pos == 1, "True positive", Status),
         Status=ifelse(tested==1 & false_lft_pos == 1 & true_pos == 0, "False positive", Status),
         Status=ifelse(!is.na(first_detection) & t > first_detection & true_pos == 1 & t <= first_detection + isolation_period,"Correctly\n isolating",Status),
         Status=ifelse(!is.na(first_detection) & t > first_detection & true_pos == 0 & t <= first_detection + isolation_period,"Incorrectly\n isolating",Status)
  )

positive_infection_statuses_lft_surveillance <-  infection_statuses_lft_surveillance %>%
  filter(!is.na(inf_time)) %>%
  left_join(indiv_key_surveillance)


## Tally
lft_repeat_tally <- infection_statuses_lft_surveillance %>% 
  filter(t == t_choose) %>%
  mutate(Status=ifelse(is.na(Status) & tested==1 & true_pos == 0 & lft_pos == 0,"True negative",Status)) %>%
  group_by(Status) %>%
  tally() %>% pivot_wider(names_from=Status,values_from=n) %>%
  mutate(`True negative/\nPost infection`=`True negative` + `Post isolation`) %>%
  mutate(i=1) %>%
  pivot_longer(-i) %>%
  filter(!(name %in% c("True negative","Post isolation","Not detected"))) %>%
  rename(Status=name) %>%
  mutate(value=value/subsamp_frac)
  
lft_repeat_tally$Status <- factor(lft_repeat_tally$Status,
                                        levels=status_levels) 
p_lft_repeat_left <- ggplot(lft_repeat_tally %>% filter(Status != "True negative/\nPost infection")) +
  geom_bar(aes(x=Status,y=value,fill=Status),stat="identity") +
  scale_fill_manual(values=linelist_color_key) +
  ylab("Count") +
  theme_classic() +
  scale_y_continuous(expand=c(0,0),limits=c(0,12000)) +
  theme(legend.position="none")+
  theme_overall +
  theme(axis.text.x=element_text(size=5,angle=45,hjust=1),
        axis.text.y=element_text(size=5))
p_lft_repeat_right <-  ggplot(lft_repeat_tally %>% filter(Status == "True negative/\nPost infection")) +
  geom_bar(aes(x=Status,y=value,fill=Status),stat="identity") +
  scale_fill_manual(values=linelist_color_key) +
  ylab("") +
  theme_classic() +
  scale_y_continuous(expand=c(0,0),limits=c(0,N)) +
  theme(legend.position="none")+
  theme_overall +
  theme(axis.text.x=element_text(size=5,angle=45,hjust=1),
        axis.text.y=element_text(size=5))

p_lft_repeat_tallies <- p_lft_repeat_left + p_lft_repeat_right + plot_layout(widths=c(3,1))


p3 <- positive_infection_statuses_lft_surveillance %>% 
  filter(inf_time >= burnin) %>%
  ggplot() + 
  geom_tile(aes(x=t,y=id,fill=Status)) +
  ggtitle(paste0(test_frequency,"-daily LFT test and isolate")) +
  theme_classic() +
  scale_fill_manual(values=linelist_color_key)+
  scale_y_continuous(expand=c(0,0)) +
  ylab("Individual") +
  scale_x_continuous(limits=c(burnin-0.5,max(times)+0.5),breaks=seq(burnin,max(times),by=10),
                     labels=seq(0,max(times)-burnin,by=10),expand=c(0,0)) +
  xlab("Time (days)")+
  theme_overall

p3_alt <- positive_infection_statuses_lft_surveillance %>% 
  filter(inf_time >= burnin) %>%
  filter(tested==1) %>%
  ggplot() + 
  geom_tile(aes(x=t,y=id,fill=Status)) +
  ggtitle(paste0(test_frequency,"-daily LFT test and isolate")) +
  theme_classic() +
  scale_fill_manual(values=linelist_color_key)+
  scale_y_continuous(expand=c(0,0)) +
  ylab("Individual") +
  scale_x_continuous(#limits=c(burnin,max(times)+0.5),
    breaks=seq(burnin,max(times),by=test_frequency*2)-1.5,
    labels=seq(0,max(times)-burnin,by=test_frequency*2),
                     expand=c(0,0)) +
  xlab("Time (days)")+
  theme_overall



## =====================
## Simulate cohort
## =====================
n_infected <- sum(inc_dat$inc)
## Get infection status at all times during the simulation for
## those who are INFECTED AT SOME POINT
## This tibble essentially gives a cohort of individuals infected on each inf_time, and the proportion
## of those cohorts that are true positives and test positive on each day of the simulation
cohort_dat <- expand_grid(inc_dat,t=times[times >= burnin]) %>% 
  dplyr::select(-per_cap_inc) %>%
  mutate(days_post_inf = t - inf_time) %>%
  rename(n=inc) %>%
  mutate(true_pos = ifelse(t >= inf_time & t < inf_time + duration_pos,1,0)) %>%
  ## Merge with test sensitivity to simulate test outcomes
  left_join(max_sens) %>%
  mutate(detect_prob_pcr=ifelse(is.na(detect_prob_pcr),0,detect_prob_pcr)) %>%
  mutate(detect_prob_lft=ifelse(is.na(detect_prob_lft),0,detect_prob_lft))
  
infection_statuses <- cohort_dat %>%
  ## If it is after infection time and before clearance time, you are a true positive
  mutate(true_pcr_pos = true_pos*rbinom(n(), size=n, prob=detect_prob_pcr), ## True positives by PCR
         false_pcr_pos = (1-true_pos)*rbinom(n(), size=n, prob=1-pcr_spec), ## False positive by PCR
         true_lft_pos = true_pos*rbinom(n(),size=n,prob=detect_prob_lft), ## True positives by LFT
         false_lft_pos = (1-true_pos)*rbinom(n(), size=n, prob=1-lft_spec)) ## False positives by LFT
positive_infection_statuses <- infection_statuses

## Enumerate out results for the negatives ie. THOSE WHO ARE NEVER INFECTED
cohort_dat_negative <- tibble(inf_time=-1,days_post_inf=-1,n=N-n_infected, true_pos=0, t=times,detect_prob_pcr=0,detect_prob_lft=0) 
negative_statuses <- cohort_dat_negative %>% 
  filter(t >= burnin) %>%
  mutate(false_pcr_pos=rbinom(n(),size=n,prob=1-pcr_spec), ## False positives by PCR
         true_pcr_pos=0, ## Always true negative PCR
         true_lft_pos=0, ## Always true negative LFT
         false_lft_pos = (1-true_pos)*rbinom(n(), size=n, prob=1-lft_spec)) ## False positives by PCR

## Combine infections and not infecteds and get overall PCR and LFT results on each day of the simulation
infection_statuses <- bind_rows(positive_infection_statuses, negative_statuses)
infection_statuses <- infection_statuses %>% 
  mutate(pcr_pos = true_pcr_pos + false_pcr_pos,
         lft_pos= true_lft_pos + false_lft_pos)

## =====================
## Simulate surveillance strategy
## =====================
## Now let's see what happens when we test every 3 days by LFT and remove positives
cohort_dat_all <- bind_rows(cohort_dat, cohort_dat_negative)

infection_statuses_surveillance <- cohort_dat_all %>% 
  mutate(tested=ifelse(t %in% seq(burnin,max(times),by=test_frequency),1,0)) %>%
  filter(tested == 1) %>%
  mutate(
    #true_pcr_pos=0,false_pcr_pos=0,pcr_pos=0,
    true_lft_pos=0,false_lft_pos=0,lft_pos=0) %>%
  mutate(#pcr_n=n,
         lft_n=n)

## For each cohort on each day, remove positives from this cohort for the remainder of time
## Same for negatives - remove false positives for remainder of time
cohorts <- unique(infection_statuses_surveillance$inf_time)
times <- unique(infection_statuses_surveillance$t)

cohort_dat_all_new <- NULL
for(index in seq_along(cohorts)){
  cohort <- cohorts[index]
  tmp_cohort <- infection_statuses_surveillance %>% filter(inf_time == cohort) %>% arrange(t)
  #tmp_cohort$true_pcr_pos[1] <- tmp_cohort$true_pos[1]*rbinom(1, size=tmp_cohort$pcr_n[1], prob=tmp_cohort$detect_prob_pcr[1])
  #tmp_cohort$false_pcr_pos[1] <- (1-tmp_cohort$true_pos[1])*rbinom(1, size=tmp_cohort$pcr_n[1], prob=1-pcr_spec)
  
  ## Simulate true positives from tested infected individuals
  tmp_cohort$true_lft_pos[1] <- tmp_cohort$true_pos[1]*rbinom(1, size=tmp_cohort$lft_n[1], prob=tmp_cohort$detect_prob_lft[1])
  ## Simulate false positives from tested not infected individuals
  tmp_cohort$false_lft_pos[1] <- (1-tmp_cohort$true_pos[1])*rbinom(1, size=tmp_cohort$lft_n[1], prob=1-lft_spec)
  
  #tmp_cohort$pcr_pos[1] <-  tmp_cohort$true_pcr_pos[1] + tmp_cohort$false_pcr_pos[1]
  tmp_cohort$lft_pos[1] <-  tmp_cohort$true_lft_pos[1] + tmp_cohort$false_lft_pos[1]
  
  for(t in 2:nrow(tmp_cohort)){
    #tmp_cohort$pcr_n[t] <- tmp_cohort$pcr_n[t-1] - tmp_cohort$pcr_pos[t-1]
    tmp_cohort$lft_n[t] <- tmp_cohort$lft_n[t-1] - tmp_cohort$lft_pos[t-1]
    
    #tmp_cohort$true_pcr_pos[t] <- tmp_cohort$true_pos[t]*rbinom(1, size=tmp_cohort$pcr_n[t], prob=tmp_cohort$detect_prob_pcr[t])
    #tmp_cohort$false_pcr_pos[t] <- (1-tmp_cohort$true_pos[t])*rbinom(1, size=tmp_cohort$pcr_n[t], prob=1-pcr_spec)
    tmp_cohort$true_lft_pos[t] <- tmp_cohort$true_pos[t]*rbinom(1, size=tmp_cohort$lft_n[t], prob=tmp_cohort$detect_prob_lft[t])
    tmp_cohort$false_lft_pos[t] <- (1-tmp_cohort$true_pos[t])*rbinom(1, size=tmp_cohort$lft_n[t], prob=1-lft_spec)
    
    #tmp_cohort$pcr_pos[t] <-  tmp_cohort$true_pcr_pos[t] + tmp_cohort$false_pcr_pos[t]
    tmp_cohort$lft_pos[t] <-  tmp_cohort$true_lft_pos[t] + tmp_cohort$false_lft_pos[t]
  }
  cohort_dat_all_new[[index]] <- tmp_cohort
}
cohort_dat_all_new <- do.call("bind_rows", cohort_dat_all_new)

## Find proportion positive on each day of surveillance
lft_surveillance <- cohort_dat_all_new %>% 
  mutate(true_pos=true_pos*lft_n) %>%
  group_by(t) %>% 
  summarize(
    tests=sum(lft_n),
    infected=sum(true_pos),
    not_infected=tests-infected, ## Everyone else is in isolation
    pos=sum(lft_pos),
    neg=tests-pos,
    true_pos=sum(true_lft_pos),
    false_pos=pos-true_pos,
    true_neg=not_infected-false_pos,
    false_neg=neg-true_neg) 

## Summary statistics
lft_surveillance_summary <- lft_surveillance %>%
  mutate(`Repeated screening prevalence`=pos/tests,
         `Repeated screening sensitivity`=true_pos/infected,
         `Repeated screening specificity`=true_neg/(true_neg+false_pos),
         `PPV repeated screening`=true_pos/pos,
         `NPV repeated screening`=true_neg/neg)


############################################
## 3. Calculate statistics
############################################
## =====================
## Calculate all statistics from the screening strategy
## =====================
prev_dat <- infection_statuses %>% 
  group_by(t) %>% 
  mutate(true_pos=true_pos*n) %>%
  summarize(tests=sum(n),
            infected=sum(true_pos),
            not_infected=N-infected,
            
            ## PCR results
            pcr_pos=sum(pcr_pos),
            pcr_neg=tests-pcr_pos,
            pcr_true_pos=sum(true_pcr_pos),
            pcr_false_pos=pcr_pos-pcr_true_pos,
            pcr_true_neg=not_infected-pcr_false_pos,
            pcr_false_neg=pcr_neg-pcr_true_neg,
            
            ## LFT results
            lft_pos=sum(lft_pos),
            lft_neg=tests-lft_pos,
            lft_true_pos=sum(true_lft_pos),
            lft_false_pos=lft_pos-lft_true_pos,
            lft_true_neg=not_infected-lft_false_pos,
            lft_false_neg=lft_neg-lft_true_neg) %>%
  mutate(pcr_check=pcr_true_neg+pcr_false_neg+pcr_true_pos+pcr_false_pos,
         lft_check=lft_true_neg+lft_false_neg+lft_true_pos+lft_false_pos)

## Get statistics
prev_dat_summary <- prev_dat %>%
  mutate(`Truth`=infected/N,
         `qPCR`=pcr_true_pos/tests, 
         `LFT`=lft_pos/tests, 
         `qPCR sensitivity`=pcr_true_pos/infected,
         `LFT sensitivity`=lft_true_pos/infected,
         `qPCR specificity`=pcr_true_neg/(pcr_true_neg+pcr_false_pos),
         `LFT specificity`=lft_true_neg/(lft_true_neg+lft_false_pos),
         `PPV qPCR`=pcr_true_pos/pcr_pos,
         `PPV LFT`=lft_true_pos/lft_pos,
         `NPV qPCR`=pcr_true_neg/pcr_neg,
         `NPV LFT`=lft_true_neg/lft_neg)


## =====================
## Positive predictive value
## =====================
PPV_dat <- prev_dat_summary %>%
  dplyr::select(-c(qPCR,LFT)) %>%
  rename(`qPCR`=`PPV qPCR`,`LFT`=`PPV LFT`) %>%
  dplyr::select(t, qPCR, LFT) %>%
   pivot_longer(-t) %>%
  rename(`Positive predictive value`=name)
PPV_dat %>% tail()
PPV_dat <- bind_rows(PPV_dat, 
                     lft_surveillance_summary %>% 
                       dplyr::select(t, `PPV repeated screening`) %>% 
                       rename(value=`PPV repeated screening`) %>% 
                       mutate(`Positive predictive value`="Regular screening"))

p_ppv <- ggplot(PPV_dat) + 
  geom_line(aes(x=t,y=value,col=`Positive predictive value`,linetype=`Positive predictive value`),size=1) + 
  geom_point(data=lft_surveillance_summary %>% rename(`Regular screening`=`PPV repeated screening`),
             aes(x=t,y=`Regular screening`),col="black",size=2) +
  scale_color_manual(values=c("LFT"=lft_col,"qPCR"=pcr_col,"Regular screening"="black")) +
  scale_linetype_manual(values=c("solid","solid","dotted")) +
  ylab("Positive predictive value") + 
  xlab("Time (days)") + 
  theme_classic() + 
  scale_y_continuous(expand=c(0,0),limits=c(0,1.01),breaks=seq(0,1,by=0.2)) +
  scale_x_continuous(limits=c(burnin-0.5,max(times)+0.5),breaks=seq(burnin,max(times),by=10),
                     labels=seq(0,max(times)-burnin,by=10),expand=c(0,0)) +
  theme(legend.position="none",legend.title=element_blank())+
  theme_overall

## =====================
## Negative predictive value
## =====================
NPV_dat <- prev_dat_summary %>%
  dplyr::select(-c(qPCR,LFT)) %>%
  rename(`qPCR`=`NPV qPCR`,`LFT`=`NPV LFT`) %>%
  dplyr::select(t, qPCR, LFT) %>%
  pivot_longer(-t) %>%
  rename(`Negative predictive value`=name)
NPV_dat %>% tail()
NPV_dat <- bind_rows(NPV_dat, 
                     lft_surveillance_summary %>% 
                       dplyr::select(t, `NPV repeated screening`) %>% 
                       rename(value=`NPV repeated screening`) %>% 
                       mutate(`Negative predictive value`="Regular screening"))
p_npv <- ggplot(NPV_dat) + 
  geom_line(aes(x=t,y=value,col=`Negative predictive value`,linetype=`Negative predictive value`),size=1) + 
  geom_point(data=lft_surveillance_summary %>% rename(`Regular screening`=`NPV repeated screening`),
             aes(x=t,y=`Regular screening`),col="black",size=2) +
  scale_color_manual(values=c("LFT"=lft_col,"qPCR"=pcr_col,"Regular screening"="black")) +
  scale_linetype_manual(values=c("solid","solid","dotted")) +
  ylab("Negative predictive value") + 
  xlab("Time (days)") + 
  theme_classic() + 
  scale_x_continuous(limits=c(burnin-0.5,max(times)+0.5),breaks=seq(burnin,max(times),by=10),
                     labels=seq(0,max(times)-burnin,by=10),expand=c(0,0)) +
  theme(legend.position="none",legend.title=element_blank())+
  theme_overall

## =====================
## Prevalence
## =====================
## Get true prevalence over time
prev_dat <- prev_dat_summary %>% 
 dplyr::select(t, `Truth`,`qPCR`,`LFT`) %>%
  pivot_longer(-t) %>% rename(Prevalence=name)

## Prevalence in screening population
prev_dat <- bind_rows(prev_dat, lft_surveillance_summary %>% 
                        dplyr::select(t, `Repeated screening prevalence`) %>%
                        rename(value=`Repeated screening prevalence`) %>%
                        mutate(Prevalence="Regular screening"))

prev_dat %>% filter(t == max(times), Prevalence=="Truth") %>% mutate(value*100)

inc_yscale <- 20

p_prev <- ggplot(prev_dat) + 
  geom_bar(data=inc_dat %>% filter(inf_time >= burnin),aes(x=inf_time,y=per_cap_inc*inc_yscale,fill="Incidence"),stat="identity",alpha=0.5,size=0.1) + 
  geom_line(aes(x=t,y=value,col=Prevalence,linetype=Prevalence),size=1) + 
  geom_point(data=prev_dat %>% filter(Prevalence=="Regular screening"),aes(x=t,y=value),size=2) +
  scale_x_continuous(limits=c(burnin-0.5,max(times)+0.5),expand=c(0,0),
                     breaks=seq(burnin,max(times),by=10),labels=seq(0,max(times)-burnin,by=10)) +
  scale_fill_manual(values=c("Incidence"="grey70")) +
  scale_color_manual(values=c("LFT"="#E37222FF","qPCR"="#2A6EBBFF",
                              "Truth"="#CD202CFF","Regular screening"="black")) +
  scale_linetype_manual(values=c("LFT"="solid","qPCR"="solid",
                                 "Truth"="solid","Regular screening"="dotted"))+
  ylab("Prevalence") + 
  xlab("Time (days)") + 
  theme_classic() + 
  theme(legend.position="top",legend.title=element_blank()) +
  scale_y_continuous(expand=c(0,0),#limits=c(0,0.02),breaks=seq(0,0.02,by=0.005),
                     sec.axis=sec_axis(trans=~.*1/inc_yscale,name="Incidence"))+
  theme_overall

############################################
## 4. Save plots
############################################
## Check all prevalence stats
p_prev_check <- prev_dat_summary %>% pivot_longer(-t) %>% ggplot() + geom_line(aes(x=t,y=value)) + facet_wrap(~name,scales="free_y")
## Check all regular screening stats
p_lft_check <- lft_surveillance_summary %>% pivot_longer(-t) %>% ggplot() + geom_line(aes(x=t,y=value)) + facet_wrap(~name,scales="free_y")

p_indiv_strategies <- (p1+theme(legend.position="none")+geom_vline(xintercept=t_choose,linetype="dashed",col="red")) +
  (p2+geom_vline(xintercept=t_choose,linetype="dashed",col="red")+
     theme(legend.position="none",axis.text.y=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank()))+
  (p3_alt+geom_vline(xintercept=t_choose,linetype="dashed",col="red")+
     theme(axis.text.y=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank())) +
  (p_pcr_left + theme(axis.title.y=element_blank(),axis.title.x=element_blank()))+
  (p_lft_left +  theme(axis.title.y=element_blank(),axis.title.x=element_blank()))+
  (p_lft_repeat_left + theme(axis.title.y=element_blank(),axis.title.x=element_blank()))+
  plot_layout(guides="collect",nrow=2, heights=c(2.5,1))


p_main2 <- (p_prev+theme(axis.title.x=element_blank()))/
  (p_ppv+theme(axis.title.x=element_blank()))/
  p_npv / p_indiv_strategies

p_main <- (p_prev+theme(axis.title.x=element_blank()))/
  (p_ppv+theme(axis.title.x=element_blank()))/
  p_npv

  
ggsave("main_p1.png",p_main, height=8,width=6,units="in",dpi=300)
ggsave("main_p2.png",p_main2, height=8,width=6,units="in",dpi=300)
ggsave("indiv_strategies.png",p_indiv_strategies, height=6,width=8,units="in",dpi=300)
