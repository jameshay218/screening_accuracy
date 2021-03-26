#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(patchwork)
library(shinyBS)

## Plot settings
lft_col <- "#E37222FF"
pcr_col <-"#2A6EBBFF"


linelist_color_key <- c("False negative"="#7D5CC6", 
                        "Not detected"="#CD202CFF", 
                        "True positive"="#69BE28FF", 
                        "Post isolation"="#7FA8D6FF", 
                        "Recovered"="#7FA8D6FF",
                        "False positive"="#F0AB00FF",
                        "Correctly isolating"="#00B2A9", 
                        "Correctly\n isolating"="#00B2A9", 
                        "Correctly isolating/positive"="#00B2A9", 
                        "Incorrectly\n isolating"="black", 
                        "Incorrectly isolating"="black", 
                        "Incorrectly isolating/negative"="black", 
                        "Incorrectly\n not isolating"="orange",
                        "Incorrectly not\n isolating/positive"="orange",
                        "True negative"="#ACADAEFF",
                        "True negative/\nRecovered"="#ACADAEFF",
                        "True negative/\nPost infection"="#ACADAEFF",
                        "Not tested"="white")

status_levels <- c("True positive","False positive","True negative","False negative","Recovered",
                   "Post isolation","Correctly\n isolating","Incorrectly\n isolating","Incorrectly\n not isolating",
                   "Not detected",
                   "True negative/\nPost infection","True negative/\nRecovered",
                   "Not tested")

theme_overall <- overall_theme <- theme(axis.text.x=element_text(size=12),
                       axis.title.x=element_text(size=14),
                       axis.text.y=element_text(size=12),
                       axis.title.y=element_text(size=16),
                       legend.text=element_text(size=14),
                       legend.title=element_text(size=14),
                       plot.title=element_text(size=14,face="bold")
                       )

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
    output$test_characteristics1 <- renderText({storage_for_plots()$test_text1})
    output$test_characteristics2 <- renderText({storage_for_plots()$test_text2})
    
    
    storage_for_plots <- reactive({
        parameters <- reactiveValuesToList(input)
        growth_rate <- parameters$growth_rate
        ## Simulation duration
        ## Cumulative incidence in the focal period (ie. between burnin and max(times))
        cumu_inc <- parameters$cumu_inc
        ## Population size
        N <- parameters$N
        ## For how long are people true infections?
        duration_pos <- parameters$duration_pos
        
        ## Burning for simulation
        burnin <- duration_pos + 25
        sim_duration <- parameters$sim_duration
        times <- 1:(burnin + sim_duration)
        
        ## =====================
        ## PCR surveillance
        ## =====================
        ## PCR sensitivity
        pcr_gamma_mode <- parameters$pcr_gamma_mode
        pcr_gamma_sd <- parameters$pcr_gamma_sd
        pcr_gamma_height <- parameters$pcr_gamma_height
        pcr_gamma_max_sens <- parameters$pcr_gamma_max_sens
        
        pcr_gamma_rate = ( pcr_gamma_mode + sqrt( pcr_gamma_mode^2 + 4*pcr_gamma_sd^2 ) ) / ( 2 * pcr_gamma_sd^2 )
        pcr_gamma_shape = 1 + pcr_gamma_mode * pcr_gamma_rate
        detect_prob_pcr <- dgamma(0:duration_pos,shape = pcr_gamma_shape,rate=pcr_gamma_rate)
        detect_prob_pcr <- detect_prob_pcr*(pcr_gamma_height/max(detect_prob_pcr))
        detect_prob_pcr <- pmin(detect_prob_pcr, pcr_gamma_max_sens)
        
        
        detect_prob_pcr_resolved <- dgamma(seq(0,duration_pos,by=0.01),shape = pcr_gamma_shape,rate=pcr_gamma_rate)
        detect_prob_pcr_resolved <- detect_prob_pcr_resolved*(pcr_gamma_height/max(detect_prob_pcr_resolved))
        detect_prob_pcr_resolved <- pmin(detect_prob_pcr_resolved, pcr_gamma_max_sens)
        
        ## =====================
        ## LFT strategy
        ## =====================
        ## LFT sensitivity
        lft_gamma_mode <- pcr_gamma_mode ## Assume peaks at same time
        lft_gamma_sd <- parameters$lft_gamma_sd
        lft_gamma_height <- parameters$lft_gamma_height
        lft_gamma_max_sens <- parameters$lft_gamma_max_sens
        
        lft_gamma_rate = ( lft_gamma_mode + sqrt( lft_gamma_mode^2 + 4*lft_gamma_sd^2 ) ) / ( 2 * lft_gamma_sd^2 )
        lft_gamma_shape = 1 + lft_gamma_mode * lft_gamma_rate
        detect_prob_lft <- dgamma(0:duration_pos,shape = lft_gamma_shape,rate=lft_gamma_rate)
        detect_prob_lft <- detect_prob_lft*(lft_gamma_height/max(detect_prob_lft))
        detect_prob_lft <- pmin(detect_prob_lft, lft_gamma_max_sens)
        
        detect_prob_lft_resolved <- dgamma(seq(0,duration_pos,by=0.01),shape = lft_gamma_shape,rate=lft_gamma_rate)
        detect_prob_lft_resolved <- detect_prob_lft_resolved*(lft_gamma_height/max(detect_prob_lft_resolved))
        detect_prob_lft_resolved <- pmin(detect_prob_lft_resolved, lft_gamma_max_sens)
        
        ## =====================
        ## Specificity and test frequency
        ## =====================
        ## PCR spec
        pcr_spec <- parameters$pcr_spec
        ## LFT spec
        lft_spec <- parameters$lft_spec
        
        ## LFT screening approach
        test_frequency <- parameters$test_frequency
        isolation_period <- parameters$isolation_period
        
        ## =====================
        ## Time-varying maximum sensitivity over infection ie. within-host model
        ## =====================
        max_sens <- tibble(days_post_inf=0:duration_pos, 
                           detect_prob_pcr=detect_prob_pcr,
                           detect_prob_lft=detect_prob_lft)
        
        
        max_sens_plot <- tibble(days_post_inf=seq(0,duration_pos,by=0.01), 
                                detect_prob_pcr=detect_prob_pcr_resolved,
                                detect_prob_lft=detect_prob_lft_resolved)
        
        
        text1 <- paste0("Overall LFT mean sensitivity: ", 
                        signif(mean(detect_prob_lft)*100,3),
                        "%; and qPCR mean sensitivity: ", 
                        signif(mean(detect_prob_pcr)*100,3),"%")
        text2 <-paste0("Day 2-12 LFT mean sensitivity: ", 
                        signif(mean(detect_prob_lft[2:max(12,length(detect_prob_lft))])*100,3),
                        "%; and qPCR mean sensitivity: ", 
                        signif(mean(detect_prob_pcr[2:max(12,length(detect_prob_pcr))])*100,3),
                       "%"
                        )
        
        
        return(list(
            "duration_pos"=duration_pos,
            "max_sens"=max_sens,
            "max_sens_dat"=max_sens_plot,
            "test_text1"=text1,
            "test_text2"=text2
            ))
    })
    
    
    output$assumptions_plot <- renderPlot({
        tmp <- storage_for_plots()
        
        p_test_sens <- tmp$max_sens_dat %>% 
            rename(`qPCR`=detect_prob_pcr,`LFT`=detect_prob_lft) %>%
            pivot_longer(-days_post_inf) %>%
            rename(Test=name) %>%
            ggplot() +
            geom_line(aes(x=days_post_inf,y=value,col=Test),size=1) +
            geom_vline(xintercept=5.5,linetype="dashed") +
            scale_color_manual(values=c("qPCR"=pcr_col,"LFT"=lft_col)) +
            scale_y_continuous(limits=c(0,1.05),breaks=seq(0,1,by=0.2),labels=seq(0,1,by=0.2)*100,expand=c(0,0)) +
            scale_x_continuous(limits=c(0,tmp$duration_pos),breaks=seq(0,tmp$duration_pos,by=2)) +
            theme_classic() + 
            theme(legend.position=c(0.8,0.8),legend.background = element_blank()) +
            overall_theme +
            ylab("Test sensitivity (%)") +
            xlab("Days post infection")
        p_test_sens
        
    },height=400,width=700)
    output$all_plots <- renderPlot({
        parameters <- reactiveValuesToList(input)
        growth_rate <- parameters$growth_rate
        ## Simulation duration
        ## Cumulative incidence in the focal period (ie. between burnin and max(times))
        cumu_inc <- parameters$cumu_inc
        ## Population size
        N <- parameters$N
        ## For how long are people true infections?
        duration_pos <- parameters$duration_pos
        
        ## Burning for simulation
        burnin <- duration_pos*2
        sim_duration <- parameters$sim_duration
        times <- 1:(burnin + sim_duration)
        
        strategy_start <- parameters$strategy_start_time + burnin
        strategy_start <- min(strategy_start, max(times))
        ## =====================
        ## Time-varying maximum sensitivity over infection ie. within-host model
        ## =====================
        max_sens <-  storage_for_plots()$max_sens
        
        ## =====================
        ## Specificity and test frequency
        ## =====================
        ## PCR spec
        pcr_spec <- parameters$pcr_spec
        ## LFT spec
        lft_spec <- parameters$lft_spec
        
        ## LFT screening approach
        test_frequency <- parameters$test_frequency
        isolation_period <- parameters$isolation_period
        
        testing_times <- seq(strategy_start,max(times),by=test_frequency)
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
        
        ############################################
        ## 2. PLOTS
        ############################################
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
            mutate(tested=ifelse(t %in% testing_times,1,0)) %>%
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
        
        PPV_stats <- PPV_dat %>% filter(t%in% range(lft_surveillance_summary$t), `Positive predictive value` == "Regular screening") %>% pull(value)
        
        p_ppv <- ggplot(PPV_dat) + 
            geom_line(aes(x=t,y=value*100,col=`Positive predictive value`,linetype=`Positive predictive value`),size=1) + 
            geom_vline(xintercept=strategy_start,linetype="dashed",size=0.5,col="grey40") +
            geom_point(data=lft_surveillance_summary %>% rename(`Regular screening`=`PPV repeated screening`),
                       aes(x=t,y=`Regular screening`*100),col="black",size=2.5) +
            scale_color_manual(values=c("LFT"=lft_col,"qPCR"=pcr_col,"Regular screening"="black")) +
            scale_linetype_manual(values=c("solid","solid","dotted")) +
            ylab("Positive predictive value (%)") + 
            xlab("Time (days)") + 
            ggtitle(paste0("Positive predictive value of screening ranges from ", signif(PPV_stats[1]*100,3), "% to ", signif(PPV_stats[2]*100,3),"%")) +
            theme_classic() + 
            scale_x_continuous(
                limits=c(burnin-0.5,max(times)+0.5),
                breaks=seq(burnin,max(times),by=test_frequency*2),
                labels=seq(burnin,max(times),by=test_frequency*2)-burnin,
                expand=c(0,0)) +
            theme(legend.position="none",legend.title=element_blank())+
            theme_overall +
            theme(axis.title.x = element_blank())
        
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
        
        NPV_stats <- NPV_dat %>% filter(t%in% range(lft_surveillance_summary$t), `Negative predictive value` == "Regular screening") %>% pull(value)
        p_npv <- ggplot(NPV_dat) + 
            geom_line(aes(x=t,y=value*100,col=`Negative predictive value`,linetype=`Negative predictive value`),size=1) + 
            geom_point(data=lft_surveillance_summary %>% rename(`Regular screening`=`NPV repeated screening`),
                       aes(x=t,y=`Regular screening`*100),col="black",size=2.5) +
            geom_vline(xintercept=strategy_start,linetype="dashed",size=0.5,col="grey40") +
            scale_color_manual(values=c("LFT"=lft_col,"qPCR"=pcr_col,"Regular screening"="black")) +
            scale_linetype_manual(values=c("solid","solid","dotted")) +
            ylab("Negative predictive value (%)") + 
            xlab("Time (days)") + 
            ggtitle(paste0("Negative predictive value of screening ranges from ", signif(NPV_stats[1]*100,3), "% to ", signif(NPV_stats[2]*100,3),"%")) +
            theme_classic() + 
            scale_x_continuous(
                limits=c(burnin-0.5,max(times)+0.5),
                breaks=seq(burnin,max(times),by=test_frequency*2),
                labels=seq(burnin,max(times),by=test_frequency*2)-burnin,
                expand=c(0,0)) +
            theme_overall +
            theme(legend.position=c(0.75,0.4),legend.title=element_blank())
        
        ## =====================
        ## Prevalence
        ## =====================
        ## Get true prevalence over time
        prev_dat <- prev_dat_summary %>% 
            dplyr::select(t, `Truth`,`qPCR`,`LFT`) %>%
            rename(`True prevalence`=Truth,`qPCR percent positive`=`qPCR`,`LFT percent positive`=LFT) %>%
            pivot_longer(-t) %>% rename(Prevalence=name)
        
        ## Prevalence in screening population
        prev_dat <- bind_rows(prev_dat, lft_surveillance_summary %>% 
                                  dplyr::select(t, `Repeated screening prevalence`) %>%
                                  rename(value=`Repeated screening prevalence`) %>%
                                  mutate(Prevalence="Regular screening"))
        
        prev_dat %>% filter(t == max(times), Prevalence=="Truth") %>% mutate(value*100)
        
        inc_yscale <- 1#10/100
        
        prev_stats <- prev_dat %>% filter(t %in% c(burnin, max(lft_surveillance_summary$t)),Prevalence=="True prevalence") %>% mutate(value=value*100) %>% pull(value)
        
        p_prev <- ggplot(prev_dat) + 
            geom_bar(data=inc_dat %>% filter(inf_time >= burnin),aes(x=inf_time,y=(per_cap_inc*100)*inc_yscale,fill="Incidence"),
                     stat="identity",alpha=0.5,size=0.1) + 
            geom_line(aes(x=t,y=value*100,col=Prevalence,linetype=Prevalence),size=1) + 
            geom_vline(xintercept=strategy_start,linetype="dashed",size=0.5,col="grey40") +
            geom_point(data=prev_dat %>% filter(Prevalence=="Regular screening"),aes(x=t,y=value*100),size=2.5) +
            scale_x_continuous(
                limits=c(burnin-0.5,max(times)+0.5),
                breaks=seq(burnin,max(times),by=test_frequency*2),
                labels=seq(burnin,max(times),by=test_frequency*2)-burnin,
                expand=c(0,0)) +
            scale_fill_manual(name=NULL,values=c("Incidence"="grey70")) +
            scale_color_manual(name=NULL,values=c("LFT percent positive"="#E37222FF","qPCR percent positive"="#2A6EBBFF",
                                                  "True prevalence"="#CD202CFF","Regular screening"="black")) +
            scale_linetype_manual(name=NULL,values=c("LFT percent positive"="solid","qPCR percent positive"="solid",
                                                     "True prevalence"="solid","Regular screening"="dotted"))+
            ylab("Prevalence or Incidence (%)") + 
            xlab("Time (days)") + 
            ggtitle(paste0("True prevalence ranges from ", signif(prev_stats[1],3), "% to ", signif(prev_stats[2],3),"%")) +
            theme_classic() + 
            theme(legend.position=c(0.75,0.75),legend.title=element_blank(),axis.title.x = element_blank()) +
            scale_y_continuous(expand=c(0,0))+
            theme_overall
        
        inc_stats <-inc_dat %>% filter(inf_time %in% range(prev_dat$t)) %>% mutate(x=per_cap_inc*10000) %>% pull(x)
        
        
        p_inc <- ggplot() + 
            geom_bar(data=inc_dat %>% filter(inf_time >= burnin),aes(x=inf_time,y=(per_cap_inc*100)*inc_yscale),fill="grey40",
                     stat="identity",alpha=0.5,size=0.1) + 
            scale_x_continuous(
                limits=c(burnin-0.5,max(times)+0.5),
                breaks=seq(burnin,max(times),by=test_frequency*2),
                labels=seq(burnin,max(times),by=test_frequency*2)-burnin,
                expand=c(0,0)) +
            ggtitle(paste0("Incidence ranges from ", signif(inc_stats[1],3), " to ", signif(inc_stats[2],3)," per 10,000 people")) +
            ylab("Incidence (%)") + 
            xlab("Time (days)") + 
            theme_classic() + 
            theme(legend.position=c(0.75,0.75),legend.title=element_blank(),axis.title.x = element_blank()) +
            scale_y_continuous(expand=c(0,0)) +
            theme_overall
        p_inc/p_prev/p_ppv/p_npv
        
    },height=1200,width=800)
    
    output$linelist_plot <- renderPlot({
        print(gc())
        parameters <- reactiveValuesToList(input)
        growth_rate <- parameters$growth_rate
        ## Simulation duration
        ## Cumulative incidence in the focal period (ie. between burnin and max(times))
        cumu_inc <- parameters$cumu_inc
        ## Population size
        N <- parameters$N
        ## For how long are people true infections?
        duration_pos <- parameters$duration_pos
        
        ## Burning for simulation
        burnin <- duration_pos*2
        sim_duration <- parameters$sim_duration
        times <- 1:(burnin + sim_duration)
        
        strategy_start <- parameters$strategy_start_time + burnin
        strategy_start <- min(strategy_start, max(times))
        ## =====================
        ## Time-varying maximum sensitivity over infection ie. within-host model
        ## =====================
        max_sens <-  storage_for_plots()$max_sens
        
        ## =====================
        ## Specificity and test frequency
        ## =====================
        ## PCR spec
        pcr_spec <- parameters$pcr_spec
        ## LFT spec
        lft_spec <- parameters$lft_spec
        
        ## LFT screening approach
        test_frequency <- parameters$test_frequency
        isolation_period <- parameters$isolation_period
        
        
        ## When to plot crosds section?
        
        testing_times <- seq(strategy_start,max(times),by=test_frequency)
        t_choose_index <- min(parameters$t_choose_index,length(testing_times))
        t_choose <- testing_times[t_choose_index]
        
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
        
        ## Subsample linelist part of simulation for speed
        subsamp_frac <- parameters$subsamp_frac
        
        
        ## Create linelist entry for each infection
        sim_dat_linelist <- inc_dat %>% dplyr::select(-per_cap_inc) %>% uncount(inc) %>% mutate(i=1:n())
        
        ## N infected overall is higher than specified cumu_inc because of burn in period
        n_infected <- nrow(sim_dat_linelist)
        
        ## Subsample the linelist to improve simulation speed
        sim_dat_linelist <- sim_dat_linelist %>% 
            sample_frac(subsamp_frac) %>%  mutate(i=1:n())
        print(nrow(sim_dat_linelist))
        print(gc())
        ## =====================
        ## Individual entries
        ## =====================
        ## Need entry for each possible time for each individual
        infection_statuses_linelist <- expand_grid(i=1:nrow(sim_dat_linelist),
                                                   t=times[times >= burnin]) %>%  ## We only care about status after the burn in period
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
        
        ## Add entries for not infected individuals
        negative_statuses_linelist <- expand_grid(i = (max(infection_statuses_linelist$i)+1):floor(N*subsamp_frac), true_pos=0, t=times[times >= burnin]) %>% 
            mutate(false_pcr_pos=rbinom(n(),size=1,prob=1-pcr_spec),
                   true_pcr_pos=0,
                   true_lft_pos=0,
                   false_lft_pos = (1-true_pos)*rbinom(n(), size=1, prob=1-lft_spec))
        
        infection_statuses_linelist <- bind_rows(infection_statuses_linelist, negative_statuses_linelist)
        
        ## Get overall positive result
        infection_statuses_linelist <- infection_statuses_linelist %>% 
            mutate(pcr_pos = as.numeric(true_pcr_pos | false_pcr_pos),
                   lft_pos=as.numeric(true_lft_pos | false_lft_pos))
        
        ## =====================
        ## Assess prevalence and get repeat screening times
        ## =====================
        ## Get prevalence and find time when prevalence is around the specified start prevalence
        prev_dat <- infection_statuses_linelist %>% 
            group_by(t) %>% 
            summarize(y=sum(true_pos)/n())
        
        ## =====================
        ## Simulate repeated screening and pull out first detection date
        ## =====================
        ## Now let's see what happens when we test every 3 days by LFT and remove positives
        infection_statuses_linelist <- infection_statuses_linelist %>% 
            mutate(tested=ifelse(t %in% testing_times,1,0))
        
        ## If test positive by LFT, remove
        first_detection_dat <- infection_statuses_linelist %>% group_by(i) %>% 
            filter(tested == 1 & lft_pos == 1) %>% 
            filter(t == min(t)) %>% 
            mutate(first_detection=t) %>% 
            mutate(first_detection_correct=ifelse(false_lft_pos == 1,"False positive","True positive")) %>%
            dplyr::select(i, first_detection,first_detection_correct)
        
        infection_statuses_surveillance_linelist <- infection_statuses_linelist %>% left_join(first_detection_dat)
        ## =====================
        ## Now tally prevalence etc over time
        ## =====================
        indiv_key <-  infection_statuses_linelist %>% 
            filter(!is.na(inf_time)) %>%
            filter(inf_time >= burnin) %>%
            dplyr::select(inf_time, i) %>%
            distinct() %>%
            arrange(inf_time) %>%
            mutate(id=1:n())
        
        ## =====================
        ## PCR
        ## =====================
        ## Get PCR status for each day of cross-sectional surveillance
        ## True negative if tested and correctly negative
        ## False positive if tested and incorrectly positive
        ## False negative if tested and incorrectly negative
        ## True positive if tested and correctly positive
        ## Recovered if more than duration_pos after infection
        infection_statuses_pcr <- infection_statuses_linelist %>% 
            filter(t >= burnin) %>%
            mutate(Status=NA,
                   Status=ifelse(true_pos == 0 & pcr_pos == 0,"True negative", Status),
                   Status=ifelse(true_pos == 0 & pcr_pos == 1,"False positive", Status),
                   Status=ifelse(true_pos == 1 & pcr_pos == 0,"False negative",Status),
                   Status=ifelse(!is.na(inf_time) & t >= inf_time + duration_pos, "Recovered", Status),
                   Status=ifelse(true_pos == 1 & pcr_pos == 1, "True positive", Status))
        infection_statuses_pcr$Status <- factor(infection_statuses_pcr$Status, levels=status_levels)
        
        ## Pull out positives only
        positive_infection_statuses_pcr <-  infection_statuses_pcr %>%
            filter(!is.na(inf_time)) %>%
            left_join(indiv_key)
        
        ## Tally PCR positives
        pcr_surveillance_tally <- infection_statuses_pcr %>% 
            filter(t == t_choose) %>%
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
        ## =====================
        ## LFT
        ## =====================
        ## True negative if tested and correctly negative
        ## False positive if tested and incorrectly positive
        ## False negative if tested and incorrectly negative
        ## True positive if tested and correctly positive
        ## Recovered if more than duration_pos after infection
        infection_statuses_lft <- infection_statuses_linelist %>%
            filter(t >= burnin) %>%
            mutate(Status=NA,
                   Status=ifelse(true_pos == 0 & lft_pos == 0,"True negative", Status),
                   Status=ifelse(true_pos == 0 & lft_pos == 1,"False positive", Status),
                   Status=ifelse(true_pos == 1 & lft_pos == 0,"False negative",Status),
                   Status=ifelse(!is.na(inf_time) & t >= inf_time + duration_pos, "Recovered", Status),
                   Status=ifelse(true_pos == 1 & lft_pos == 1, "True positive", Status))
        ## Pull out positives only
        positive_infection_statuses_lft <-  infection_statuses_lft %>%
            filter(!is.na(inf_time)) %>%
            left_join(indiv_key)
        
        ## Tally LFT positives
        lft_surveillance_tally <- infection_statuses_lft %>% 
            filter(t == t_choose) %>%
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
        ## =====================
        ## LFT repeated testing strategy
        ## =====================
        indiv_key_surveillance <-  infection_statuses_surveillance_linelist %>% 
            filter(!is.na(inf_time)) %>%
            filter(inf_time >= burnin) %>%
            dplyr::select(inf_time, i) %>%
            distinct() %>%
            arrange(inf_time) %>%
            mutate(id=1:n())
        
        ## Not tested if tested == 0
        ## True negative if tested and correctly negative
        ## False positive if tested and incorrectly positive
        ## False negative if tested and incorrectly negative
        ## True positive if tested and correctly positive
        ## Post isolation if was detected, and now recovered
        ## Not detected if was infected, it's after the positive period and was never detected
        ## Correctly isolating if was true positive upon detection and within isolation period
        ## Incorrectly isolating if was false positive upon detection and within isolation period
        infection_statuses_lft_surveillance <- infection_statuses_surveillance_linelist %>%
            filter(t>= burnin) %>%
            mutate(Status=NA,
                   Status=ifelse(tested ==0, "Not tested",Status),
                   Status=ifelse(tested==1 & true_pos == 0 & lft_pos == 0,"True negative", Status),
                   Status=ifelse(tested==1 & true_pos == 0 & lft_pos == 1,"False positive", Status),
                   Status=ifelse(tested==1 & true_pos == 1 & lft_pos == 0,"False negative",Status),
                   Status=ifelse(tested==1 & true_pos == 1 & lft_pos == 1, "True positive", Status),
                   Status=ifelse(!is.na(first_detection) & t >= first_detection + isolation_period, "Post isolation", Status),
                   Status=ifelse(!is.na(inf_time) & t >= inf_time + duration_pos & is.na(first_detection), "Not detected", Status),
                   Status=ifelse(!is.na(first_detection) & t > first_detection & first_detection_correct=="True positive" & t <= first_detection + isolation_period,"Correctly\n isolating",Status),
                   Status=ifelse(!is.na(first_detection) & t > first_detection & first_detection_correct=="False positive" & t <= first_detection + isolation_period,"Incorrectly\n isolating",Status)
            )
        ## Pull out only those who have been infected post burn in
        positive_infection_statuses_lft_surveillance <-  infection_statuses_lft_surveillance %>%
            filter(!is.na(inf_time), inf_time >= burnin) %>%
            left_join(indiv_key_surveillance)
        
        
        ## Tally
        lft_repeat_tally <- infection_statuses_lft_surveillance %>% 
            filter(t == t_choose) %>%
            group_by(Status) %>%
            tally() %>% 
            full_join(tibble(Status=status_levels)) %>%
            mutate(n=ifelse(is.na(n),0,n)) %>%
            pivot_wider(names_from=Status,values_from=n) %>%
            mutate(`True negative/\nPost infection`=`True negative` + `Post isolation`) %>%
            mutate(i=1) %>%
            pivot_longer(-i) %>%
            filter(!(name %in% c("True negative","Post isolation","Not detected"))) %>%
            rename(Status=name) %>%
            mutate(value=value/subsamp_frac)
        
        lft_repeat_tally$Status <- factor(lft_repeat_tally$Status,
                                          levels=status_levels) 
        
        ## Tally over time
        lft_repeat_tally2 <- infection_statuses_lft_surveillance %>% 
            group_by(Status, t) %>%
            tally() %>% 
            full_join(tibble(Status=status_levels)) %>%
            mutate(n=ifelse(is.na(n),0,n)) %>%
            mutate(n=floor(n/subsamp_frac)) ## Rescale numbers, as was only subsampling before
        
        ## How many individuals not in isolation who should be?
        lft_repeat_tally2_not_isolating <- infection_statuses_lft_surveillance %>% filter(true_pos == 1, Status != "True positive", 
                                                                                          Status!="Correctly\n isolating",Status != "Post isolation", 
                                                                                          Status != "Incorrectly\n isolating") %>%
            group_by(t) %>% tally() %>% mutate(Status="Incorrectly\n not isolating")  %>%
            mutate(n=floor(n/subsamp_frac)) ## Rescale numbers, as was only subsampling before
        
        lft_repeat_tally2 <- bind_rows(lft_repeat_tally2, lft_repeat_tally2_not_isolating)
        
        lft_repeat_tally2$Status <- factor(lft_repeat_tally2$Status,
                                           levels=status_levels)
        
        bind_rows(lft_surveillance_tally, pcr_surveillance_tally, lft_repeat_tally) %>% filter(!(Status %in% c("True negative/\nPost infection","Recovered","True negative/\nRecovered","Not tested"))) %>%
            pull(value) %>% max() -> ymax_counts
        ymax_counts <- ceiling(ymax_counts/1000)*1000
        
        p_pcr_left <- ggplot(pcr_surveillance_tally %>% filter(Status != "True negative/\nRecovered")) +
            geom_bar(aes(x=Status,y=value,fill=Status),stat="identity",col="grey10") +
            scale_fill_manual(values=linelist_color_key) +
            ylab("Count") +
            theme_classic() +
            scale_y_continuous(expand=c(0,0),limits=c(0,ymax_counts),breaks=seq(0,ymax_counts,by=floor(ymax_counts/10))) +
            theme(legend.position="none")+
            theme_overall+
            theme(axis.text.x=element_text(size=12,angle=45,hjust=1),
                  axis.text.y=element_text(size=12))
        p_lft_left <- ggplot(lft_surveillance_tally %>% filter(Status != "True negative/\nRecovered")) +
            geom_bar(aes(x=Status,y=value,fill=Status),stat="identity",col="grey10") +
            scale_fill_manual(values=linelist_color_key) +
            ylab("Count") +
            theme_classic() +
            scale_y_continuous(expand=c(0,0),limits=c(0,ymax_counts),breaks=seq(0,ymax_counts,by=floor(ymax_counts/10))) +
            theme(legend.position="none")+
            theme_overall+
            theme(axis.text.x=element_text(size=12,angle=45,hjust=1),
                  axis.text.y=element_text(size=12))
        p_lft_repeat_left <- ggplot(lft_repeat_tally %>% filter(!(Status %in% c("True negative/\nPost infection","Recovered","True negative/\nRecovered","Not tested","Incorrectly\n not isolating")))) +
            geom_bar(aes(x=Status,y=value,fill=Status),stat="identity",col="grey10") +
            scale_fill_manual(values=linelist_color_key) +
            ylab("Count") +
            theme_classic() +
            scale_y_continuous(expand=c(0,0),limits=c(0,ymax_counts),breaks=seq(0,ymax_counts,by=floor(ymax_counts/10))) +
            theme(legend.position="none")+
            theme_overall +
            theme(axis.text.x=element_text(size=12,angle=45,hjust=1),
                  axis.text.y=element_text(size=12))
        
        
        p1 <- positive_infection_statuses_pcr %>% 
            filter(inf_time >= burnin) %>%
            ggplot() + 
            geom_tile(aes(x=t,y=id,fill=Status)) +
            theme_classic() +
            ggtitle("Daily PCR surveillance") +
            scale_fill_manual(values=linelist_color_key) +
            scale_y_continuous(expand=c(0,0)) +
            ylab("Individual") +
            scale_x_continuous(limits=c(burnin-0.5,max(times)+1.5),expand=c(0,0),
                               breaks=seq(burnin,max(times),by=5),
                               labels=seq(0,max(times)-burnin,by=5)) +
            xlab("Time (days)")+
            theme_overall+
            theme(legend.position="none")
        
        
        
        p2 <- positive_infection_statuses_lft %>% 
            filter(inf_time >= burnin) %>%
            ggplot() + geom_tile(aes(x=t,y=id,fill=Status)) +
            theme_classic() +
            ggtitle("Daily LFT surveillance") +
            scale_fill_manual(values=linelist_color_key) +
            scale_y_continuous(expand=c(0,0)) +
            ylab("Individual") +
            scale_x_continuous(limits=c(burnin-0.5,max(times)+1.5),expand=c(0,0),
                               breaks=seq(burnin,max(times),by=5),
                               labels=seq(0,max(times)-burnin,by=5)) +
            xlab("Time (days)")+
            theme_overall+
            theme(legend.position="none")
        
        positive_infection_statuses_lft_surveillance %>% 
            mutate(recent = inf_time >= strategy_start - 7) %>% ## Get infections within 7 days of starting strategy
            dplyr::select(id,t,inf_time, Status,recent) %>% 
            filter(recent == TRUE) %>%
            filter(inf_time == min(inf_time)) %>% ## Get earliest infections
            filter(id == min(id)) %>%
            pull(id) %>%
            unique() -> lowest_id
        
        
        p3_alt <- positive_infection_statuses_lft_surveillance %>% 
            filter(tested==1) %>%
            ggplot() + 
            geom_tile(aes(x=t,y=id,fill=Status)) +
            geom_hline(yintercept=lowest_id,col="yellow",size=1,linetype="dashed") +
            ggtitle(paste0(test_frequency,"-day LFT test and isolate")) +
            theme_classic() +
            scale_fill_manual(values=linelist_color_key)+
            scale_y_continuous(expand=c(0,0)) +
            ylab("Individual") +
            scale_x_continuous(limits=c(strategy_start-1.5,max(times)+1.5),
                               breaks=seq(strategy_start,max(times),by=5),
                               labels=seq(strategy_start,max(times),by=5)-burnin,
                               expand=c(0,0)) +
            xlab("Time (days)")+
            theme_overall 
        
        p_isolations2 <- lft_repeat_tally2 %>% 
            filter(Status == "Incorrectly\n not isolating" | Status == "Correctly\n isolating") %>% 
            filter(t >= strategy_start) %>% pivot_wider(names_from=Status,values_from=n) %>%
            mutate(ratio = `Incorrectly\n not isolating`/(`Incorrectly\n not isolating`+`Correctly\n isolating`)) %>% 
            ggplot() + geom_line(aes(x=t,y=ratio))
        
        p_isolations3 <- lft_repeat_tally2 %>%
            filter(Status == "Incorrectly\n not isolating" | Status == "Correctly\n isolating") %>% 
            ggplot() + 
            geom_bar(aes(x=t,y=n,fill=Status),stat="identity")
        
        
        
        p_isolating <- lft_repeat_tally2 %>% 
            filter(Status %in% c("Correctly\n isolating","Incorrectly\n isolating","False positive","True positive","Incorrectly\n not isolating")) %>%
            pivot_wider(names_from=Status,values_from=n,values_fill=list(n=0)) %>%
            mutate(`False positive`=ifelse(is.na(`False positive`), 0, `False positive`),
                   `True positive`=ifelse(is.na(`True positive`), 0, `True positive`)) %>%
            mutate("Correctly isolating"=`Correctly\n isolating` + `True positive`,
                   "Incorrectly isolating"=`Incorrectly\n isolating` + `False positive`,
                   "Incorrectly not\n isolating/positive"=`Incorrectly\n not isolating`) %>%
            pivot_longer(-t) %>%
            filter((name %in% c("Correctly isolating","Incorrectly isolating"))) %>% #,"Incorrectly not\n isolating/positive"))) %>%
            ggplot() +
            geom_bar(aes(x=t,y=value,fill=name),stat="identity",col="grey10") +
            theme_classic() +
            scale_fill_manual(values=linelist_color_key)+
            scale_y_continuous(expand=expansion(mult=c(0,0.05)),breaks=seq(0,1000000,by=5000)) +
            ylab("") +
            scale_x_continuous(limits=c(strategy_start-1,max(times)+0.5),
                               breaks=seq(strategy_start,max(times),by=test_frequency*2),
                               labels=seq(strategy_start,max(times),by=test_frequency*2)-burnin,
                               expand=c(0,0)) +
            xlab("Time (days)")+
            ggtitle("Number of individuals in isolation over time under routine testing") +
            theme_overall +
            theme(legend.position="none",
                  legend.title=element_blank())
        
        p_positives <- lft_repeat_tally2 %>% 
            filter(Status %in% c("Correctly\n isolating","Incorrectly\n isolating","False positive","True positive")) %>%
            filter(t %in% seq(strategy_start,max(times),by=test_frequency)) %>%
            pivot_wider(names_from=Status,values_from=n) %>%
            mutate(`False positive`=ifelse(is.na(`False positive`), 0, `False positive`),
                   `True positive`=ifelse(is.na(`True positive`), 0, `True positive`)) %>%
            mutate("Correctly isolating"=`Correctly\n isolating` + `True positive`,
                   "Incorrectly isolating"=`Incorrectly\n isolating` + `False positive`) %>%
            pivot_longer(-t) %>%
            mutate(value = ifelse(is.na(value),0,value)) %>%
            filter((name %in% c("True positive","False positive"))) %>%
            full_join(expand_grid(name=c("False positive","True positive"),t=seq(burnin,max(times),by=1),value=0)) %>%
            group_by(name,t) %>%
            summarize(value=sum(value)) %>%
            ggplot() +
            geom_bar(aes(x=t,y=value,fill=name),stat="identity",col="grey10") +
            theme_classic() +
            scale_fill_manual(values=linelist_color_key)+
            scale_y_continuous(expand=expansion(mult=c(0,0.05)),breaks=seq(0,1000000,by=2000)) +
            ylab("Count") +
            ggtitle("Number of positives over time under routine testing") +
            scale_x_continuous(limits=c(strategy_start-1,max(times)+0.5),
                               breaks=seq(strategy_start,max(times),by=test_frequency*2),
                               labels=seq(0,max(times)-strategy_start,by=test_frequency*2),
                               expand=c(0,0)) +
            xlab("Days of regular testing")+
            theme_overall +
            theme(legend.position="none",
                  legend.title=element_blank())
        
        
    
        ##################################
        ## Infection annotation text
        ##################################
        
        prev_time <- strategy_start - burnin
        
        ## Infections in week preceding screening
        precede_infections <- inc_dat %>% filter(inf_time >= strategy_start-prev_time, inf_time < strategy_start) %>% summarize(x=sum(inc)) %>% pull(x)
         
        ## Infections after screening initiated
        initiated_infections <- inc_dat %>% filter(inf_time >= strategy_start) %>% summarize(x=sum(inc)) %>% pull(x)
        
        ## Number of infections in the final week
        final_infections <- inc_dat %>% filter(inf_time >= max(inf_time) - 7) %>% summarize(x=sum(inc)) %>% pull(x)
        
        ## Number of infections in the screening period less one week
        screen_period_infections <- inc_dat %>% filter(inf_time < max(inf_time) - 7, inf_time >= strategy_start) %>% summarize(x=sum(inc)) %>% pull(x)
        
        precede_infections_text <- paste0(signif(precede_infections,2), " infections preceding the start of routine screening by up to ", prev_time, " days.")
        initiated_infections_text <- paste0(signif(initiated_infections,2)," infections after initiation of routine screening.")
        final_infections_text <- paste0(signif(final_infections,2), " infections in final week of simulation.")
        screen_infections_text <- paste0(signif(screen_period_infections,2), " infections between start of routine screening and final week.")
        all_infection_text <- paste0(precede_infections_text,"\n",
                                     initiated_infections_text,"\n",
                                     final_infections_text,"\n",
                                     screen_infections_text,"\n")
        
        p_infection_text <- ggplot() + annotate("text",x=4,y=25,size=4.5,label=all_infection_text) + theme_void()
        
        ##################################
        ## Detected annotation text
        ##################################
        ## Missed just before start
        missed_in_preceding <- infection_statuses_lft_surveillance %>%
            filter(is.na(first_detection) & !is.na(inf_time), inf_time >= strategy_start - prev_time, inf_time < strategy_start)  %>% 
            dplyr::select(i, inf_time) %>% 
            distinct() %>% 
            tally() %>%
            mutate(n=n/subsamp_frac) %>% pull(n)
        
        ## Missed during surveillance
        missed_during_surveillance <- infection_statuses_lft_surveillance %>%
            filter(is.na(first_detection) & !is.na(inf_time), inf_time >= strategy_start)  %>% 
            dplyr::select(i, inf_time) %>% 
            distinct() %>% 
            tally() %>%
            mutate(n=n/subsamp_frac) %>%
            pull(n)
        
        ## Missed before final week
        missed_before_end <- infection_statuses_lft_surveillance %>%
            filter(is.na(first_detection) & !is.na(inf_time), inf_time >= strategy_start, inf_time < max(t)-7)  %>% 
            dplyr::select(i, inf_time) %>% 
            distinct() %>% 
            tally() %>%
            mutate(n=n/subsamp_frac)%>%
            pull(n)
        
        
        ## Missed before final week
        missed_final_week <- infection_statuses_lft_surveillance %>%
            filter(is.na(first_detection) & !is.na(inf_time), inf_time >= max(t)-7)  %>% 
            dplyr::select(i, inf_time) %>% 
            distinct() %>% 
            tally() %>%
            mutate(n=n/subsamp_frac)%>%
            pull(n)
        
        precede_missed_text <- paste0(missed_in_preceding, " missed infections preceding the start of routine screening by up to ", prev_time, " days.")
        initiated_missed_text <- paste0(missed_during_surveillance," missed infections after initiation of routine screening.")
        final_missed_text <- paste0(missed_final_week, " missed infections in final week of simulation.")
        screen_missed_text <- paste0(missed_before_end, " missed infections between start of routine screening and final week.")
        all_missed_text <- paste0(precede_missed_text,"\n",
                                  initiated_missed_text,"\n",
                                  final_missed_text,"\n",
                                  screen_missed_text,"\n")
        
        p_missed_text <- ggplot() + annotate("text",x=4,y=25,size=4.5,label=all_missed_text) + theme_void()
        
        layout <- "
        AAABBB
        CCDDEE
        CCDDEE
        CCDDEE
        CCDDEE
        CCDDEE
        FFGGHH
        FFGGHH
        IIIJJJ
        IIIJJJ
        "
        
        p_indiv_strategies <- p_infection_text + 
            p_missed_text +
            (p1+theme(legend.position="none")+geom_vline(xintercept=t_choose,linetype="dashed",col="red")) +
            (p2+geom_vline(xintercept=t_choose,linetype="dashed",col="red")+
                 theme(legend.position="none",axis.text.y=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank()))+
            (p3_alt+geom_vline(xintercept=t_choose,linetype="dashed",col="red")+
                 theme(axis.text.y=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank())) +
            (p_pcr_left + theme(axis.title.y=element_blank(),axis.title.x=element_blank()))+
            (p_lft_left +  theme(axis.title.y=element_blank(),axis.title.x=element_blank()))+
            (p_lft_repeat_left + theme(axis.title.y=element_blank(),axis.title.x=element_blank()))+
            p_positives + 
            p_isolating + 
            plot_layout(guides="collect",design=layout)
        
            p_indiv_strategies
    },height=1200,width=1200)

})
