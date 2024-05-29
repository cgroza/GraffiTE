####################################### load libraries $###################################################################

library(ggplot2)
library(cowplot)
library(reshape2)
library(dplyr)

####################################### load data ########################################################################

h_all <- read.table("human_sim100.eval")
m_all <-read.table("maize_100.eval")

####################################### define functions ################################################################

####################################### 
## general function to format tables ##
#######################################
# name columns
formatSim <- function(sim, len=250){
  names(sim)<-c("chr", "start", "consID", "REF", "ALT", "suppNB", "supp_vec",
                  "consSvlen", "svimSvLen", "svimID", "unkSvLen", "RMhits",
                  "RMfrags", "RMlen", "RMname", "RMclass", "RMcov", "L15PINV",
                  "SVAonly", "simSvLen", "simID", "simType",  "simSvLensigned",
                  "raw", "replicates", "mode")
# here we mark as negative if GraffiTE gives nhits > 1 and it is not a 5PINV.
one_hit<-ifelse(sim$RMhits > 1 & (is.na(sim$L15PINV) | sim$L15PINV == "None") & sim$raw == "FP", "TN", 
                 ifelse(sim$RMhits > 1 & (is.na(sim$L15PINV) | sim$L15PINV == "None") & sim$raw == "TP", "FN", 
                        sim$raw)
                 )

RMclass2<-ifelse(sim$RMhits == 1 & grepl("AluY", sim$RMname), "SINE/AluY", 
                   ifelse(sim$RMhits == 1 & sim$RMclass == "SINE/Alu", "SINE/Alu_other", sim$RMclass))

# attach the new test column to the main table
sim<-cbind(sim, one_hit, RMclass2)
# filter further on call >= len as positives
one_hit_len<-ifelse(sim$one_hit == "FP" & abs(sim$svimSvLen) < len, "TN",
                ifelse(sim$one_hit == "TP" & abs(sim$svimSvLen) < len, "FN", sim$one_hit))
sim<-cbind(sim, one_hit_len)
# reorder the TE classes to make it prettier in graphs
#sim$RMclass2<-factor(sim$RMclass2, levels = c(""))
#
simF<-melt(sim, measure.vars = c("raw", "one_hit", "one_hit_len"), variable.name = "filter", value.name = "test")
return(simF)
}


#####################################
# function for F1 score calculation #
#####################################

# not in use but keep for archive
# F1sim<- function(testcol, name){
# TP<-length(testcol[testcol == "TP"])
# TN<-length(testcol[testcol == "TN"])
# FP<-length(testcol[testcol == "FP"])
# FN<-length(testcol[testcol == "FN"])
# TPR<-TP/(TP+FN) # Recall
# TNR<-TN/(TN+FP) # Specificity
# PPV<-TP/(TP+FP) # Precision, Positive Predictive Value
# F1<-(2*PPV*TPR)/(PPV+TPR)
# outrow<-c(TP, TN, FP, FN, TPR, TNR, PPV, F1)
# names(outrow)<-c("TP", "TN", "FP", "FN", "TPR", "TNR", "PPV", "F1")
# return(outrow)
# }

F1sim_v2<-function(table, len=250){
  Ftable<-formatSim(table, len = len)
  Ftable %>% group_by(mode, filter) %>%
    summarise(
      mode = unique(mode),
      TP = length(test[test == "TP"]),
      TN = length(test[test == "TN"]),
      FP = length(test[test == "FP"]),
      FN = length(test[test == "FN"]),
    ) -> out
  out$TPR = out$TP/(out$TP+out$FN) # Recall
  out$TNR = out$TN/(out$TN+out$FP) # Specificity
  out$PPV = out$TP/(out$TP+out$FP) # Precision, Positive Predictive Value
  out$F1 = (2*out$PPV*out$TPR)/(out$PPV+out$TPR)
  return(out)
  }

####################################### analyses  #######################################################################

### human ###

human<-F1sim_v2(h_all)

# plot test results
human$filter<-factor(human$filter, levels = c("raw", "one_hit", "one_hit_len"), labels = c("raw", "1 hit", "1 hit; >250bp"))
human$mode<-factor(human$mode, levels = c("GT-sv", "GT-sn-ONT", "GT-sn-HIFI", "GT-svsn-ONT", "GT-svsn-HIFI"))
hres<-ggplot(human, aes(x = as.numeric(PPV), y = as.numeric(TPR), col = filter, shape = mode))+
  geom_line(aes(group = mode), linetype = "dashed", color = "grey")+
  geom_point(size = 2)+
  scale_color_manual(values = c("salmon", "red2", "firebrick"))+
  scale_shape_manual(values = c(16, 17, 15, 2, 0))+
  coord_cartesian(xlim = c(0.45, 1), ylim = c(0.65, 1))+
  xlab("precision")+
  ylab("recall")+
  theme_classic()

# sim landscape distro: ### FIX!

h_all_Fall<-formatSim(h_all)
# keep only the sim info of 1 mode, because these are the same for each mode.
h_all_F<-h_all_Fall[h_all_Fall$mode == "GT-sv",]
h_all_F$simClass<-ifelse(grepl("AluY", h_all_F$simType),"AluY",
                 ifelse(grepl("L1", h_all_F$simType), "LINE1", 
                        ifelse(grepl("SVA", h_all_F$simType),"SVA","bgSV")))

h_all_F$simClass<-factor(h_all_F$simClass, levels = c("AluY", "LINE1", "SVA", "bgSV"))
h_size_plot<-ggplot(h_all_F)+
  geom_density(aes(x = simSvLen, fill = simClass), alpha = .5)+
  scale_fill_manual(values = c("royalblue", "steelblue1", "lightblue", "grey50"))+
  geom_vline(xintercept = 250, col = "firebrick", linetype = "dashed")+
  scale_x_log10()+
  xlab("simulated variants length")+
  theme_classic()


### MAIZE ###

### maize ###
maize<-F1sim_v2(m_all, len = 500)

# plot test results
maize$filter<-factor(maize$filter, levels = c("raw", "one_hit", "one_hit_len"), labels = c("raw", "1 hit", "1 hit; >500bp"))
maize$mode<-factor(maize$mode, levels = c("GT-sv", "GT-sn-ONT", "GT-sn-HIFI", "GT-svsn-ONT", "GT-svsn-HIFI"))
mres<-ggplot(maize, aes(x = as.numeric(PPV), y = as.numeric(TPR), col = filter, shape = mode))+
  geom_line(aes(group = mode), linetype = "dashed", color = "grey")+
  geom_point(size = 2)+
  scale_color_manual(values = c("salmon", "red2", "firebrick"))+
  scale_shape_manual(values = c(16, 17, 15, 2, 0))+
  coord_cartesian(xlim = c(0.45, 1), ylim = c(0.65, 1))+
  xlab("precision")+
  ylab("recall")+
  theme_classic()

# maize
m_all_Fall<-formatSim(m_all)
# keep only the sim info of 1 mode, because these are the same for each mode.
m_all_F<-m_all_Fall[m_all_Fall$mode == "GT-sv",]

simClass2<-ifelse(grepl("Copia", m_all_F$simType),"LTR/Copia",
                  ifelse(grepl("Gypsy", m_all_F$simType), "LTR/Ty3", 
                         ifelse(grepl("DTA", m_all_F$simType),"TIR","bgSV")))
m_all_F$simClass<-simClass2
m_all_F$simClass<-factor(m_all_F$simClass, levels = c("LTR/Copia", "LTR/Ty3", "TIR", "bgSV"))
m_size_plot<-ggplot(m_all_F)+
  geom_density(aes(x = simSvLen, fill = simClass), alpha = .5)+
  scale_fill_manual(values = c("green4", "green1", "salmon", "grey50"))+
  geom_vline(xintercept = 500, col = "firebrick", linetype = "dashed")+
  scale_x_log10()+
  xlab("simulated variants length")+
  theme_classic()

sims<-cowplot::plot_grid(h_size_plot, m_size_plot, hres, mres)
