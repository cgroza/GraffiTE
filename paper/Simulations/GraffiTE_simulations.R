####################################### load libraries $###################################################################

library(ggplot2)
library(cowplot)
library(reshape2)

######################################## data prep #######################################################################
# for each directory, need to remove duplicates due to the bedtools intersect with TE or mappability tracks:
# for dir in $folders ; do sort -u -k20,20 $dir/3_TSD_search/$dir.intersect.csv >> output.results.csv; done 
# the only cases where they can be more than 40 lines are the false positives not in a bgSV

####################################### load data ########################################################################
svimSim <- read.table("~/Documents/GraffiTE/sim/sim_ouputs/svim.sim.100.chr22.results.csv")
snifOntSim <- read.table("~/Documents/GraffiTE/sim/sim_ouputs/snif-ont.sim.100.chr22.results.csv")
svimsnifOnt <- read.table("~/Documents/GraffiTE/sim/sim_ouputs/svim-snif-ont.sim.100.chr22.results.csv")
snifHifiSim <- read.table("~/Documents/GraffiTE/sim/sim_ouputs/snif-hifi.sim.100.chr22.results.csv")
svimsnifHifi <- read.table("~/Documents/GraffiTE/sim/sim_ouputs/svim-snif-hifi.sim.100.chr22.results.csv")

####################################### define functions ################################################################

####################################### 
## general function to format tables ##
#######################################
# name columns
formatSim <- function(sim){
  names(sim)<-c("chr", "start", "end", "consID", "suppNB", "supp_vec",
                  "consSvlen", "svimSvLen", "svimID", "unkSvLen", "RMhits",
                  "RMfrags", "RMlen", "RMname", "RMclass", "RMcov", "L15PINV",
                  "SVAonly", "simSvLen", "simID", "simType",  "simSvLensigned",
                  "test", "TEoverlap", "k24Map")

# here we mark as negative if GraffiTE gives nhits > 1 and it is not a 5PINV.
test1hit<-ifelse(sim$RMhits > 1 & (is.na(sim$L15PINV) | sim$L15PINV == "None") & sim$test == "FP", "TN", 
                 ifelse(sim$RMhits > 1 & (is.na(sim$L15PINV) | sim$L15PINV == "None") & sim$test == "TP", "FN", 
                        sim$test)
                 )

RMclass2<-ifelse(sim$RMhits == 1 & grepl("AluY", sim$RMname), "SINE/AluY", 
                   ifelse(sim$RMhits == 1 & sim$RMclass == "SINE/Alu", "SINE/Alu_other", sim$RMclass))

# attach the new test column to the main table
sim<-cbind(sim, test1hit, RMclass2)
# filter further on call >= 250 as positives
testLen<-ifelse(sim$test1hit == "FP" & abs(sim$svimSvLen) < 250, "TN",
                ifelse(sim$test1hit == "TP" & abs(sim$svimSvLen) < 250, "FN", sim$test1hit))
sim<-cbind(sim, testLen)
# reorder the TE classes to make it prettier in graphs
sim$RMclass2<-factor(sim$RMclass2, levels = c(""))
# 
sim$k24Map[sim$k24Map == "."] <- NA
sim$k24Map<-as.numeric(sim$k24Map)
return(sim)
}


#####################################
# function for F1 score calculation #
#####################################
F1sim<- function(testcol, name){
TP<-length(testcol[testcol == "TP"])
TN<-length(testcol[testcol == "TN"])
FP<-length(testcol[testcol == "FP"])
FN<-length(testcol[testcol == "FN"])
TPR<-TP/(TP+FN) # Recall
TNR<-TN/(TN+FP) # Specificity
PPV<-TP/(TP+FP) # Precision, Positive Predictive Value
F1<-(2*PPV*TPR)/(PPV+TPR)
outrow<-c(TP, TN, FP, FN, TPR, TNR, PPV, F1)
names(outrow)<-c("TP", "TN", "FP", "FN", "TPR", "TNR", "PPV", "F1")
return(outrow)
}




####################################### analyses  #######################################################################

# create enhanced tables
Sv<-formatSim(svimSim)
SnO<-formatSim(snifOntSim)
SnSvO<-formatSim(svimsnifOnt)
SnH<-formatSim(snifHifiSim)
SnSvH<-formatSim(svimsnifHifi)

# store scores for different filters and modules:
scores<-as.data.frame(
          cbind(
            c("svim", "svim", "svim",
              "sniffles2-ont", "sniffles2-ont", "sniffles2-ont",
              "svim-sniffles2-ont", "svim-sniffles2-ont",  "svim-sniffles2-ont",
              "sniffles2-hifi", "sniffles2-hifi", "sniffles2-hifi",
              "svim-sniffles2-hifi", "svim-sniffles2-hifi", "svim-sniffles2-hifi"),
            c("all", "1hit","1hit+250bp", "all", "1hit","1hit+250bp" ,"all", "1hit" ,"1hit+250bp","all", "1hit","1hit+250bp","all", "1hit","1hit+250bp"),
          rbind(F1sim(Sv$test, "test"),
            F1sim(Sv$test1hit, "test1hit"),
            F1sim(Sv$testLen, "1hit+250bp"),
            F1sim(SnO$test, "test"),
            F1sim(SnO$test1hit, "test1hit"),
            F1sim(SnO$testLen, "1hit+250bp"),
            F1sim(SnSvO$test, "test"),
            F1sim(SnSvO$test1hit, "test1hit"),
            F1sim(SnSvO$testLen, "1hit+250bp"),
            F1sim(SnH$test, "test"),
            F1sim(SnH$test1hit, "test1hit"),
            F1sim(SnH$testLen, "1hit+250bp"),
            F1sim(SnSvH$test, "test"),
            F1sim(SnSvH$test1hit, "test1hit"),
            F1sim(SnSvH$testLen, "1hit+250bp")
            )
          )
        )
names(scores)<-c("sim", "filter", "TP", "TN", "FP", "FN", "TPR", "TNR", "PPV", "F1")
# make the table ggplot-compatible
m.scores<-melt(scores, id.vars = c("sim", "filter", "TP", "TN", "FP", "FN"), value.name = "score", variable.name = "type")
# reorder factors
m.scores$sim<-factor(m.scores$sim, levels = c("svim", "sniffles2-ont", "sniffles2-hifi", "svim-sniffles2-ont", "svim-sniffles2-hifi"))
m.scores$filter<-factor(m.scores$filter, levels = c("all", "1hit", "1hit+250bp"))
# make human friendly labels
m.scores$type<-factor(m.scores$type, labels = c("Recall", "Specificity", "Precision", "F1-score"))
# calculate mean F1 and s.d.
round(mean(as.numeric(scores[scores$filter == "all",]$F1)),3)
round(sd(as.numeric(scores[scores$filter == "all",]$F1)),3)
# calculate mean recall and s.d.
round(mean(as.numeric(scores[scores$filter == "all",]$TPR)),3)
round(sd(as.numeric(scores[scores$filter == "all",]$TPR)),3)
# plot
ggplot(m.scores)+
  geom_bar(aes(y = as.numeric(score), x = sim, fill = type),stat = "identity", position = position_dodge2())+
  scale_fill_manual(values = c('#91cd55', '#80aaaf', '#f17546', '#dd0028'))+
  facet_wrap(~filter, nrow = 3)+
  scale_y_continuous(breaks = c(0,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00))+
  ylab("")+
  xlab("simulation")+
  coord_cartesian(ylim = c(0.7,1))+
  theme_minimal_hgrid()