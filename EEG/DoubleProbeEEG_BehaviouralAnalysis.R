# Analysis of behavioural data from DoubleProbeEEG - SRC

#### install necessary packages -----
library("ggplot2")
library("dplyr")
library("circular")
library("pracma")
library("afex")
library("BayesFactor")


path <- "/Users/user/Desktop/Experiments/Nick/DoubleProbe/EEG/data/allEEG.csv"


EEG <- read.csv(path, header=TRUE, as.is=TRUE); EEG <- dplyr::select(EEG, -X)

EEG$rdif1 <- seq(1,dim(EEG)[1],by=1); EEG$rdif2      <- seq(1,dim(EEG)[1],by=1) #make a new column called rdif2
EEG$cond  <- seq(1,dim(EEG)[1],by=1); 
for(i in seq(24,26,by=1)) EEG[,i] <- 0
EEG <- dplyr::select(EEG, -angs_3, -angs_4, -locs_1, -locs_2, -locs_3, -locs_4, -dist_1, -dist_2, -prec_1, -prec_2)

sublist <- unique(EEG$subID)

for(i in seq(1,dim(EEG)[1]))
  if(EEG[i,11] == 1 & EEG[i,8] == 1) EEG[i,16] <- "cued1"   
for(i in seq(1,dim(EEG)[1]))
  if(EEG[i,11] == 1 & EEG[i,8] == 2) EEG[i,16] <- "cued2"
for(i in seq(1,dim(EEG)[1]))
  if(EEG[i,11] == -2) EEG[i,16] <- "neutral"

EEG <- dplyr::mutate(EEG, rdif1 = as.circular(theta_1 - rad(angs_1), units="radians")) #don't change this!!
EEG <- dplyr::mutate(EEG, rdif2 = as.circular(theta_2 - rad(angs_2), units="radians")) #don't change this!!

EEG <- dplyr::filter(EEG, time_1 != Inf & time_2 != Inf)

neutral <- vector("numeric", length = length(sublist))
cued    <- vector("numeric", length = length(sublist))
uncued  <- vector("numeric", length = length(sublist))

accs <- data.frame(neutral, cued, uncued, neutral, cued, uncued, stringsAsFactors = FALSE)
colnames(accs) <- c("n1", "c1", "u1", "n2", "c2", "u2")

count <- 1
for(i in sublist){
  temp <- dplyr::filter(EEG, subID == i)
  tempneut  <- dplyr::filter(temp, cond == "neutral")
  tempcued1 <- dplyr::filter(temp, cond == "cued1")
  tempcued2 <- dplyr::filter(temp, cond == "cued2")
  accs[i,1] <- 1/sd(tempneut[,14]);# E1precs[count,1] <- mean(abs(tempneut[,14]))
  accs[i,2] <- 1/sd(tempcued1[,14]);# E1precs[count,2] <- mean(abs(tempcued1[,14]))
  accs[i,3] <- 1/sd(tempcued2[,14]);# E1precs[count,3] <- mean(abs(tempcued2[,14]))
  accs[i,4] <- 1/sd(tempneut[,15]);# E1precs[count,4] <- mean(abs(tempneut[,15]))
  accs[i,5] <- 1/sd(tempcued2[,15]);# E1precs[count,5] <- mean(abs(tempcued2[,15]))
  accs[i,6] <- 1/sd(tempcued1[,15]);# E1precs[count,6] <- mean(abs(tempcued1[,15]))
  count <- count + 1
}

accuracy <- summarise_each(accs, funs(mean)); accuracy <- as.data.frame(t(accuracy)); colnames(accuracy)[1] <- "accuracy"
se <- function(x) sd(x)/sqrt(length(x))

#this makes data frames for each experiment for easy plotting
accuracy$condition  <- seq(1, dim(accuracy)[1], by=1);
accuracy$probe      <- seq(1, dim(accuracy)[1], by=1);
accuracy$sem        <- t(summarise_each(accs, funs(se)));

accuracy[c(1,2,3),]$probe <- "first"; accuracy[c(4,5,6),]$probe <- "second"
accuracy[c(1,4),]$condition <- "neutral";
accuracy[c(2,5),]$condition <- "cued";
accuracy[c(3,6),]$condition <- "uncued";
accuracy$condition <- as.factor(accuracy$condition); accuracy$probe <- as.factor(accuracy$probe)
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

ggplot(accuracy, aes(x=probe, y=accuracy)) + #generate ggplot, plotting probe and precision
  ggtitle("DoubleProbeEEG - accuracy") +
  geom_bar(aes(fill=condition), stat="identity", color="black", position=position_dodge()) + #create the bars and fill by cue condition
  scale_fill_manual(values = c("neutral" = "#ffffbf","cued" = "#fc8d59","uncued" = "#91bfdb")) + #manually set the colours for the bars
  xlab("Probe") + ylab("Accuracy (1/SD)") +  #set axis labels
  #plot the error bars
  geom_errorbar(aes(x=probe, ymin=accuracy-sem, ymax = accuracy+sem, fill = condition), stat="identity", position=position_dodge(0.9), width = 0.1) +
  scale_colour_manual(values = c("cued" = "#000000" , "neutral" = "#000000", "uncued" = "#000000")) + #set colour of errs manually
  theme(panel.background = element_rect(fill = "white"), panel.grid.minor = element_blank(), axis.text = element_text(size = 14))
#+ ylim(0,1.6)
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


id <- seq(1, 6*length(sublist), by = 1)
acc <- data.frame(id, stringsAsFactors = FALSE); colnames(acc)[1] <- "id"

acc$accuracy   <- seq(1, dim(accs)[2], by = 1)
acc$probe      <- seq(1, dim(accs)[2], by = 1)
acc$cue        <- seq(1, dim(accs)[2], by = 1)

count <-  1; sub <- 1 #we have to reset subject numbers starting from 1 so it doesn't interfere with anova function later
for(i in sublist){
  temp <- dplyr::filter(EEG, subID == i)
  tempneut  <- dplyr::filter(temp, cond == "neutral")
  tempcued1 <- dplyr::filter(temp, cond == "cued1")
  tempcued2 <- dplyr::filter(temp, cond == "cued2")
  acc[count,2] <- 1/sd(tempneut[,14]) ; acc[count,1] <- i; count <- count + 1
  acc[count,2] <- 1/sd(tempcued1[,14]); acc[count,1] <- i; count <- count + 1
  acc[count,2] <- 1/sd(tempcued2[,14]); acc[count,1] <- i; count <- count + 1
  acc[count,2] <- 1/sd(tempneut[,15]) ; acc[count,1] <- i; count <- count + 1
  acc[count,2] <- 1/sd(tempcued2[,15]); acc[count,1] <- i; count <- count + 1
  acc[count,2] <- 1/sd(tempcued1[,15]); acc[count,1] <- i; count <- count + 1
  sub <- sub + 1
} #put accs into long format

acc$probe[acc$probe %in% c(1,2,3)] <- "first"; acc$probe[acc$probe %in% c(4,5,6)] <- "second"
acc$cue[acc$cue %in% c(1,4)] <- "neutral"
acc$cue[acc$cue %in% c(2,5)] <- "cued"    
acc$cue[acc$cue %in% c(3,6)] <- "uncued"  
acc$cue <- as.factor(acc$cue); acc$probe <- as.factor(acc$probe); acc$id <- as.factor(acc$id)

aov <- aov_ez('id', 'accuracy', data = acc, within = c("probe", "cue"))
nice(aov, es = 'pes')
bf_aov <- anovaBF(accuracy ~id + cue*probe, data = acc, whichRandom = 'id', iterations = 200000, progress = TRUE, whichModels = "withmain")
bf_aov[1] #main effect of probe
bf_aov[2] #main effect of cue
lmbf_full  <- lmBF(accuracy ~ id + cue + probe + cue:probe #+
                   #id:cue + id:probe
                   ,data = acc, whichRandom = 'id', iterations = 200000)
lmbf_mains <- lmBF(accuracy ~ id + cue + probe #+
                   #id:cue + id:probe
                   ,data = acc, whichRandom = 'id', iterations = 200000)
lmbf_full/lmbf_mains #evidence for model with cue:probe interaction, vs. just main effects

acc_NeutralUncued <- dplyr::filter(acc, cue == "neutral" | cue == "uncued")
acc_NeutralUncued$id <- as.factor(acc_NeutralUncued$id)
acc_NeutralUncued$cue <- as.character(acc_NeutralUncued$cue);     acc_NeutralUncued$cue <- as.factor(acc_NeutralUncued$cue)
acc_NeutralUncued$probe <- as.character(acc_NeutralUncued$probe); acc_NeutralUncued$probe <- as.factor(acc_NeutralUncued$probe)
aov_NU <- aov_ez('id', 'accuracy', data = acc_NeutralUncued, within = c("probe", "cue"))
nice(aov_NU, es = 'pes')

bf_NU <- anovaBF(accuracy ~id + cue*probe, data = acc_NeutralUncued, whichRandom = 'id', iterations = 200000, progress = TRUE, whichModels = "withmain")
lmbf_NU_full <- lmBF(accuracy ~ id + cue + probe + cue:probe,
                        data = acc_NeutralUncued, whichRandom = 'id', iterations = 200000)
lmbf_NU_main <- lmBF(accuracy ~ id + cue + probe,
                        data = acc_NeutralUncued, whichRandom = 'id', iterations = 200000)
lmbf_NU_full/lmbf_NU_main # evidence for interaction of cue and probe in neutral vs. uncued trials in Experiment 1

acc_NeutralCued <- dplyr::filter(acc, cue == "neutral" | cue == "cued")
acc_NeutralCued$id <- as.factor(acc_NeutralCued$id)
acc_NeutralCued$cue <- as.character(acc_NeutralCued$cue);     acc_NeutralCued$cue <- as.factor(acc_NeutralCued$cue)
acc_NeutralCued$probe <- as.character(acc_NeutralCued$probe); acc_NeutralCued$probe <- as.factor(acc_NeutralCued$probe)
aov_NC <- aov_ez('id', 'accuracy', data = acc_NeutralCued, within = c("probe", "cue"))
nice(aov_NC, es = 'pes')
bf_NC <- anovaBF(accuracy ~id + cue*probe, data = acc_NeutralCued, whichRandom = 'id', iterations = 200000, progress = TRUE, whichModels = "withmain")
# using bf_E1_NC[7]/bf_E1_NC[4] as the probe:cue interaction, get BF 1.09. probe:cue + id against id (bf_E1_NC[3]) gives BF = 0.593

bf_t_NC1 <- ttestBF(x = accs$n1, y = accs$c1, paired = TRUE); bf_tE1_NC2 <- ttestBF(x = accs$n2, y = accs$c2, paired = TRUE)
t_NC1    <- t.test(x = accs$n1, y = accs$c1, paired = TRUE) ; tE1_NC2    <- t.test( x = accs$n2, y = accs$c2, paired = TRUE)
bf_t_NU1 <- ttestBF(x = accs$n1, y = accs$u1, paired = TRUE); bf_tE1_NU2 <- ttestBF(x = accs$n2, y = accs$u2, paired = TRUE)
t_NU1    <- t.test(x = accs$n1, y = accs$u1, paired = TRUE) ; tE1_NU2    <- t.test( x = accs$n2, y = accs$u2, paired = TRUE)
