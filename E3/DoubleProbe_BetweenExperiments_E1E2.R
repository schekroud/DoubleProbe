# Between-Experiments Analysis of E1 and E2 accuracy effects - SRC

#### install necessary packages -----
library("ggplot2")
library("dplyr")
library("circular")
library("pracma")
library("afex")
library("BayesFactor")


path1 <- "/Users/user/Desktop/DoubleProbe3/data/E1.csv"
path2 <- "/Users/user/Desktop/DoubleProbe3/data/E2.csv"

E1 <- read.csv(path1, header=TRUE, as.is=TRUE); E1 <- dplyr::select(E1, -X)
E2 <- read.csv(path2, header=TRUE, as.is=TRUE); E2 <- dplyr::select(E2, -X)


E1$rdif1 <- seq(1,dim(E1)[1],by=1); E1$rdif2      <- seq(1,dim(E1)[1],by=1) #make a new column called rdif2
E1$cond  <- seq(1,dim(E1)[1],by=1); E1$experiment <- seq(1,dim(E1)[1],by=1)
for(i in seq(24,26,by=1)) E1[,i] <- 0
E1[,27] <- 1
E1 <- dplyr::select(E1, -angs_3, -angs_4, -locs_1, -locs_2, -locs_3, -locs_4, -dist_1, -dist_2, -prec_1, -prec_2)

E2$rdif1 <- seq(1,dim(E2)[1],by=1); E2$rdif2      <- seq(1,dim(E2)[1],by=1) #make a new column called rdif2
E2$cond  <- seq(1,dim(E2)[1],by=1); E2$experiment <- seq(1,dim(E2)[1],by=1)
for(i in seq(24,26,by=1)) E2[,i] <- 0
E2[,27] <- 2
E2 <- dplyr::select(E2, -angs_3, -angs_4, -locs_1, -locs_2, -locs_3, -locs_4, -dist_1, -dist_2, -prec_1, -prec_2)

sublist1 <- unique(E1$subID); sublist2 <- unique(E2$subID)

data <- dplyr::bind_rows(E1,E2)
data <- data[order(data$experiment),]

#columns: pord = 8, cues = 11
#pord 1 = probed first, 2 = probed second
#cues 1 = cued, -2 = neutral
for(i in seq(1,dim(data)[1]))
  if(data[i,11] == 1 & data[i,8] == 1) data[i,16] <- "cued1"   
for(i in seq(1,dim(data)[1]))
  if(data[i,11] == 1 & data[i,8] == 2) data[i,16] <- "cued2"
for(i in seq(1,dim(data)[1]))
  if(data[i,11] == -2) data[i,16] <- "neutral"

#write.csv(data, file = "/Users/user/Desktop/DoubleProbe3/data/Preprocessed_E1E2.csv", sep=',', eol='\n')

#data <- dplyr::mutate(data, rdif1 = as.circular(theta_1 - rad(angs_1), units="radians")) #don't change this!!
#data <- dplyr::mutate(data, rdif2 = as.circular(theta_2 - rad(angs_2), units="radians")) #don't change this!!

#if want degrees...
data <- dplyr::mutate(data, rdif1 = as.circular(resp_1 - angs_1, units="degrees")) #don't change this!!
data <- dplyr::mutate(data, rdif2 = as.circular(resp_2 - angs_2, units="degrees")) #don't change this!!


neutral <- vector("numeric", length = length(sublist1))
cued    <- vector("numeric", length = length(sublist1))
uncued  <- vector("numeric", length = length(sublist1))

E1precs <- data.frame(neutral, cued, uncued, neutral, cued, uncued, stringsAsFactors = FALSE)
E2precs <- data.frame(neutral, cued, uncued, neutral, cued, uncued, stringsAsFactors = FALSE)
#------------------------------------------------------------------------------------------------------------------------
colnames(E1precs) <- c("n1", "c1", "u1", "n2", "c2", "u2")
colnames(E2precs) <- c("n1", "c1", "u1", "n2", "c2", "u2")

count <- 1
for(i in sublist1){
  temp <- dplyr::filter(data, subID == i & experiment==1)
  tempneut  <- dplyr::filter(temp, cond == "neutral")
  tempcued1 <- dplyr::filter(temp, cond == "cued1")
  tempcued2 <- dplyr::filter(temp, cond == "cued2")
  E1precs[count,1] <- mean(abs(tempneut[,14]))
  E1precs[count,2] <- mean(abs(tempcued1[,14]))
  E1precs[count,3] <- mean(abs(tempcued2[,14]))
  E1precs[count,4] <- mean(abs(tempneut[,15]))
  E1precs[count,5] <- mean(abs(tempcued2[,15]))
  E1precs[count,6] <- mean(abs(tempcued1[,15]))
  count <- count + 1
} # put Experiment 1 accuracy data into matrix of Participant x Condition
count <- 1
for(i in sublist2){
  temp <- dplyr::filter(data, subID == i & experiment == 2)
  tempneut  <- dplyr::filter(temp, cond == "neutral")
  tempcued1 <- dplyr::filter(temp, cond == "cued1")
  tempcued2 <- dplyr::filter(temp, cond == "cued2")
  E2precs[count,1] <- mean(abs(tempneut[,14]))
  E2precs[count,2] <- mean(abs(tempcued1[,14]))
  E2precs[count,3] <- mean(abs(tempcued2[,14]))
  E2precs[count,4] <- mean(abs(tempneut[,15]))
  E2precs[count,5] <- mean(abs(tempcued2[,15]))
  E2precs[count,6] <- mean(abs(tempcued1[,15]))
  count <- count + 1
} # put Experiment 2 accuracy data into matrix of Participant x Condition

E1accuracy <- summarise_each(E1precs, funs(mean)); E1accuracy <- as.data.frame(t(E1accuracy)); colnames(E1accuracy)[1] <- "accuracy"
E2accuracy <- summarise_each(E2precs, funs(mean)); E2accuracy <- as.data.frame(t(E2accuracy)); colnames(E2accuracy)[1] <- "accuracy"
se <- function(x) sd(x)/sqrt(length(x))

#this makes data frames for each experiment for easy plotting
E1accuracy$condition  <- seq(1, dim(E1accuracy)[1], by=1);
E1accuracy$probe      <- seq(1, dim(E1accuracy)[1], by=1);
E1accuracy$experiment <- seq(1, dim(E1accuracy)[1], by=1);
E1accuracy$sem        <- t(summarise_each(E1precs, funs(se)));

E2accuracy$condition  <- seq(1, dim(E2accuracy)[1], by=1);
E2accuracy$probe      <- seq(1, dim(E2accuracy)[1], by=1);
E2accuracy$experiment <- seq(1, dim(E2accuracy)[1], by=1);
E2accuracy$sem        <- t(summarise_each(E2precs, funs(se)))

E1accuracy[c(1,2,3),]$probe <- "first"; E1accuracy[c(4,5,6),]$probe <- "second"
E1accuracy[c(1,4),]$condition <- "neutral";
E1accuracy[c(2,5),]$condition <- "cued";
E1accuracy[c(3,6),]$condition <- "uncued";
E1accuracy$experiment <- 1
E1accuracy$condition <- as.factor(E1accuracy$condition); E1accuracy$probe <- as.factor(E1accuracy$probe)

E2accuracy[c(1,2,3),]$probe <- "first"; E2accuracy[c(4,5,6),]$probe <- "second"
E2accuracy[c(1,4),]$condition <- "neutral";
E2accuracy[c(2,5),]$condition <- "cued";
E2accuracy[c(3,6),]$condition <- "uncued";
E2accuracy$experiment <- 2
E2accuracy$condition <- as.factor(E2accuracy$condition); E2accuracy$probe <- as.factor(E2accuracy$probe)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#these plots should be accurate (SRC: 7/11/16) -- check!
ggplot(E1accuracy, aes(x=probe, y=accuracy)) + #generate ggplot, plotting probe and precision
  ggtitle("Experiment 1 - accuracy") +
  geom_bar(aes(fill=condition), stat="identity", color="black", position=position_dodge()) + #create the bars and fill by cue condition
  scale_fill_manual(values = c("neutral" = "#ffffbf","cued" = "#fc8d59","uncued" = "#91bfdb")) + #manually set the colours for the bars
  xlab("Probe") + ylab("Accuracy (1/SD)") +  #set axis labels
  #plot the error bars
  geom_errorbar(aes(x=probe, ymin=accuracy-sem, ymax = accuracy+sem, fill = condition), stat="identity", position=position_dodge(0.9), width = 0.1) +
  scale_colour_manual(values = c("cued" = "#000000" , "neutral" = "#000000", "uncued" = "#000000")) + #set colour of errs manually
  #geom_point(data=plotprecs, aes(x=probe, y=precision, fill = cue), position = position_dodge(0.9)) + #plot individual people
  theme(panel.background = element_rect(fill = "white"), panel.grid.minor = element_blank(), axis.text = element_text(size = 14)) + ylim(0,1.6)

ggplot(E2accuracy, aes(x=probe, y=accuracy)) + #generate ggplot, plotting probe and precision
  ggtitle("Experiment 2 - accuracy") +
  geom_bar(aes(fill=condition), stat="identity", color="black", position=position_dodge()) + #create the bars and fill by cue condition
  scale_fill_manual(values = c("neutral" = "#ffffbf","cued" = "#fc8d59","uncued" = "#91bfdb")) + #manually set the colours for the bars
  xlab("Probe") + ylab("Accuracy (1/SD)") +  #set axis labels
  #plot the error bars
  geom_errorbar(aes(x=probe, ymin=accuracy-sem, ymax = accuracy+sem, fill = condition), stat="identity", position=position_dodge(0.9), width = 0.1) +
  scale_colour_manual(values = c("cued" = "#000000" , "neutral" = "#000000", "uncued" = "#000000")) + #set colour of errs manually
  #geom_point(data=plotprecs, aes(x=probe, y=precision, fill = cue), position = position_dodge(0.9)) + #plot individual people
  theme(panel.background = element_rect(fill = "white"), panel.grid.minor = element_blank(), axis.text = element_text(size = 14)) + ylim(0,1.6)
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# here we're going to put E1accuracy and E2accuracy into long format data frames
id <- seq(1, 6*length(sublist1), by = 1)
E1acc <- data.frame(id, stringsAsFactors = FALSE); colnames(E1acc)[1] <- "id"
E2acc <- data.frame(id, stringsAsFactors = FALSE); colnames(E2acc)[1] <- "id"

E1acc$accuracy   <- seq(1, dim(E1precs)[2], by = 1)
E1acc$probe      <- seq(1, dim(E1precs)[2], by = 1)
E1acc$cue        <- seq(1, dim(E1precs)[2], by = 1)
E1acc$experiment <- 1

E2acc$accuracy   <- seq(1, dim(E2precs)[2], by = 1)
E2acc$probe      <- seq(1, dim(E2precs)[2], by = 1)
E2acc$cue        <- seq(1, dim(E2precs)[2], by = 1)
E2acc$experiment <- 2

count <-  1; sub <- 1 #we have to reset subject numbers starting from 1 so it doesn't interfere with anova function later
for(i in sublist1){
  temp <- dplyr::filter(data, experiment == 1 & subID == i)
  tempneut  <- dplyr::filter(temp, cond == "neutral")
  tempcued1 <- dplyr::filter(temp, cond == "cued1")
  tempcued2 <- dplyr::filter(temp, cond == "cued2")
  E1acc[count,2] <- 1/sd(tempneut[,14]) ; E1acc[count,1] <- sub; count <- count + 1
  E1acc[count,2] <- 1/sd(tempcued1[,14]); E1acc[count,1] <- sub; count <- count + 1
  E1acc[count,2] <- 1/sd(tempcued2[,14]); E1acc[count,1] <- sub; count <- count + 1
  E1acc[count,2] <- 1/sd(tempneut[,15]) ; E1acc[count,1] <- sub; count <- count + 1
  E1acc[count,2] <- 1/sd(tempcued2[,15]); E1acc[count,1] <- sub; count <- count + 1
  E1acc[count,2] <- 1/sd(tempcued1[,15]); E1acc[count,1] <- sub; count <- count + 1
  sub <- sub + 1
} #put E1precs into long format

count <-  1; # DO NOT RESET SUB, IT KEEPS PARTICIPANTS FROM THE EXPERIMENTS WITH SEPARATE ID NUMBERS!!!
for(i in sublist2){
  temp <- dplyr::filter(data, experiment == 2 & subID == i)
  tempneut  <- dplyr::filter(temp, cond == "neutral")
  tempcued1 <- dplyr::filter(temp, cond == "cued1")
  tempcued2 <- dplyr::filter(temp, cond == "cued2")
  E2acc[count,2] <- 1/sd(tempneut[,14]) ; E2acc[count,1] <- sub; count <- count + 1
  E2acc[count,2] <- 1/sd(tempcued1[,14]); E2acc[count,1] <- sub; count <- count + 1
  E2acc[count,2] <- 1/sd(tempcued2[,14]); E2acc[count,1] <- sub; count <- count + 1
  E2acc[count,2] <- 1/sd(tempneut[,15]) ; E2acc[count,1] <- sub; count <- count + 1
  E2acc[count,2] <- 1/sd(tempcued2[,15]); E2acc[count,1] <- sub; count <- count + 1
  E2acc[count,2] <- 1/sd(tempcued1[,15]); E2acc[count,1] <- sub; count <- count + 1
  sub <- sub + 1
} #put E2precs into long format

E1acc$probe[E1acc$probe %in% c(1,2,3)] <- "first"; E1acc$probe[E1acc$probe %in% c(4,5,6)] <- "second"
E2acc$probe[E2acc$probe %in% c(1,2,3)] <- "first"; E2acc$probe[E2acc$probe %in% c(4,5,6)] <- "second"

E1acc$cue[E1acc$cue %in% c(1,4)] <- "neutral"; E2acc$cue[E2acc$cue %in% c(1,4)] <- "neutral"
E1acc$cue[E1acc$cue %in% c(2,5)] <- "cued"   ; E2acc$cue[E2acc$cue %in% c(2,5)] <- "cued"   
E1acc$cue[E1acc$cue %in% c(3,6)] <- "uncued" ; E2acc$cue[E2acc$cue %in% c(3,6)] <- "uncued"
E1acc$cue <- as.factor(E1acc$cue); E1acc$probe <- as.factor(E1acc$probe); E1acc$id <- as.factor(E1acc$id); E1acc$experiment <- as.factor(E1acc$experiment)

aov_E1 <- aov_ez('id', 'accuracy', data = E1acc, within = c("probe", "cue"))
nice(aov_E1, es = 'pes')
bf_E1 <- anovaBF(accuracy ~id + cue*probe, data = E1acc, whichRandom = 'id', iterations = 200000, progress = TRUE, whichModels = "withmain")
bf_E1[1] #main effect of probe
bf_E1[2] #main effect of cue
lmbf_E1_full  <- lmBF(accuracy ~ id + cue + probe + cue:probe #+
                                 #id:cue + id:probe
                                 ,data = E1acc, whichRandom = 'id', iterations = 200000)
lmbf_E1_mains <- lmBF(accuracy ~ id + cue + probe #+
                                 #id:cue + id:probe
                                 ,data = E1acc, whichRandom = 'id', iterations = 200000)
lmbf_E1_full/lmbf_E1_mains #evidence for model with cue:probe interaction, vs. just main effects
#when models dont include id interactions with fixed effects,anovaBF and lmBF functions are equivalent doing model comparisons
bf_E1[4]/bf_E1[3] # evidence for model including cue:probe interaction vs. just main effects


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#comparing lmBF and anovaBF
anovabf_E1 <- anovaBF(accuracy ~  id + cue + probe + cue:probe + #fixed effects
                                  id:cue + id:probe, #random effects
                                  data = E1acc, whichRandom = 'id', iterations = 200000) #all models using anovaBF, full model = model 4
lmbf_E1_full   <- lmBF(accuracy ~ id + cue + probe + cue:probe + #fixed effects
                                  id:cue + id:probe, #random effects
                                  data = E1acc, whichRandom = 'id', iterations = 200000) #full model using lmBF
lmbf_E1_id <- lmBF(accuracy ~ id, data = E1acc, whichRandom = 'id', iterations = 200000) #model of just id (between subjects variance, random factor)
lmbf_E1_ModelProbe <- lmBF(accuracy ~ id + probe, data = E1acc, whichRandom = 'id', iterations = 200000)
lmbf_E1_ModelProbeCue <- lmBF(accuracy ~ id + probe + cue #+
                                         #id:probe + id:cue
                                        ,data = E1acc, whichRandom = 'id', iterations = 200000)
anovabf_E1[4] # full model, compared to model of between-subjects variance (as per function denominator!)
lmbf_E1_full/lmbf_E1_id # full model compared to between subjects variance, these two change by a multiple of 5, but they're huge huge huge anyways

anovabf_E1[1] # main effect of probe (i.e. probe+id vs id, as per anovaBF function denominator for BF calculation)
lmbf_E1_ModelProbe/lmbf_E1_id # model of probe + id vs just id (i.e. main effect of probe w/o id)

anovabf_E1[4]/anovabf_E1[3]
lmbf_E1_full/lmbf_E1_ModelProbeCue

#these bayesfactors for computing the same main effect (probe) are within the variance limits of each other
#as long as you do the same comparison (i.e. probe + id vs. id) like the anovaBF function does, suggesting that they're equivalent
# when the model comparisons are appropriate
# suggests that it doesnt matter which we use, but stick to anovaBF as it compares to between-subjects variance as a default, so the bayesfactors are more useful, i think
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


E1_NeutralUncued <- dplyr::filter(E1acc, cue == "neutral" | cue == "uncued")
E1_NeutralUncued$id <- as.factor(E1_NeutralUncued$id)
E1_NeutralUncued$cue <- as.character(E1_NeutralUncued$cue);     E1_NeutralUncued$cue <- as.factor(E1_NeutralUncued$cue)
E1_NeutralUncued$probe <- as.character(E1_NeutralUncued$probe); E1_NeutralUncued$probe <- as.factor(E1_NeutralUncued$probe)
aov_E1_NU <- aov_ez('id', 'accuracy', data = E1_NeutralUncued, within = c("probe", "cue"))
nice(aov_E1_NU, es = 'pes')
bf_E1_NU <- anovaBF(accuracy ~id + cue*probe, data = E1_NeutralUncued, whichRandom = 'id', iterations = 200000, progress = TRUE, whichModels = "withmain")

E1_NeutralCued <- dplyr::filter(E1acc, cue == "neutral" | cue == "cued")
E1_NeutralCued$id <- as.factor(E1_NeutralCued$id)
E1_NeutralCued$cue <- as.character(E1_NeutralCued$cue);     E1_NeutralCued$cue <- as.factor(E1_NeutralCued$cue)
E1_NeutralCued$probe <- as.character(E1_NeutralCued$probe); E1_NeutralCued$probe <- as.factor(E1_NeutralCued$probe)
aov_E1_NC <- aov_ez('id', 'accuracy', data = E1_NeutralCued, within = c("probe", "cue"))
nice(aov_E1_NC, es = 'pes')
bf_E1_NC <- anovaBF(accuracy ~id + cue*probe, data = E1_NeutralCued, whichRandom = 'id', iterations = 200000, progress = TRUE, whichModels = "withmain")
# using bf_E1_NC[7]/bf_E1_NC[4] as the probe:cue interaction, get BF 1.09. probe:cue + id against id (bf_E1_NC[3]) gives BF = 0.593

bf_tE1_NC1 <- ttestBF(x = E1precs$n1, y = E1precs$c1, paired = TRUE); bf_tE1_NC2 <- ttestBF(x = E1precs$n2, y = E1precs$c2, paired = TRUE)
bf_tE1_NU1 <- ttestBF(x = E1precs$n1, y = E1precs$u1, paired = TRUE); bf_tE1_NU2 <- ttestBF(x = E1precs$n2, y = E1precs$u2, paired = TRUE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
E2acc$cue <- as.factor(E2acc$cue); E2acc$probe <- as.factor(E2acc$probe); E2acc$id <- as.factor(E2acc$id); E2acc$experiment <- as.factor(E2acc$experiment)
aov_E2 <- aov_ez('id', 'accuracy', data = E2acc, within = c("probe", "cue"))
nice(aov_E2, es = 'pes')
bf_E2 <- anovaBF(accuracy ~id + cue*probe, data = E2acc, whichRandom = 'id', iterations = 200000, progress = TRUE, whichModels = "withmain")

E2_NeutralUncued <- dplyr::filter(E2acc, cue == "neutral" | cue == "uncued")
E2_NeutralUncued$id <- as.factor(E2_NeutralUncued$id)
E2_NeutralUncued$cue <- as.character(E2_NeutralUncued$cue);     E2_NeutralUncued$cue <- as.factor(E2_NeutralUncued$cue)
E2_NeutralUncued$probe <- as.character(E2_NeutralUncued$probe); E2_NeutralUncued$probe <- as.factor(E2_NeutralUncued$probe)
aov_E2_NU <- aov_ez('id', 'accuracy', data = E2_NeutralUncued, within = c("probe", "cue"))
nice(aov_E2_NU, es = 'pes')
bf_E2_NU <- anovaBF(accuracy ~id + cue*probe, data = E2_NeutralUncued, whichRandom = 'id', iterations = 200000, progress = TRUE, whichModels = "withmain")


E2_NeutralCued <- dplyr::filter(E2acc, cue == "neutral" | cue == "cued")
E2_NeutralCued$id <- as.factor(E2_NeutralCued$id)
E2_NeutralCued$cue <- as.character(E2_NeutralCued$cue);     E2_NeutralCued$cue <- as.factor(E2_NeutralCued$cue)
E2_NeutralCued$probe <- as.character(E2_NeutralCued$probe); E2_NeutralCued$probe <- as.factor(E2_NeutralCued$probe)
aov_E2_NC <- aov_ez('id', 'accuracy', data = E2_NeutralCued, within = c("probe", "cue"))
nice(aov_E2_NC, es = 'pes')
bf_E2_NC <- anovaBF(accuracy ~id + cue*probe, data = E2_NeutralCued, whichRandom = 'id', iterations = 200000, progress = TRUE, whichModels = "withmain")


bf_tE2_NC1 <- ttestBF(x = E2precs$n1, y = E2precs$c1, paired = TRUE); bf_tE2_NC2 <- ttestBF(x = E2precs$n2, y = E2precs$c2, paired = TRUE)
bf_tE2_NU1 <- ttestBF(x = E2precs$n1, y = E2precs$u1, paired = TRUE); bf_tE2_NU2 <- ttestBF(x = E2precs$n2, y = E2precs$u2, paired = TRUE)

# # # # Mixed Effects Models for between-experiment analysis # # # #

combined <- as.data.frame(bind_rows(E1acc, E2acc)) #combines both dataframes together
combined_NeutralCued   <- as.data.frame(bind_rows(E1acc, E2acc))
combined_NeutralUncued <- as.data.frame(bind_rows(E1acc, E2acc))
combined_Neutral       <- as.data.frame(bind_rows(E1acc, E2acc))
combined_NeutralFirst  <- as.data.frame(bind_rows(E1acc, E2acc))
combined_NeutralSecond <- as.data.frame(bind_rows(E1acc, E2acc))
combined_NeutralCuedFirst <- as.data.frame(bind_rows(E1acc,E2acc))
combined_NeutralUncuedFirst <- as.data.frame(bind_rows(E1acc,E2acc))
combined_NeutralCuedSecond <- as.data.frame(bind_rows(E1acc,E2acc))
combined_NeutralUncuedSecond <- as.data.frame(bind_rows(E1acc,E2acc))

combined$id  <- as.factor(combined$id) ; combined$probe      <- as.factor(combined$probe)
combined$cue <- as.factor(combined$cue); combined$experiment <- as.factor(combined$experiment)

mixedEffects    <- aov_ez('id', 'accuracy', data = combined, within = c("probe", "cue"), between = c("experiment"))
nice(mixedEffects, es = 'pes')

#bf_mixedEffects <- anovaBF(accuracy ~ id + cue*probe*experiment, data = combined, whichRandom = "id", iterations = 200000, progress = TRUE, whichModels = "withmain")
# bf_mixedEffects[4]/bf_mixedEffects[3] # <- evidence favouring a model that has both main effects of probe/cue, and probe:cue interaction, versus without the interaction
#bf_mixedEffects[18]/bf_mixedEffects[17]


# mixed effect anovas between cueing conditions
combined_NeutralCued <- dplyr::filter(combined_NeutralCued, cue == "neutral" | cue == "cued")
combined_NeutralCued$id <- as.factor(combined_NeutralCued$id) ; combined_NeutralCued$probe <- as.factor(combined_NeutralCued$probe)
combined_NeutralCued$cue  <- as.factor(combined_NeutralCued$cue) ; combined_NeutralCued$experiment <- as.factor(combined_NeutralCued$experiment)

combined_NeutralUncued <- dplyr::filter(combined_NeutralUncued, cue == "neutral" | cue == "uncued")
combined_NeutralUncued$id <- as.factor(combined_NeutralUncued$id) ; combined_NeutralUncued$probe <- as.factor(combined_NeutralUncued$probe)
combined_NeutralUncued$cue  <- as.factor(combined_NeutralUncued$cue) ; combined_NeutralUncued$experiment <- as.factor(combined_NeutralUncued$experiment)

combined_Neutral <- dplyr::filter(combined_Neutral, cue == "neutral")
combined_Neutral$id <- as.factor(combined_Neutral$id) ; combined_Neutral$probe <- as.factor(combined_Neutral$probe)
combined_Neutral$cue  <- as.factor(combined_Neutral$cue) ; combined_Neutral$experiment <- as.factor(combined_Neutral$experiment)

combined_NeutralFirst <- dplyr::filter(combined_NeutralFirst, cue == 'neutral' & probe == 'first')
combined_NeutralFirst$cue <- as.character(combined_NeutralFirst$cue); combined_NeutralFirst$cue <- as.factor(combined_NeutralFirst$cue)
combined_NeutralFirst$probe <- as.character(combined_NeutralFirst$probe); combined_NeutralFirst$probe <- as.factor(combined_NeutralFirst$probe)
combined_NeutralFirst$experiment <- as.factor(combined_NeutralFirst$experiment); combined_NeutralFirst$id <- as.factor(combined_NeutralFirst$id)

combined_NeutralSecond <- dplyr::filter(combined_NeutralSecond, cue == 'neutral' & probe == 'second')
combined_NeutralSecond$cue <- as.character(combined_NeutralSecond$cue); combined_NeutralSecond$cue <- as.factor(combined_NeutralSecond$cue)
combined_NeutralSecond$probe <- as.character(combined_NeutralSecond$probe); combined_NeutralSecond$probe <- as.factor(combined_NeutralSecond$probe)
combined_NeutralSecond$experiment <- as.factor(combined_NeutralSecond$experiment); combined_NeutralSecond$id <- as.factor(combined_NeutralSecond$id)

combined_NeutralCuedFirst <- dplyr::filter(combined_NeutralCuedFirst, (cue == "neutral" | cue == "cued") & probe == "first")
combined_NeutralCuedFirst$id <- as.factor(combined_NeutralCuedFirst$id) ; combined_NeutralCuedFirst$probe <- as.factor(combined_NeutralCuedFirst$probe)
combined_NeutralCuedFirst$cue  <- as.factor(combined_NeutralCuedFirst$cue) ; combined_NeutralCuedFirst$experiment <- as.factor(combined_NeutralCuedFirst$experiment)

combined_NeutralUncuedFirst <- dplyr::filter(combined_NeutralUncuedFirst, (cue == "neutral" | cue == "uncued") & probe == "first")
combined_NeutralUncuedFirst$id <- as.factor(combined_NeutralUncuedFirst$id) ; combined_NeutralUncuedFirst$probe <- as.factor(combined_NeutralUncuedFirst$probe)
combined_NeutralUncuedFirst$cue  <- as.factor(combined_NeutralUncuedFirst$cue) ; combined_NeutralUncuedFirst$experiment <- as.factor(combined_NeutralUncuedFirst$experiment)

combined_NeutralCuedSecond <- dplyr::filter(combined_NeutralCuedSecond, (cue == "neutral" | cue == "cued") & probe == "second")
combined_NeutralCuedSecond$id <- as.factor(combined_NeutralCuedSecond$id) ; combined_NeutralCuedSecond$probe <- as.factor(combined_NeutralCuedSecond$probe)
combined_NeutralCuedSecond$cue  <- as.factor(combined_NeutralCuedSecond$cue) ; combined_NeutralCuedSecond$experiment <- as.factor(combined_NeutralCuedSecond$experiment)

combined_NeutralUncuedSecond <- dplyr::filter(combined_NeutralUncuedSecond, (cue == "neutral" | cue == "uncued") & probe == "second")
combined_NeutralUncuedSecond$id <- as.factor(combined_NeutralUncuedSecond$id) ; combined_NeutralUncuedSecond$probe <- as.factor(combined_NeutralUncuedSecond$probe)
combined_NeutralUncuedSecond$cue  <- as.factor(combined_NeutralUncuedSecond$cue) ; combined_NeutralUncuedSecond$experiment <- as.factor(combined_NeutralUncuedSecond$experiment)

mixedEffects_Neutral             <- aov_ez('id', 'accuracy', data = combined_Neutral,             within = c("probe")       , between = "experiment")
mixedEffects_NeutralFirst        <- aov_ez('id', 'accuracy', data = combined_NeutralFirst                                   , between = "experiment")
mixedEffects_NeutralSecond       <- aov_ez('id', 'accuracy', data = combined_NeutralSecond                                  , between = "experiment")
mixedEffects_NeutralCued         <- aov_ez('id', 'accuracy', data = combined_NeutralCued,         within = c("probe", "cue"), between = "experiment")
mixedEffects_NeutralUncued       <- aov_ez('id', 'accuracy', data = combined_NeutralUncued,       within = c("probe", "cue"), between = "experiment")
mixedEffects_NeutralCuedFirst    <- aov_ez('id', 'accuracy', data = combined_NeutralCuedFirst,    within = "cue"            , between = "experiment")
mixedEffects_NeutralCuedSecond   <- aov_ez('id', 'accuracy', data = combined_NeutralCuedSecond,   within = "cue"            , between = "experiment")
mixedEffects_NeutralUncuedFirst  <- aov_ez('id', 'accuracy', data = combined_NeutralUncuedFirst,  within = "cue"            , between = "experiment")
mixedEffects_NeutralUncuedSecond <- aov_ez('id', 'accuracy', data = combined_NeutralUncuedSecond, within = "cue"            , between = "experiment")

nice(mixedEffects,                     es = 'pes')
nice(mixedEffects_Neutral,             es = 'pes')
nice(mixedEffects_NeutralFirst,        es = 'pes')
nice(mixedEffects_NeutralSecond,       es = 'pes')
nice(mixedEffects_NeutralCued,         es = 'pes')
nice(mixedEffects_NeutralUncued,       es = 'pes')
nice(mixedEffects_NeutralCuedFirst,    es = 'pes')
nice(mixedEffects_NeutralCuedSecond,   es = 'pes')
nice(mixedEffects_NeutralUncuedFirst,  es = 'pes')
nice(mixedEffects_NeutralUncuedSecond, es = 'pes')