# Analysis of Experiment 3 for the DoubleProbe JEP:HPP paper - SRC
#### install necessary packages -----
library("ggplot2")
library("dplyr")
library("circular")
library("pracma")
library("afex")

### load in data ----
#fileList <- list.files(path= "/Users/user/Desktop/DoubleProbe3/data/csv/", pattern = ".csv") #change path to your data directory
#dataFiles <- list(NULL) #creates blank list
#count <- 1
#for(i in fileList){
#  path <- paste0("/Users/user/Desktop/DoubleProbe3/data/csv/",i) #this generates a filename to access a file
#  dataFiles[[count]] <- read.csv(path, header=TRUE,as.is=TRUE) #this says that the first item in this list the file above, and reads in the csv file you have
#  count <- count + 1 #moves onto the next
#}

### collating data, adding response accuracy and conditions ----
#df <- do.call("rbind",dataFiles) #creates a dataframe by binding all the objects in the list dataFiles into one thing called df
#df <- df[order(df$subID),] #orders the rows in the dataframe by subID (subject ID in this case)
#df <- dplyr::select(df, -angs_3, -angs_4)

#if have the allData.csv file, can ignore above..
path <- "/Users/user/Desktop/Experiments/Nick/DoubleProbe/E3/data/ExpDat/allData_E3.csv" #change your path here!!!
df   <- read.csv(path, header=TRUE, as.is=TRUE)
df <- dplyr::select(df, -X)

df$rdif1 <- seq(1,dim(df)[1],by=1) #make a new column called rdif1
df$rdif2 <- seq(1,dim(df)[1],by=1) #make a new column called rdif2
df$cond  <- seq(1,dim(df)[1],by=1) #make a new column called cond

for(i in seq(18,20,by=1)) df[,i] <- 0

for(i in seq(1,dim(df)[1]))
  if(df[i,12] == 1 & df[i,13]==1) df[i,20] <- "cued1"   #column 12 is cues, column 9 is probe order (pord)
for(i in seq(1,dim(df)[1]))
  if(df[i,12] == 1 & df[i,13]==2) df[i,20] <- "cued2"
for(i in seq(1,dim(df)[1]))
  if(df[i,12] == 0) df[i,20] <- "neutral"


df <- dplyr::mutate(df, rdif1 = as.circular(theta_1 - rad(angs_1), units="radians")) #don't change this!!
df <- dplyr::mutate(df, rdif2 = as.circular(theta_2 - rad(angs_2), units="radians")) #don't change this!!

sublist <- c(2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 22, 24)
neutral1 <- vector("numeric", length=length(sublist))
cued1 <- vector("numeric", length=length(sublist))
uncued1 <- vector("numeric", length=length(sublist))
neutral2 <- vector("numeric", length=length(sublist))
cued2 <- vector("numeric", length=length(sublist))
uncued2 <- vector("numeric", length=length(sublist))
#------------------------------------------------------------------------------------------------------------------------
precs_old    <- data.frame(neutral1, cued1, uncued1, neutral2, cued2, uncued2, stringsAsFactors = FALSE)
precs_new    <- data.frame(neutral1, cued1, uncued1, neutral2, cued2, uncued2, stringsAsFactors = FALSE)
#------------------------------------------------------------------------------------------------------------------------
colnames(precs_old) <- c("n1", "c1", "u1", "n2", "c2", "u2")
colnames(precs_new) <- c("n1", "c1", "u1", "n2", "c2", "u2")

count <- 1
for(i in sublist){
  tempsub   <- dplyr::filter(df, subID == i)
  tempblock <- dplyr::filter(tempsub, btype==0)
  tempneut  <- dplyr::filter(tempblock, cond == "neutral")
  tempcued1 <- dplyr::filter(tempblock, cond == "cued1")
  tempcued2 <- dplyr::filter(tempblock, cond == "cued2")
  
  precs_old[count,1]    <- 1/sd(tempneut[,18]) #rdif 1 = column 18
  precs_old[count,2]    <- 1/sd(tempcued1[,18])
  precs_old[count,3]    <- 1/sd(tempcued2[,18])
  precs_old[count,4]    <- 1/sd(tempneut[,19]) #rdif 2 - column 19
  precs_old[count,5]    <- 1/sd(tempcued2[,19])
  precs_old[count,6]    <- 1/sd(tempcued1[,19])
  count <- count+1
}
count <- 1
for(i in sublist){
  tempblock <- dplyr::filter(df, btype==1)
  tempsub   <- dplyr::filter(tempblock, subID == i)
  tempneut  <- dplyr::filter(tempsub, cond == "neutral")
  tempcued1 <- dplyr::filter(tempsub, cond == "cued1")
  tempcued2 <- dplyr::filter(tempsub, cond == "cued2")
  
  precs_new[count,1]    <- 1/sd(tempneut[,18])
  precs_new[count,2]    <- 1/sd(tempcued1[,18])
  precs_new[count,3]    <- 1/sd(tempcued2[,18])
  precs_new[count,4]    <- 1/sd(tempneut[,19])
  precs_new[count,5]    <- 1/sd(tempcued2[,19])
  precs_new[count,6]    <- 1/sd(tempcued1[,19])
  count <- count+1
}

#------------------------------------------------------------------------------------------------------------------------
prec_old     <- summarise_each(precs_old, funs(mean));  prec_old     <- as.data.frame(t(prec_old))
prec_new     <- summarise_each(precs_new, funs(mean));  prec_new     <- as.data.frame(t(prec_new))
#------------------------------------------------------------------------------------------------------------------------
colnames(prec_old)[1] <- "precision"; 
colnames(prec_new)[1] <- "precision"; 
#------------------------------------------------------------------------------------------------------------------------
# reformat data into new dataframe for plotting 
se <- function(x) sd(x)/sqrt(length(x))

prec_old$condition <- seq(1, dim(prec_old)[1],by=1); prec_old$probe <- seq(1, dim(prec_old)[1], by=1); prec_old$experiment <- 1
prec_old$sem <- t(summarise_each(precs_old, funs(se)))
#------------------------------------------------------------------------------------------------------------------------
prec_new$condition <- seq(1, dim(prec_new)[1],by=1); prec_new$probe <- seq(1, dim(prec_new)[1], by=1); prec_new$experiment <- 2
prec_new$sem <- t(summarise_each(precs_new, funs(se)))
#------------------------------------------------------------------------------------------------------------------------
# cycle through values and label them with conditions
prec_old[c(1,2,3),]$probe     <- "first"; prec_old[c(4,5,6),]$probe    <- "second"
prec_old[c(1,4),]$condition     <- "neutral"; prec_old[c(2,5),]$condition     <- "cued"; prec_old[c(3,6),]$condition     <- "uncued"
prec_old$condition <- as.factor(prec_old$condition); prec_old$probe <- as.factor(prec_old$probe)
#------------------------------------------------------------------------------------------------------------------------
prec_new[c(1,2,3),]$probe     <- "first"; prec_new[c(4,5,6),]$probe    <- "second"
prec_new[c(1,4),]$condition     <- "neutral"; prec_new[c(2,5),]$condition     <- "cued"; prec_new[c(3,6),]$condition     <- "uncued"
prec_new$condition <- as.factor(prec_new$condition); prec_new$probe <- as.factor(prec_new$probe)
#------------------------------------------------------------------------------------------------------------------------
# ggplots of the first block type

setEPS() #comment/uncomment these lines to save high-res images of ggplots that are generated
postscript("/Users/user/Desktop/DoubleProbe3/NoOrder_accuracy_n=20.eps")
ggplot(prec_old, aes(x=probe, y=precision)) + #generate ggplot, plotting probe and precision
  ggtitle("spatial cue - accuracy, n=20") +
  geom_bar(aes(fill=condition), stat="identity", color="black", position=position_dodge()) + #create the bars and fill by cue condition
  scale_fill_manual(values = c("neutral" = "#ffffbf","cued" = "#fc8d59","uncued" = "#91bfdb")) + #manually set the colours for the bars
  xlab("Probe") + ylab("Accuraccy (1/SD)") +  #set axis labels
  #plot the error bars
  geom_errorbar(aes(x=probe, ymin=precision-sem, ymax = precision+sem, fill = condition), stat="identity", position=position_dodge(0.9), width = 0.1) +
  scale_colour_manual(values = c("cued" = "#000000" , "neutral" = "#000000", "uncued" = "#000000")) + #set colour of errs manually
  #geom_point(data=plotprecs, aes(x=probe, y=precision, fill = cue), position = position_dodge(0.9)) + #plot individual people
  theme(panel.background = element_rect(fill = "white"), panel.grid.minor = element_blank(), axis.text = element_text(size = 14)) + ylim(0,2.0)
dev.off()

#------------------------------------------------------------------------------------------------------------------------
# ggplots of the second block type

setEPS()
postscript("/Users/user/Desktop/DoubleProbe3/OrderInfo_accuracy_n=20.eps")
ggplot(prec_new, aes(x=probe, y=precision)) + #generate ggplot, plotting probe and precision
  ggtitle("spatio-temporal cue - accuracy, n=20") +
  geom_bar(aes(fill=condition), stat="identity", color="black", position=position_dodge()) + #create the bars and fill by cue condition
  scale_fill_manual(values = c("neutral" = "#ffffbf","cued" = "#fc8d59","uncued" = "#91bfdb")) + #manually set the colours for the bars
  xlab("Probe") + ylab("Accuraccy (1/SD)") +  #set axis labels
  #plot the error bars
  geom_errorbar(aes(x=probe, ymin=precision-sem, ymax = precision+sem, fill = condition), stat="identity", position=position_dodge(0.9), width = 0.1) +
  scale_colour_manual(values = c("cued" = "#000000" , "neutral" = "#000000", "uncued" = "#000000")) + #set colour of errs manually
  #geom_point(data=plotprecs, aes(x=probe, y=precision, fill = cue), position = position_dodge(0.9)) + #plot individual people
  theme(panel.background = element_rect(fill = "white"), panel.grid.minor = element_blank(), axis.text = element_text(size = 14)) + ylim(0,2.0)
dev.off()
#------------------------------------------------------------------------------------------------------------------------


##### not for plotting, but actually for ANOVAs
#plotsubs <- c(2,3,4,5,6,7,8,9,11,12,13,14,15,16,17)
plotsubs <- c(1:20) #just using plotsubs now instead of sublist as we only have 20 in the dframes being used, plus only using 20 subs
id <- seq(1,6*length(plotsubs),by=1)
plotPrecsOld <- data.frame(id, stringsAsFactors = FALSE)
plotPrecsOld$Precision  <- seq(1, dim(precs_old)[2], by=1)
plotPrecsOld$Cue        <- seq(1, dim(precs_old)[2], by=1) 
plotPrecsOld$Probe      <- seq(1, dim(precs_old)[2], by=1)
plotPrecsOld$Experiment <- seq(1, dim(precs_old)[2], by=1)

count <- 1
for(i in plotsubs){
  plotPrecsOld[count,2] <- precs_old[i,1]; plotPrecsOld[count,1] <- i; count <- count + 1
  plotPrecsOld[count,2] <- precs_old[i,2]; plotPrecsOld[count,1] <- i; count <- count + 1
  plotPrecsOld[count,2] <- precs_old[i,3]; plotPrecsOld[count,1] <- i; count <- count + 1
  plotPrecsOld[count,2] <- precs_old[i,4]; plotPrecsOld[count,1] <- i; count <- count + 1
  plotPrecsOld[count,2] <- precs_old[i,5]; plotPrecsOld[count,1] <- i; count <- count + 1
  plotPrecsOld[count,2] <- precs_old[i,6]; plotPrecsOld[count,1] <- i; count <- count + 1
} #put precs_old into long format -> plotPrecsOld is made

plotPrecsOld$Experiment <- 1
plotPrecsOld$Probe[plotPrecsOld$Probe == 1 | plotPrecsOld$Probe == 2 | plotPrecsOld$Probe == 3] <- "first"
plotPrecsOld$Probe[plotPrecsOld$Probe == 4 | plotPrecsOld$Probe == 5 | plotPrecsOld$Probe == 6] <- "second"
plotPrecsOld$Cue[plotPrecsOld$Cue == 1 | plotPrecsOld$Cue == 4] <- "neutral"
plotPrecsOld$Cue[plotPrecsOld$Cue == 2 | plotPrecsOld$Cue == 5] <- "cued"
plotPrecsOld$Cue[plotPrecsOld$Cue == 3 | plotPrecsOld$Cue == 6] <- "uncued"


#------------------------------------------------------------------------------------------------------------------------
plotPrecsNew <- data.frame(id, stringsAsFactors = FALSE)
plotPrecsNew$Precision  <- seq(1, dim(precs_new)[2], by=1)
plotPrecsNew$Cue        <- seq(1, dim(precs_new)[2], by=1) 
plotPrecsNew$Probe      <- seq(1, dim(precs_new)[2], by=1)
plotPrecsNew$Experiment <- seq(1, dim(precs_new)[2], by=1)

count <- 1
for(i in plotsubs){
  plotPrecsNew[count,2] <- precs_new[i,1]; plotPrecsNew[count,1] <- i; count <- count + 1
  plotPrecsNew[count,2] <- precs_new[i,2]; plotPrecsNew[count,1] <- i; count <- count + 1
  plotPrecsNew[count,2] <- precs_new[i,3]; plotPrecsNew[count,1] <- i; count <- count + 1
  plotPrecsNew[count,2] <- precs_new[i,4]; plotPrecsNew[count,1] <- i; count <- count + 1
  plotPrecsNew[count,2] <- precs_new[i,5]; plotPrecsNew[count,1] <- i; count <- count + 1
  plotPrecsNew[count,2] <- precs_new[i,6]; plotPrecsNew[count,1] <- i; count <- count + 1
} #put precs_new into long format -> plotPrecsNew is made

plotPrecsNew$Experiment <- 2
plotPrecsNew$Probe[plotPrecsNew$Probe == 1 | plotPrecsNew$Probe == 2 | plotPrecsNew$Probe == 3] <- "first"
plotPrecsNew$Probe[plotPrecsNew$Probe == 4 | plotPrecsNew$Probe == 5 | plotPrecsNew$Probe == 6] <- "second"
plotPrecsNew$Cue[plotPrecsNew$Cue == 1 | plotPrecsNew$Cue == 4] <- "neutral"
plotPrecsNew$Cue[plotPrecsNew$Cue == 2 | plotPrecsNew$Cue == 5] <- "cued"
plotPrecsNew$Cue[plotPrecsNew$Cue == 3 | plotPrecsNew$Cue == 6] <- "uncued"

plotprecs <- as.data.frame(bind_rows(plotPrecsOld, plotPrecsNew))
plotprecs$id <- as.factor(plotprecs$id)      ; plotprecs$Cue <- as.factor(plotprecs$Cue);
plotprecs$Probe <- as.factor(plotprecs$Probe); plotprecs$Experiment <- as.factor(plotprecs$Experiment)

plotPrecsOld$id <- as.factor(plotPrecsOld$id); plotPrecsOld$Cue <- as.factor(plotPrecsOld$Cue);
plotPrecsOld$Probe <- as.factor(plotPrecsOld$Probe); plotPrecsOld$Experiment <- as.factor(plotPrecsOld$Experiment);

plotPrecsNew$id <- as.factor(plotPrecsNew$id); plotPrecsNew$Cue <- as.factor(plotPrecsNew$Cue);
plotPrecsNew$Probe <- as.factor(plotPrecsNew$Probe); plotPrecsNew$Experiment <- as.factor(plotPrecsNew$Experiment)

cohens_d <- function(x,y){
  m1 <- base::mean(x); l1 = length(x); variance1 <- stats::var(x); 
  m2 <- base::mean(y); l2 = length(y); variance2 <- stats::var(y)
  sigma <- (sqrt((l1-1)*variance1 + (l2-1)*variance2)/(l1 + l2 - 2))
  d <- (m1-m2)/sigma
  return(d)
}

aov_old <- aov_ez('id', 'Precision', plotPrecsOld, within = c("Cue", "Probe"))              ; nice(aov_old, es='pes')
aov_new <- aov_ez('id', 'Precision', plotPrecsNew, within = c("Cue", "Probe"))              ; nice(aov_new, es='pes')
aov_all <- aov_ez('id', 'Precision', plotprecs   , within = c("Cue", "Probe", "Experiment")); nice(aov_all, es='pes')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

plotPrecsOld_NeutralCued   <- dplyr::filter(plotPrecsOld, Cue != "uncued")
plotPrecsOld_NeutralCued$Cue <- as.character(plotPrecsOld_NeutralCued$Cue); plotPrecsOld_NeutralCued$Cue <- as.factor(plotPrecsOld_NeutralCued$Cue)

plotPrecsOld_NeutralUncued <- dplyr::filter(plotPrecsOld, Cue != "cued")
plotPrecsOld_NeutralUncued$Cue <- as.character(plotPrecsOld_NeutralUncued$Cue); plotPrecsOld_NeutralUncued$Cue <- as.factor(plotPrecsOld_NeutralUncued$Cue)

plotPrecsNew_NeutralCued   <- dplyr::filter(plotPrecsNew, Cue != "uncued")
plotPrecsNew_NeutralCued$Cue <- as.character(plotPrecsNew_NeutralCued$Cue); plotPrecsNew_NeutralCued$Cue <- as.factor(plotPrecsNew_NeutralCued$Cue)

plotPrecsNew_NeutralUncued <- dplyr::filter(plotPrecsNew, Cue != "cued")
plotPrecsNew_NeutralUncued$Cue <- as.character(plotPrecsNew_NeutralUncued$Cue); plotPrecsNew_NeutralUncued$Cue <- as.factor(plotPrecsNew_NeutralUncued$Cue)


plotPrecs_NeutralCued   <- as.data.frame(bind_rows(plotPrecsOld_NeutralCued,plotPrecsNew_NeutralCued))
plotPrecs_NeutralUncued <- as.data.frame(bind_rows(plotPrecsOld_NeutralUncued,plotPrecsNew_NeutralUncued))
plotPrecs_NeutralCued$Experiment <- as.factor(plotPrecs_NeutralCued$Experiment)
plotPrecs_NeutralUncued$Experiment <- as.factor(plotPrecs_NeutralUncued$Experiment)


aov_neutralcued     <- aov_ez('id', 'Precision', plotPrecs_NeutralCued, within = c("Cue","Probe", "Experiment")); nice(aov_neutralcued,     es = 'pes')
aov_neutralcued_old <- aov_ez('id', 'Precision', plotPrecsOld_NeutralCued, within = c("Cue", "Probe"))          ; nice(aov_neutralcued_old, es = 'pes')
aov_neutralcued_new <- aov_ez('id', 'Precision', plotPrecsNew_NeutralCued, within = c("Cue", "Probe"))          ; nice(aov_neutralcued_new, es = 'pes')

n1c1_old <- t.test(precs_old$n1, precs_old$c1, paired=TRUE); n1c1old_effsize <- cohens_d(precs_old$n1, precs_old$c1)
n2c2_old <- t.test(precs_old$n2, precs_old$c2, paired=TRUE); n2c2old_effsize <- cohens_d(precs_old$n2, precs_old$c2)
n1u1_old <- t.test(precs_old$n1, precs_old$u1, paired=TRUE); n1u1old_effsize <- cohens_d(precs_old$n1, precs_old$u1)
n2u2_old <- t.test(precs_old$n2, precs_old$u2, paired=TRUE); n2u2old_effsize <- cohens_d(precs_old$n2, precs_old$u2)

aov_neutraluncued     <- aov_ez('id', 'Precision', plotPrecs_NeutralUncued, within = c("Cue","Probe", "Experiment")); nice(aov_neutraluncued, es = 'pes')
bf_neutraluncued      <- anovaBF(Precision ~ id + Cue*Probe*Experiment, data = plotPrecs_NeutralUncued, whichRandom = 'id', iterations = 200000, progress = TRUE, whichModels = "withmain")

aov_neutraluncued_old <- aov_ez('id', 'Precision', plotPrecsOld_NeutralUncued, within = c("Cue", "Probe")); nice(aov_neutraluncued_old, es = 'pes')
aov_neutraluncued_new <- aov_ez('id', 'Precision', plotPrecsNew_NeutralUncued, within = c("Cue", "Probe")); nice(aov_neutraluncued_new, es = 'pes')

n1c1_new <- t.test(precs_new$n1, precs_new$c1, paired=TRUE); n1c1_new_effsize <- cohens_d(precs_new$n1, precs_new$c1)
n2c2_new <- t.test(precs_new$n2, precs_new$c2, paired=TRUE); n2c2_new_effsize <- cohens_d(precs_new$n2, precs_new$c2)
n1u1_new <- t.test(precs_new$n1, precs_new$u1, paired=TRUE); n1u1_new_effsize <- cohens_d(precs_new$n1, precs_new$u1)
n2u2_new <- t.test(precs_new$n2, precs_new$u2, paired=TRUE); n2u2_new_effsize <- cohens_d(precs_new$n2, precs_new$u2)

plotprecs_neutraluncued_second <- dplyr::filter(plotprecs, Probe == "second" & (Cue == "neutral" | Cue == "uncued"))
plotprecs_neutraluncued_second$Cue <- as.character(plotprecs_neutraluncued_second$Cue); plotprecs_neutraluncued_second$Cue <- as.factor(plotprecs_neutraluncued_second$Cue)
plotprecs_neutraluncued_second$Probe <- as.character(plotprecs_neutraluncued_second$Probe); plotprecs_neutraluncued_second$Probe <- as.factor(plotprecs_neutraluncued_second$Probe)

aov_neutraluncued_second <- aov_ez('id', 'Precision', plotprecs_neutraluncued_second, within = c("Cue", "Experiment"))
bf_neutraluncued_second <- anovaBF(Precision ~ id + Cue*Experiment,data = plotprecs_neutraluncued_second, whichRandom = "id", iterations = 200000, progress = TRUE, whichModels = "all")

plotprecs_neutraluncued_first <- dplyr::filter(plotprecs, Probe == "first" & (Cue == "neutral" | Cue == "uncued"))
plotprecs_neutraluncued_first$Cue <- as.character(plotprecs_neutraluncued_first$Cue); plotprecs_neutraluncued_first$Cue <- as.factor(plotprecs_neutraluncued_first$Cue)
plotprecs_neutraluncued_first$Probe <- as.character(plotprecs_neutraluncued_first$Probe); plotprecs_neutraluncued_first$Probe <- as.factor(plotprecs_neutraluncued_first$Probe)

aov_neutraluncued_first <- aov_ez('id', 'Precision', plotprecs_neutraluncued_first, within = c("Cue", "Experiment"))
bf_neutraluncued_first  <- anovaBF(Precision ~ id + Cue*Experiment,data = plotprecs_neutraluncued_first, whichRandom = "id", iterations = 200000, progress = TRUE, whichModels = "all")

#------------------------------------------------------------------------------------------------------------------------
# Bayesian statistics
library("BayesFactor")

bf_all  <- anovaBF(Precision ~ id + Cue*Probe*Experiment, data = plotprecs,  whichRandom = "id", iterations = 200000, progress=TRUE, whichModels = "withmain")
bf_old  <- anovaBF(Precision ~ id + Cue*Probe,            data = plotPrecsOld, whichRandom = "id", iterations = 200000, progress=TRUE, whichModels = "withmain")
bf_new  <- anovaBF(Precision ~ id + Cue*Probe,            data = plotPrecsNew, whichRandom = "id", iterations = 200000, progress=TRUE, whichModels = "withmain")

bf_all[3]/bf_all[8] #bayes factor for argument against an Experiment effect
#i.e. yields BF of ~5.7, strong evidence that there is no effect of Experiment on Precision (performance no significantly better or worse between block types)
# a model without an effect of experiment fares 6.58 times better

bf_neutralcued_old     <- anovaBF(Precision ~ id + Cue*Probe,        data = plotPrecsOld_NeutralCued,   whichRandom="id", iterations=200000, progress=TRUE, whichModels = "withmain")
bf_neutraluncued_old   <- anovaBF(Precision ~ id + Cue*Probe,       data = plotPrecsOld_NeutralUncued, whichRandom="id", iterations=200000, progress=TRUE, whichModels = "withmain")
bf_neutralcued_new     <- anovaBF(Precision ~ id + Cue*Probe, data = plotPrecsNew_NeutralCued, whichRandom="id", iterations=200000, progress=TRUE, whichModels = "withmain")
bf_neutraluncued_new   <- anovaBF(Precision ~ id + Cue*Probe, data = plotPrecsNew_NeutralUncued, whichRandom="id", iterations=200000, progress=TRUE, whichModels = "withmain")


lmbf_old_full <- lmBF(Precision ~ id + Cue + Probe + Cue:Probe, data = plotPrecsOld, whichRandom= 'id', iterations = 200000)
lmbf_old_main <- lmBF(Precision ~ id + Cue + Probe,             data = plotPrecsOld, whichRandom= 'id', iterations = 200000)
lmbf_old_full/lmbf_old_main #evidence for interaction of cue and probe factors in Experiment 3 blocktype = order not cued

lmbf_old_NC_full <- lmBF(Precision ~ id + Cue + Probe + Cue:Probe,   data = plotPrecsOld_NeutralCued, whichRandom = 'id', iterations  = 200000)
lmbf_old_NC_main <- lmBF(Precision ~ id + Cue + Probe,               data = plotPrecsOld_NeutralCued, whichRandom = 'id', iterations  = 200000)
lmbf_old_NC_full/lmbf_old_NC_main # evidence for cue:probe interaction when comparing cue and probe in only neutral and cued trials, blocktype = order not cued

lmbf_old_NU_full <- lmBF(Precision ~ id + Cue + Probe + Cue:Probe,  data = plotPrecsOld_NeutralUncued, whichRandom = 'id', iterations  = 200000)
lmbf_old_NU_main <- lmBF(Precision ~ id + Cue + Probe,              data = plotPrecsOld_NeutralUncued, whichRandom = 'id', iterations  = 200000)
lmbf_old_NU_full/lmbf_old_NU_main # evidence for cue:probe interaction when comparing cue and probe in only neutral and cued trials, blocktype = order not cued


lmbf_new_full <- lmBF(Precision ~ id + Cue + Probe + Cue:Probe, data = plotPrecsNew, whichRandom= 'id', iterations = 200000)
lmbf_new_main <- lmBF(Precision ~ id + Cue + Probe,             data = plotPrecsNew, whichRandom= 'id', iterations = 200000)
lmbf_new_full/lmbf_new_main #evidence for interaction of cue and probe factors in Experiment 3 blocktype = order not cued

lmbf_new_NC_full <- lmBF(Precision ~ id + Cue + Probe + Cue:Probe,   data = plotPrecsNew_NeutralCued, whichRandom = 'id', iterations  = 200000)
lmbf_new_NC_main <- lmBF(Precision ~ id + Cue + Probe,               data = plotPrecsNew_NeutralCued, whichRandom = 'id', iterations  = 200000)
lmbf_new_NC_full/lmbf_new_NC_main # evidence for cue:probe interaction when comparing cue and probe in only neutral and cued trials, blocktype = order not cued

lmbf_new_NU_full <- lmBF(Precision ~ id + Cue + Probe + Cue:Probe,  data = plotPrecsNew_NeutralUncued, whichRandom = 'id', iterations  = 200000)
lmbf_new_NU_main <- lmBF(Precision ~ id + Cue + Probe,              data = plotPrecsNew_NeutralUncued, whichRandom = 'id', iterations  = 200000)
lmbf_new_NU_full/lmbf_new_NU_main # evidence for cue:probe interaction when comparing cue and probe in only neutral and cued trials, blocktype = order not cued


bf_n1c1_old <- ttestBF(x=precs_old$c1, y=precs_old$n1, paired=TRUE)
bf_n1u1_old <- ttestBF(x=precs_old$u1, y=precs_old$n1, paired=TRUE)
bf_n2c2_old <- ttestBF(x=precs_old$c2, y=precs_old$n2, paired=TRUE)
bf_n2u2_old <- ttestBF(x=precs_old$u2, y=precs_old$n2, paired=TRUE)


bf_n1c1_new <- ttestBF(x=precs_new$c1, y=precs_new$n1, paired=TRUE)
bf_n1u1_new <- ttestBF(x=precs_new$u1, y=precs_new$n1, paired=TRUE)
bf_n2c2_new <- ttestBF(x=precs_new$c2, y=precs_new$n2, paired=TRUE)
bf_n2u2_new <- ttestBF(x=precs_new$u2, y=precs_new$n2, paired=TRUE)