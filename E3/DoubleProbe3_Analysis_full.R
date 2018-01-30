# Analysis of Experiment 3 for the DoubleProbe JEP:HPP paper - SRC
#### install necessary packages -----

libs <- c('pracma', 'tidyverse', 'ggplot2', 'circular', 'circular', 'afex', 'BayesFactor', 'magrittr')
invisible(lapply(libs, require, character.only = T))
theme_set(theme_bw() +
            theme(axis.text   = element_text(size = 12),
                  axis.title  = element_text(size = 12),
                  title       = element_text(size = 12),
                  legend.text = element_text(size = 12)
            ))
#load in the data
path <- "/Users/user/Desktop/Experiments/Nick/DoubleProbe/E3/data/ExpDat/allData_E3.csv" #change your path here!!!
df   <- read.csv(path, header=TRUE, as.is=TRUE)
df <- dplyr::select(df, -X)

df %<>% dplyr::mutate(rdif1 = 0, rdif2 = 0, cond = 0) #set to zero if present, create if not present already
df %<>% dplyr::mutate(cond = ifelse(cues==1 & pord==1, 'cued1', cond),
                      cond = ifelse(cues==1 & pord==2, 'cued2', cond),
                      cond = ifelse(cues==0, 'neutral', cond))

df %<>% dplyr::mutate(rdif1 = as.circular(theta_1 - rad(angs_1), units="radians")) %>%#don't change this!!
  dplyr::mutate(rdif2 = as.circular(theta_2 - rad(angs_2), units="radians")) #don't change this!!

sublist <- c(2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 22, 24)
se <- function(x) sd(x)/sqrt(length(x))
precs        <- df %>% dplyr::group_by(subID, btype, cond) %>% summarise_at(c("rdif1","rdif2"),funs(1/sd(.)))
precs        %<>% dplyr::filter(subID %in% sublist)
precs_old    <- precs %>% dplyr::filter(btype ==0) %>% dplyr::group_by(cond) %>% dplyr::summarise_at(c('rdif1', 'rdif2'), .funs=mean) %>% melt() %>% dplyr::rename(mean = value)
precs_old_se <- precs %>% dplyr::filter(btype ==0) %>% dplyr::group_by(cond) %>% dplyr::summarise_at(c('rdif1', 'rdif2'), .funs=se)   %>% melt() %>% dplyr::rename(sem  = value)
precs_new    <- precs %>% dplyr::filter(btype ==1) %>% dplyr::group_by(cond) %>% dplyr::summarise_at(c('rdif1', 'rdif2'), .funs=mean) %>% melt() %>% dplyr::rename(mean = value)
precs_new_se <- precs %>% dplyr::filter(btype ==1) %>% dplyr::group_by(cond) %>% dplyr::summarise_at(c('rdif1', 'rdif2'), .funs=se)   %>% melt() %>% dplyr::rename(sem  = value)

plot.precs.old <- precs_old %>% dplyr::left_join(precs_old_se)
plot.precs.new <- precs_new %>% dplyr::left_join(precs_new_se)

plot.precs.old %<>%
  dplyr::mutate(variable = ifelse(variable=='rdif1', 'first', 'second')) %>% dplyr::rename(probe = variable) %>%
  dplyr::mutate(cond = ifelse(cond == 'cued1' & probe == 'first' , 'cued'  , cond),
                cond = ifelse(cond == 'cued2' & probe == 'first' , 'uncued', cond),
                cond = ifelse(cond == 'cued1' & probe == 'second', 'uncued', cond),
                cond = ifelse(cond == 'cued2' & probe == 'second', 'cued'  , cond))
plot.precs.new %<>%
  dplyr::mutate(variable = ifelse(variable=='rdif1', 'first', 'second')) %>% dplyr::rename(probe = variable) %>%
  dplyr::mutate(cond = ifelse(cond == 'cued1' & probe == 'first' , 'cued'  , cond),
                cond = ifelse(cond == 'cued2' & probe == 'first' , 'uncued', cond),
                cond = ifelse(cond == 'cued1' & probe == 'second', 'uncued', cond),
                cond = ifelse(cond == 'cued2' & probe == 'second', 'cued'  , cond))

#ggplot old blocktype
plot.precs.old %>% 
  ggplot(aes(x = probe, y = mean)) +
  geom_bar(aes(fill=cond), stat = 'identity', color = 'black', position = position_dodge()) +
  scale_fill_manual(values = c("neutral" = "#bdbdbd","cued" = "#1464f4","uncued" = "#deebf7")) +
  geom_errorbar(aes(x = probe, ymin=mean-sem, ymax=mean+sem, fill=cond), position = position_dodge(.9), width = .2) +
  scale_colour_manual(values = c("cued" = "#000000" , "neutral" = "#000000", "uncued" = "#000000")) +
  labs(x = 'Probe', y = "Accuracy (1/SD)", title = 'order not cued')
ggsave(filename='~/Desktop/Experiments/Nick/DoubleProbe/E3/figures/NoOrder_accuracy.eps', width = 12, height = 8)

plot.precs.new %>% 
  ggplot(aes(x = probe, y = mean)) +
  geom_bar(aes(fill=cond), stat = 'identity', color = 'black', position = position_dodge()) +
  scale_fill_manual(values = c("neutral" = "#bdbdbd","cued" = "#1464f4","uncued" = "#deebf7")) +
  geom_errorbar(aes(x = probe, ymin=mean-sem, ymax=mean+sem, fill=cond), position = position_dodge(.9), width = .2) +
  scale_colour_manual(values = c("cued" = "#000000" , "neutral" = "#000000", "uncued" = "#000000")) +
  labs(x = 'Probe', y = "Accuracy (1/SD)", title = 'order cued')
ggsave(filename='~/Desktop/Experiments/Nick/DoubleProbe/E3/figures/OrderCue_accuracy.eps', width = 12, height = 8)


#------------------------------------------------------------------------------------------------------------------------
# quickly check uncued costs in both block types, and compare
cohens_d <- function(x,y){
  m1 <- base::mean(x); l1 = length(x); variance1 <- stats::var(x); 
  m2 <- base::mean(y); l2 = length(y); variance2 <- stats::var(y)
  sigma <- (sqrt((l1-1)*variance1 + (l2-1)*variance2)/(l1 + l2 - 2))
  d <- (m1-m2)/sigma
  return(d)
}
t.precs.old <- precs %>% melt(., id.vars = c('subID', 'cond', 'btype')) %>%
  dplyr::filter(btype == 0)                                             %>%
  dplyr::select(subID, cond, variable,value)                            %>%
  dplyr::rename(probe = variable)                                       %>% 
  dplyr::mutate(probe = ifelse(probe == 'rdif1', 'first', 'second'))    %>%
  acast(subID~probe+cond, value.var = 'value')                          %>%
  as.data.frame()

t.precs.new <- precs %>% melt(., id.vars = c('subID', 'cond', 'btype')) %>%
  dplyr::filter(btype == 1)                                             %>%
  dplyr::select(subID, cond, variable, value)                           %>%
  dplyr::rename(probe = variable)                                       %>%
  dplyr::mutate(probe = ifelse(probe == 'rdif1', 'first', 'second'))    %>%
  acast(subID~probe+cond, value.var = 'value')                          %>%
  as.data.frame()

t.precs.old %<>% dplyr::rename(c1 = first_cued1, n1 = first_neutral, u1 = first_cued2, c2 = second_cued2, n2 = second_neutral, u2 = second_cued1)
t.precs.new %<>% dplyr::rename(c1 = first_cued1, n1 = first_neutral, u1 = first_cued2, c2 = second_cued2, n2 = second_neutral, u2 = second_cued1)

  
uncuedcost_old <- t.precs.old %>% dplyr::select(-c1, -c2) %>% dplyr::mutate(cost1 = n1-u1, cost2 = n2-u2) %>% dplyr::select(cost1,cost2)
uncuedcost_new <- t.precs.new %>% dplyr::select(-c1, -c2) %>% dplyr::mutate(cost1 = n1-u1, cost2 = n2-u2) %>% dplyr::select(cost1,cost2) 

t.test(uncuedcost_old$cost1, uncuedcost_old$cost2, paired = TRUE); ttestBF(x = uncuedcost_old$cost1, y = uncuedcost_old$cost2, paired = TRUE); cohens_d(uncuedcost_old$cost1, uncuedcost_old$cost2)
t.test(uncuedcost_new$cost1, uncuedcost_new$cost2, paired = TRUE); ttestBF(x = uncuedcost_new$cost1, y = uncuedcost_new$cost2, paired = TRUE); cohens_d(uncuedcost_new$cost1, uncuedcost_new$cost2)

t.test(t.precs.old$n1, t.precs.old$n2, paired = TRUE)
t.test(t.precs.new$n1, t.precs.new$n2, paired = TRUE)
#

# create long format dframes for ANOVAs on the experiment
aov.precs <- melt(precs, id.vars = 1:3) %>%
  dplyr::rename(probe = variable, accuracy = value, cue = cond) %>%
  dplyr::mutate(probe = ifelse(probe == 'rdif1', 'first', 'second'),
                cue  = ifelse(cue == 'cued1' & probe == 'first', 'cued', cue),
                cue  = ifelse(cue == 'cued2' & probe == 'second', 'cued', cue),
                cue  = ifelse(cue == 'cued1' & probe == 'second', 'uncued', cue),
                cue  = ifelse(cue == 'cued2' & probe == 'first', 'uncued', cue)) %>%
  dplyr::mutate(cue = as.factor(cue), probe = as.factor(probe), btype = as.factor(btype), subID = as.factor(subID))

aov.precs.old <- aov.precs %>% dplyr::filter(btype == 0) %>% dplyr::select(-btype)
aov.precs.new <- aov.precs %>% dplyr::filter(btype == 1) %>% dplyr::select(-btype)

aov_old <- aov_ez('subID', 'accuracy', aov.precs.old, within = c("cue", "probe"))              ; nice(aov_old, es='pes')
aov_new <- aov_ez('subID', 'accuracy', aov.precs.new, within = c("cue", "probe"))              ; nice(aov_new, es='pes')
aov_all <- aov_ez('subID', 'accuracy', aov.precs    , within = c("cue", "probe", "btype")); nice(aov_all, es='pes')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ANOVA of cueing benefits and costs (i.e. removing one level of cue factor)

#create dataframes for the analysis
aov.precs.neutralcued   <- aov.precs %>% dplyr::filter(cue != 'uncued') %>% droplevels()
aov.precs.neutraluncued <- aov.precs %>% dplyr::filter(cue != 'cued')   %>% droplevels()

#analysis of neutral vs cued across both block types
aov_neutralcued <- aov.precs.neutralcued %>% aov_ez('subID', 'accuracy', ., within = c('cue', 'probe', 'btype')); nice(aov_neutralcued, es = 'pes')
bf_neutralcued  <- anovaBF(accuracy ~ subID + cue*probe*btype, data = aov.precs.neutralcued, whichRandom = 'subID', iterations = 200000, progress = T, whichModels = 'withmain')

#analysis of neutral vs uncued across both block types
aov_neutraluncued <- aov.precs.neutraluncued %>% aov_ez('subID','accuracy', ., within = c('cue','probe', 'btype')); nice(aov_neutraluncued, es = 'pes')
bf_neutraluncued <- anovaBF(accuracy ~ subID + cue*probe*btype, data = aov.precs.neutraluncued, whichRandom = 'subID', iterations = 200000, progress=T, whichModels = 'withmain')

#analysis of just the old block type
aov_neutralcued_old   <- aov.precs.neutralcued %>% dplyr::filter(btype == 0) %>% droplevels() %>% aov_ez('subID', 'accuracy', ., within = c('cue', 'probe')); nice(aov_neutralcued_old  , es = 'pes')
bf_neutralcued_old    <- aov.precs.neutralcued %>%  

aov_neutraluncued_old <- aov.precs.neutraluncued %>% dplyr::filter(btype == 0) %>% droplevels() %>% aov_ez('subID', 'accuracy', ., within = c('cue', 'probe')); nice(aoc_neutraluncued_old, es = 'pes')
bf_neutraluncued_old  <-

#ttests of the old block type, and effect sizes
n1c1_old <- t.test(t.precs.old$n1, t.precs.old$c1, paired=TRUE); n1c1old_effsize <- cohens_d(t.precs.old$n1, t.precs.old$c1)
n2c2_old <- t.test(t.precs.old$n2, t.precs.old$c2, paired=TRUE); n2c2old_effsize <- cohens_d(t.precs.old$n2, t.precs.old$c2)
n1u1_old <- t.test(t.precs.old$n1, t.precs.old$u1, paired=TRUE); n1u1old_effsize <- cohens_d(t.precs.old$n1, t.precs.old$u1)
n2u2_old <- t.test(t.precs.old$n2, t.precs.old$u2, paired=TRUE); n2u2old_effsize <- cohens_d(t.precs.old$n2, t.precs.old$u2)

#bayesian t-tests for old block type
bf.n1c1.old <- ttestBF(t.precs.old$n1, t.precs.old$c1, paired = T)
bf.n2c2.old <- ttestBF(t.precs.old$n2, t.precs.old$c2, paired = T)
bf.n1u1.old <- ttestBF(t.precs.old$n1, t.precs.old$u1, paired = T)
bf.n2u2.old <- ttestBF(t.precs.old$n2, t.precs.old$u2, paired = T)


# # same analysis pipelines for the new block type
aov_neutralcued_new   <- aov.precs.neutralcued %>% dplyr::filter(btype == 1) %>% droplevels() %>% aov_ez('subID', 'accuracy', ., within = c('cue', 'probe')); nice(aov_neutralcued_new, es='pes')
bf_neutralcued_new    <- aov.precs.neutralcued %>% dplyr::filter(btype == 1) %>% droplevels() %>% # remove old block data, and drop levels of the factor
  anovaBF(accuracy ~ subID + cue*probe, data = ., whichRandom = 'subID', iterations = 200000, progress = T, whichModels = 'withmain') # run bayesian anova

aov_neutraluncued_new <- aov.precs.neutraluncued %>% dplyr::filter(btype == 1) %>% droplevels() %>% aov_ez('subID', 'accuracy', ., within = c('cue', 'probe')); nice(aov_neutraluncued_new, es = 'pes')
bf_neutraluncued_new  <- aov.precs.neutraluncued %>% dplyr::filter(btype == 1) %>% droplevels() %>% # remove old block data, and drop levels of the factor
  anovaBF(accuracy ~ subID + cue*probe, data = ., whichRandom = 'subID', iterations = 200000, progress = T, whichModels = 'withmain') # run bayesian anova

#ttests of the new block type, and effect sizes
n1c1_new <- t.test(t.precs.new$n1, t.precs.new$c1, paired = T); n1c1_new_effsize <- cohens_d(t.precs.new$n1, t.precs.new$c1)
n2c2_new <- t.test(t.precs.new$n2, t.precs.new$c2, paired = T); n2c2_new_effsize <- cohens_d(t.precs.new$n2, t.precs.new$c2)
n1u1_new <- t.test(t.precs.new$n1, t.precs.new$u1, paired = T); n1u1_new_effsize <- cohens_d(t.precs.new$n1, t.precs.new$u1)
n2u2_new <- t.test(t.precs.new$n2, t.precs.new$u2, paired = T); n2u2_new_effsize <- cohens_d(t.precs.new$n2, t.precs.new$u2)

#bayesian t-tests for new block type
bf.n1c1.new <- ttestBF(t.precs.new$n1, t.precs.new$c1, paired = T)
bf.n2c2.new <- ttestBF(t.precs.new$n2, t.precs.new$c2, paired = T)
bf.n1u1.new <- ttestBF(t.precs.new$n1, t.precs.new$u1, paired = T)
bf.n2u2.new <- ttestBF(t.precs.new$n2, t.precs.new$u2, paired = T)

#------------------------------------------------------------------------------------------------------------------------
# Bayesian statistics (for getting bayes factors for interactions/main effects)

bf_all  <- anovaBF(accuracy ~ subID + cue*probe*btype, data = aov.precs    , whichRandom = "subID", iterations = 200000, progress = T, whichModels = "withmain")
bf_old  <- anovaBF(accuracy ~ subID + cue*probe      , data = aov.precs.old, whichRandom = "subID", iterations = 200000, progress = T, whichModels = "withmain")
bf_new  <- anovaBF(accuracy ~ subID + cue*probe      , data = aov.precs.new, whichRandom = "subID", iterations = 200000, progress = T, whichModels = "withmain")

bf_all[3]/bf_all[8] #bayes factor for argument against an Experiment effect
#i.e. yields BF of ~5.7, strong evidence that there is no effect of Experiment on Precision (performance no significantly better or worse between block types)
# a model without an effect of experiment fares 6.58 times better

lmbf_old_full <- lmBF(accuracy ~ subID + cue + probe + cue:probe, data = aov.precs.old, whichRandom= 'subID', iterations = 200000)
lmbf_old_main <- lmBF(accuracy ~ subID + cue + probe,             data = aov.precs.old, whichRandom= 'subID', iterations = 200000)
lmbf_old_full/lmbf_old_main #evidence for interaction of cue and probe factors in Experiment 3 blocktype = order not cued

# calculate bayes factor for cue*probe interaction in cued and neutral trials when order was not cued (blocktype == 0)
lmbf_old_NC_full <- aov.precs.neutralcued %>% dplyr::filter(btype == 0) %>% droplevels() %>% #make the dataframe for the analysis (take only old block type)
  lmBF(accuracy ~ subID + cue + probe + cue:probe, data = ., whichRandom = 'subID', iterations  = 200000) #run the lmBF
lmbf_old_NC_main <- aov.precs.neutralcued %>% dplyr::filter(btype == 0) %>% droplevels() %>% #make the dataframe for the analysis (take only old block type)
  lmBF(accuracy ~ subID + cue + probe, data = ., whichRandom = 'subID', iterations  = 200000) #run the lmBF
lmbf_old_NC_full/lmbf_old_NC_main # evidence for cue:probe interaction when comparing cue and probe in only neutral and cued trials, blocktype = order not cued


#calculate bayes factor for cue*probe interaction in neutral and uncued items only, when order was not cued (blocktype == 0)
lmbf_old_NU_full <- aov.precs.neutraluncued %>% dplyr::filter(btype == 0) %>% droplevels() %>% #make the dataframe for the analysis (take only old block type)
  lmBF(accuracy ~ subID + cue + probe + cue:probe,  data = ., whichRandom = 'subID', iterations  = 200000)
lmbf_old_NU_main <- aov.precs.neutraluncued %>% dplyr::filter(btype == 0) %>% droplevels() %>% #create datarame for the analysis and pipe it into the lmBF
  lmBF(accuracy ~ subID + cue + probe, data = ., whichRandom = 'subID', iterations  = 200000)
lmbf_old_NU_full/lmbf_old_NU_main # evidence for cue:probe interaction when comparing cue and probe in only neutral and uncued items, blocktype = order not cued


#calculate bayes factor for cue*probe interaction where order was cued (blocktype == 1)
lmbf_new_full <- lmBF(accuracy ~ subID + cue + probe + cue:probe, data = aov.precs.new, whichRandom= 'subID', iterations = 200000)
lmbf_new_main <- lmBF(accuracy ~ subID + cue + probe,             data = aov.precs.new, whichRandom= 'subID', iterations = 200000)
lmbf_new_full/lmbf_new_main #evidence for interaction of cue and probe factors in Experiment 3 blocktype = order not cued

lmbf_new_NC_full <- aov.precs.neutralcued %>% dplyr::filter(btype == 1) %>% droplevels() %>% #create dataframe to pipe into lmBF (drop old blocktype)
  lmBF(accuracy ~ subID + cue + probe + cue:probe, data = ., whichRandom = 'subID', iterations  = 200000)
lmbf_new_NC_main <- aov.prevs.neutralcued %>% dplyr::filter(btype == 1) %>% droplebels() %>% #create dataframe to pipe into lmBF (drop old blocktype)
  lmBF(accuracy ~ subID + cue + probe, data = ., whichRandom = 'subID', iterations  = 200000)
lmbf_new_NC_full/lmbf_new_NC_main # evidence for cue:probe interaction when comparing cue and probe in only neutral and cued trials, blocktype = order  cued

lmbf_new_NU_full <- aov.precs.neutraluncued %>% dplyr::filter(btype == 1) %>% droplevels() %>% #create dataframe to pipe into lmBF (drop old blocktype)
  lmBF(accuracy ~ subID + cue + probe + cue:probe,  data = ., whichRandom = 'subID', iterations  = 200000)
lmbf_new_NU_main <- aov.precs.neutraluncued %>% dplyr::filter(btype == 1) %>% droplevels() %>% #create dataframe to pipe into lmBF (drop old blocktype)
  lmBF(accuracy ~ subID + cue + probe, data = ., whichRandom = 'subID', iterations  = 200000)
lmbf_new_NU_full/lmbf_new_NU_main # evidence for cue:probe interaction when comparing cue and probe in only neutral and uncued trials, blocktype = order  cued