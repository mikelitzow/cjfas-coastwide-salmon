# This script uses mixed-effects models to test for changing SST effects 
# on salmon by era, fitting separate models to each species, combining data over all three regions.
# Uses the approach outlined by Zuur et al. 2009.

library(MuMIn)
library(dplyr)
library(ggplot2)
library(nlme)
library(lmtest)
library(tidyverse)
options(max.print = 99999)


run.dat <- read.csv("coastwide salmon data.csv", row.names = 1)

# make sure that era is a factor
run.dat$era <- as.factor(run.dat$era)

# Based on the results in Mueter et al. 2002 CFAS, Mueter et al. 2005 TAFS and the relationships in this correlation plot,
# we'll define the months to use *three different ways*. 
# Support for the non-stationary hypothesis will be supported for each set of months,
# and the best set of months (SST time window) will be selected with AICc model selection. 
# (Note that there are three SST groupings for sockeye, reflecting their more complex life history, 
# but only 2 for pink & chum).
# 
# SST1: 
# Jan-Oct of BY + 1 for chum and pink,
# Jan BY + 1 through May BY + 3 for AK sockeye.
# Jan BY + 1 through May BY + 2 for southern sockeye.
# 
# SST2:
# Jan BY + 1 through Oct BY + 2/3 for sockeye, 
# with Alaskan sockeye values determined by the average age at ocean entry for each run.
# 
# SST3 (summer only, corresponding to values in Mueter et al. 2002):
# April-July for Washington / BC / SEAK, May-August for remainder of GOA, June-September for EBS.
# Pink and Chum BY+1, South sockeye BY+2, AK sockeye BY+2/3 based on average age of ocean entry.

# SST1, SST2, and SST3 correspond with k1, k2 and k3 in Table 1 of the ms.

# first, define the full model with all fixed effects
formula.full <- 
as.formula(log(recruits/spawners) ~ 1 + stock:spawners +  sst1:region + sst1:region:era)

# now select best random structure - all models throughout include autocorrelated residuals
# fixed effects model
p1 <- gls(formula.full, method = "REML", 
correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Pink"))

# random intercept for stocks
p2 <- lme(formula.full, random = ~1 | stock,
method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Pink"))

# full interaction with intercept
p3 <- lme(formula.full, random = ~1 + sst1:region:era | stock,
method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Pink"),
control=lmeControl(opt="optim", maxIter = 200, msMaxIter = 200))

# reduced interaction with intercept
p4 <- lme(formula.full, random = ~1 + sst1:region | stock,
method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Pink"),
control=lmeControl(opt="optim", maxIter = 200, msMaxIter = 200))

# full interaction, no intercept
p5 <- lme(formula.full, random = ~ -1 + sst1:region:era | stock,
method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Pink"),
control=lmeControl(opt="optim", maxIter = 200, msMaxIter = 200))

# reduced interaction, no intercept
p6 <- lme(formula.full, random = ~-1 + sst1:region | stock,
method = "REML", correlation=corAR1(form= ~1 | stock), data = filter(run.dat, species=="Pink"))

# compare with AICc
AIC <- AICc(p1, p2, p3, p4, p5, p6) 
AIC <- AIC %>%
mutate(model=rownames(AIC)) %>%
mutate(dAICc=AICc-min(AICc)) %>%
arrange(AICc)
AIC

# Model p2 (random intercepts random structure) is the best model by far. 
# Look at the residuals:

qqnorm(p2, ~ranef(., level=1))

# pretty good 

# Now the fixed effects structure.

# full model is p2 from above - but changing to ML estimation
pfull <- lme(formula.full, random = ~1 | stock,
method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Pink"))

formula.reduced <- 
as.formula(log(recruits/spawners) ~ 1 + stock:spawners + sst1:region)

p1 <- lme(formula.reduced, random = ~1 | stock,
method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Pink"))

anova(pfull, p1)

# The sst1:region:era term is strongly supported. We'll compare with the model without autocorrelated residuals.

# alternate model without AR(1) residuals
p2 <- lme(formula.full, random = ~1 | stock,
          method = "ML", data = filter(run.dat, species=="Pink"))

anova(pfull, p2)

# So there isn't any big difference in the AIC scores, but p2 is better, 
# though there is no support in the likellihood ration test
# for retaining the autocorrelated structure. 
# For uniformity with the other spp., where AR(1) residuals are well supported, 
# we'll retain that structure here. Also, will save as p.sst1 for model comparison.

p.sst1 <- pfull

# Significant sst1:era effects for GOA and South.

# And now the plot.
pplot <- as.data.frame(intervals(pfull, which="fixed")$fixed[c(37:42),])

pplot$term <- rownames(pplot)

pplot$term <- reorder(pplot$term, pplot$est.)

ggplot(pplot, aes(term, est.)) + geom_errorbar(aes(ymax=upper, ymin=lower), width=0.2) + geom_point(size=3) + coord_flip() +
  xlab("") + ylab("Coefficient") +  ggtitle("Pink salmon") + geom_hline(yintercept = 0) + 
  theme(axis.text.y  = element_text(size=10))


# Continue looking at SST1 - now with Sockeye salmon

# Same routine - begin with random effects. 

formula.full <- 
  as.formula(log(recruits/spawners) ~ 1 + stock:spawners + sst1:region + sst1:region:era)

# fixed effects model
s1 <- gls(formula.full, method = "REML", 
          correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"))

# random intercept for stocks
s2 <- lme(formula.full, random = ~1 | stock,
          method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"))

# full interaction with intercept
s3 <- lme(formula.full, random = ~1 + sst1:region:era | stock,
          method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"),
          control=lmeControl(opt="optim", maxIter = 200, msMaxIter = 200))

# reduced interaction with intercept
s4 <- lme(formula.full, random = ~1 + sst1:region | stock,
          method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"),
          control=lmeControl(opt="optim", maxIter = 200, msMaxIter = 200))

# full interaction, no intercept
s5 <- lme(formula.full, random = ~ -1 + sst1:region:era | stock,
          method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"),
          control=lmeControl(opt="optim", maxIter = 200, msMaxIter = 200))

# reduced interaction, no intercept
s6 <- lme(formula.full, random = ~-1 + sst1:region | stock,
          method = "REML", correlation=corAR1(form= ~1 | stock), data = filter(run.dat, species=="Sockeye"))

AIC <- AICc(s1, s2, s3, s4, s5, s6) 
AIC <- AIC %>%
  mutate(model=rownames(AIC)) %>%
  mutate(dAICc=AICc-min(AICc)) %>%
  arrange(AICc)
AIC

# Again, model 2 massively better. Here are the random effect residuals.
qqnorm(s2, ~ranef(., level=1))

# Now the hypothesis test via comparing fixed effects.

# full model from above - but changing to ML estimation
sfull <- lme(formula.full, random = ~1 | stock,
             method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"))

formula.reduced <- 
  as.formula(log(recruits/spawners) ~ 1 + stock:spawners + sst1:region)

s1 <- lme(formula.reduced, random = ~1 | stock,
          method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"))

anova(sfull, s1)

# The sst1:region:era term is supported. Compare with the model without autocorrelated residuals.

s2 <- lme(formula.full, random = ~1 | stock,
method = "ML", data = filter(run.dat, species=="Sockeye"))

anova(sfull, s2)

# In this case the AR(1) errors are far better. Save ML version for comparison below.
s.sst1 <- sfull

# And the plot.

splot <- as.data.frame(intervals(sfull, which="fixed")$fixed[c(32:37),])

splot$term <- rownames(splot)

splot$term <- reorder(splot$term, splot$est.)

ggplot(splot, aes(term, est.)) + geom_errorbar(aes(ymax=upper, ymin=lower), width=0.2) + geom_point(size=3) + coord_flip() +
xlab("") + ylab("Coefficient") +  ggtitle("Sockeye salmon") + geom_hline(yintercept = 0) + 
theme(axis.text.y  = element_text(size=10))

# Now SST1 and chum salmon

# Random effects selection.

formula.full <- 
as.formula(log(recruits/spawners) ~ 1 + stock:spawners + sst1:region + sst1:region:era)

# fixed effects model
c1 <- gls(formula.full, method = "REML", 
correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Chum"))

# random intercept for stocks
c2 <- lme(formula.full, random = ~1 | stock,
method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Chum"))

# full interaction with intercept - DOESN'T FIT!
  # c3 <- lme(formula.full, random = ~1 + sst1:region:era | stock,
  #           method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Chum"),
  #           control=lmeControl(opt="optim", maxIter = 200, msMaxIter = 200))
  
  # reduced interaction with intercept
  c4 <- lme(formula.full, random = ~1 + sst1:region | stock,
            method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Chum"),
            control=lmeControl(opt="optim", maxIter = 200, msMaxIter = 200))

# full interaction, no intercept
c5 <- lme(formula.full, random = ~ -1 + sst1:region:era | stock,
          method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Chum"),
          control=lmeControl(opt="optim", maxIter = 200, msMaxIter = 200))

# reduced interaction, no intercept
c6 <- lme(formula.full, random = ~-1 + sst1:region | stock,
          method = "REML", correlation=corAR1(form= ~1 | stock), data = filter(run.dat, species=="Chum"))

AIC <- AICc(c1, c2, c4, c5, c6) 
AIC <- AIC %>%
  mutate(model=rownames(AIC)) %>%
  mutate(dAICc=AICc-min(AICc)) %>%
  arrange(AICc)
AIC

# Model 2 again the best.

qqnorm(c2, ~ranef(., level=1))

# Now select fixed effects.

# full model - changing to ML estimation
cfull <- lme(formula.full, random = ~1 | stock,
             method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Chum"))

formula.reduced <- 
  as.formula(log(recruits/spawners) ~ 1 + stock:spawners + sst1:region)

c1 <- lme(formula.reduced, random = ~1 | stock,
          method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Chum"))

anova(cfull, c1)

# Nonstationary model supported. Compare different error structures.

cfull <- lme(formula.full, random = ~1 | stock,
             method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Chum"))

c2 <- lme(formula.full, random = ~1 | stock,
          method = "ML", data = filter(run.dat, species=="Chum"))

anova(cfull, c2)

# Again, the AR(1) errors are far better!
  
# Save the ML version.
c.sst1 <- cfull

# And plot.

cplot <- as.data.frame(intervals(cfull, which="fixed")$fixed[c(23:28),])

cplot$term <- rownames(cplot)

cplot$term <- reorder(cplot$term, cplot$est.)

ggplot(cplot, aes(term, est.)) + geom_errorbar(aes(ymax=upper, ymin=lower), width=0.2) + geom_point(size=3) + coord_flip() +
  xlab("") + ylab("Coefficient") +  ggtitle("Chum salmon") + geom_hline(yintercept = 0) + 
  theme(axis.text.y  = element_text(size=10))

# Now repeat with SST2. This SST grouping is only tested for sockeye.

# Same routine - begin with random effects. 

formula.full <- 
as.formula(log(recruits/spawners) ~ 1 + stock:spawners + sst2:region + sst2:region:era)

# fixed effects model
s1 <- gls(formula.full, method = "REML", 
correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"))

# random intercept for stocks
s2 <- lme(formula.full, random = ~1 | stock,
method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"))

# full interaction with intercept
s3 <- lme(formula.full, random = ~1 + sst2:region:era | stock,
method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"),
control=lmeControl(opt="optim", maxIter = 200, msMaxIter = 200))

# reduced interaction with intercept
s4 <- lme(formula.full, random = ~1 + sst2:region | stock,
method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"),
control=lmeControl(opt="optim", maxIter = 200, msMaxIter = 200))

# full interaction, no intercept
s5 <- lme(formula.full, random = ~ -1 + sst2:region:era | stock,
method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"),
control=lmeControl(opt="optim", maxIter = 200, msMaxIter = 200))

# reduced interaction, no intercept
s6 <- lme(formula.full, random = ~-1 + sst2:region | stock,
method = "REML", correlation=corAR1(form= ~1 | stock), data = filter(run.dat, species=="Sockeye"))

AIC <- AICc(s1, s2, s3, s4, s5, s6) 
AIC <- AIC %>%
mutate(model=rownames(AIC)) %>%
mutate(dAICc=AICc-min(AICc)) %>%
arrange(AICc)
AIC

# Same random structure is best.

qqnorm(s2, ~ranef(., level=1))

# Now the hypothesis test via comparing fixed effects.

# full model from above - but changing to ML estimation
sfull <- lme(formula.full, random = ~1 | stock,
method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"))

formula.reduced <- 
as.formula(log(recruits/spawners) ~ 1 + stock:spawners + sst2:region)
s1 <- lme(formula.reduced, random = ~1 | stock,
method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"))

anova(sfull, s1)

# Nonstationary model supported, but not as strongly. Save for model comparison.

s.sst2 <- sfull

# Finally, SST3.

# Pink salmon
# Begin by fitting the random effects. 

formula.full <- 
as.formula(log(recruits/spawners) ~ 1 + stock:spawners + sst3:region + sst3:region:era)

# fixed effects model
p1 <- gls(formula.full, method = "REML", 
correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Pink"))

# random intercept for stocks
p2 <- lme(formula.full, random = ~1 | stock,
method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Pink"))

# full interaction with intercept
p3 <- lme(formula.full, random = ~1 + sst3:region:era | stock,
method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Pink"),
control=lmeControl(opt="optim", maxIter = 200, msMaxIter = 200))

# reduced interaction with intercept
p4 <- lme(formula.full, random = ~1 + sst3:region | stock,
method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Pink"),
control=lmeControl(opt="optim", maxIter = 200, msMaxIter = 200))

# full interaction, no intercept
p5 <- lme(formula.full, random = ~ -1 + sst3:region:era | stock,
method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Pink"),
control=lmeControl(opt="optim", maxIter = 200, msMaxIter = 200))

# reduced interaction, no intercept
p6 <- lme(formula.full, random = ~-1 + sst3:region | stock,
method = "REML", correlation=corAR1(form= ~1 | stock), data = filter(run.dat, species=="Pink"))

AIC <- AICc(p1, p2, p3, p4, p5, p6) 
AIC <- AIC %>%
mutate(model=rownames(AIC)) %>%
mutate(dAICc=AICc-min(AICc)) %>%
arrange(AICc)
AIC

# Random intercepts random structure is the best model by far. Look at the residuals:

qqnorm(p2, ~ranef(., level=1))

# Now the fixed effects structure.

# full model is p2 from above - but changing to ML estimation
pfull <- lme(formula.full, random = ~1 | stock,
method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Pink"))

formula.reduced <- 
as.formula(log(recruits/spawners) ~ 1 + stock:spawners + sst3:region)

p1 <- lme(formula.reduced, random = ~1 | stock,
method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Pink"))

anova(pfull, p1)

# The sst3:region:era term is also strongly supported. Save best model for comparison.

p.sst3 <- pfull

# Sockeye salmon

# Same routine - begin with random effects. 

formula.full <- 
as.formula(log(recruits/spawners) ~ 1 + stock:spawners + sst3:region + sst3:region:era)

# fixed effects model
s1 <- gls(formula.full, method = "REML", 
correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"))

# random intercept for stocks
s2 <- lme(formula.full, random = ~1 | stock,
method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"))

# full interaction with intercept
s3 <- lme(formula.full, random = ~1 + sst3:region:era | stock,
method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"),
control=lmeControl(opt="optim", maxIter = 200, msMaxIter = 200))

# reduced interaction with intercept
s4 <- lme(formula.full, random = ~1 + sst3:region | stock,
method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"),
control=lmeControl(opt="optim", maxIter = 200, msMaxIter = 200))

# full interaction, no intercept
s5 <- lme(formula.full, random = ~ -1 + sst3:region:era | stock,
method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"),
control=lmeControl(opt="optim", maxIter = 200, msMaxIter = 200))

# reduced interaction, no intercept
s6 <- lme(formula.full, random = ~-1 + sst3:region | stock,
method = "REML", correlation=corAR1(form= ~1 | stock), data = filter(run.dat, species=="Sockeye"))

AIC <- AICc(s1, s2, s3, s4, s5, s6) 
AIC <- AIC %>%
mutate(model=rownames(AIC)) %>%
mutate(dAICc=AICc-min(AICc)) %>%
arrange(AICc)
AIC

# Same random structure is best.

qqnorm(s2, ~ranef(., level=1))

# Now the hypothesis test via comparing fixed effects.

# full model from above - but changing to ML estimation
sfull <- lme(formula.full, random = ~1 | stock,
method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"))

formula.reduced <- 
as.formula(log(recruits/spawners) ~ 1 + stock:spawners + sst3:region)

s1 <- lme(formula.reduced, random = ~1 | stock,
method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"))

anova(sfull, s1)

# Nonstationary model supported, but more weakly. Save for model comparison.

s.sst3 <- sfull

# Chum salmon

# Random effects.

formula.full <- 
as.formula(log(recruits/spawners) ~ 1 + stock:spawners + sst3:region + sst3:region:era)

# fixed effects model
c1 <- gls(formula.full, method = "REML", 
correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Chum"))

# random intercept for stocks
c2 <- lme(formula.full, random = ~1 | stock,
method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Chum"))

# full interaction with intercept - DOESN'T FIT!
  # c3 <- lme(formula.full, random = ~1 + sst3:region:era | stock,
  #           method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Chum"),
  #           control=lmeControl(opt="optim", maxIter = 200, msMaxIter = 200))
  
  # reduced interaction with intercept - DOESN'T FIT!
  # c4 <- lme(formula.full, random = ~1 + sst3:region | stock,
  #          method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Chum"),
  #          control=lmeControl(opt="optim", maxIter = 200, msMaxIter = 200))
  
  # full interaction, no intercept
  c5 <- lme(formula.full, random = ~ -1 + sst3:region:era | stock,
            method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Chum"),
            control=lmeControl(opt="optim", maxIter = 200, msMaxIter = 200))

# reduced interaction, no intercept
c6 <- lme(formula.full, random = ~-1 + sst3:region | stock,
          method = "REML", correlation=corAR1(form= ~1 | stock), data = filter(run.dat, species=="Chum"))

AIC <- AICc(c1, c2, c4, c5, c6) 
AIC <- AIC %>%
  mutate(model=rownames(AIC)) %>%
  mutate(dAICc=AICc-min(AICc)) %>%
  arrange(AICc)
AIC

# Model 2 again.

qqnorm(c2, ~ranef(., level=1))

# And fixed effects.

# full model - changing to ML estimation
cfull <- lme(formula.full, random = ~1 | stock,
             method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Chum"))

formula.reduced <- 
  as.formula(log(recruits/spawners) ~ 1 + stock:spawners + sst3:region)
c1 <- lme(formula.reduced, random = ~1 | stock,
          method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Chum"))

anova(cfull, c1)

# Nonstationary model supported. Save.

c.sst3 <- cfull

# Model comparison for support of different possible SST groupings.

# Use AICc selection to compare SST1 and SST3 for each species (and SST2 for sockeye).

# First, pink.

pink.aic <- AICc(p.sst1, p.sst3) %>%
  mutate(dAICc=AICc-min(AICc), model=c("p.sst1",  "p.sst3")) %>%
  arrange(dAICc) %>%
  print()

sock.aic <- AICc(s.sst1, s.sst2, s.sst3) %>%
  mutate(dAICc=AICc-min(AICc), model=c("s.sst1", "s.sst2", "s.sst3")) %>%
  arrange(dAICc) %>%
  print()

chum.aic <- AICc(c.sst1, c.sst3) %>%
  mutate(dAICc=AICc-min(AICc), model=c("c.sst1", "c.sst3")) %>%
  arrange(dAICc) %>%
  print()

# Finally, combine into a plot for the paper, using the best-supported sst range for each. 
# Also, get the t tables for each best-supported model.

# fit the best model for each species using REML
# sst3 is the best for pinks
formula.full <- 
  as.formula(log(recruits/spawners) ~ 1 + stock:spawners + sst3:region + sst3:region:era)
pbest <- lme(formula.full, random = ~1 | stock,
             method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Pink"))

# sst1 is the best for chums and sockeye
formula.full <- 
  as.formula(log(recruits/spawners) ~ 1 + stock:spawners + sst1:region + sst1:region:era)
sbest <- lme(formula.full, random = ~1 | stock,
             method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"))
cbest <- lme(formula.full, random = ~1 | stock,
             method = "REML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Chum"))

# the t-tables
summary(pbest)$tTable[37:42,]

summary(sbest)$tTable[32:37,]

summary(cbest)$tTable[23:28,]

# And plot.

plot <- data.frame(species=rep(c("Pink", "Sockeye", "Chum"), each=6), region=rep(c("EBS", "GOA", "South"), each=2),
                   era=c("Before 1988/89", "After 1988/89"), LCI=NA, effect=NA, UCI=NA)

plot[1,4:6] <- intervals(pbest, which="fixed")$fixed[37,] # pink EBS early 
plot[2,4:6] <- intervals(pbest, which="fixed")$fixed[37,2]+ 
  intervals(pbest, which="fixed")$fixed[40,] # pink EBS late

plot[3,4:6] <- intervals(pbest, which="fixed")$fixed[38,] # pink GOA early
plot[4,4:6] <- intervals(pbest, which="fixed")$fixed[38,2]+ 
  intervals(pbest, which="fixed")$fixed[41,] # pink GOA late

plot[5,4:6] <- intervals(pbest, which="fixed")$fixed[39,] # pink South early
plot[6,4:6] <- intervals(pbest, which="fixed")$fixed[39,2]+ 
  intervals(pbest, which="fixed")$fixed[42,] # pink South late

plot[7,4:6] <- intervals(sbest, which="fixed")$fixed[32,] # sockeye EBS early
plot[8,4:6] <- intervals(sbest, which="fixed")$fixed[32,2] + 
  intervals(sbest, which="fixed")$fixed[35,] # sockeye EBS late

plot[9,4:6] <- intervals(sbest, which="fixed")$fixed[33,] # sockeye GOA early
plot[10,4:6] <- intervals(sbest, which="fixed")$fixed[33,2] + 
  intervals(sbest, which="fixed")$fixed[36,] # sockeye GOA late

plot[11,4:6] <- intervals(sbest, which="fixed")$fixed[34,] # sockeye South early
plot[12,4:6] <- intervals(sbest, which="fixed")$fixed[34,2] + 
  intervals(sbest, which="fixed")$fixed[37,] # sockeye South late

plot[13,4:6] <- intervals(cbest, which="fixed")$fixed[23,] # chum EBS early
plot[14,4:6] <- intervals(cbest, which="fixed")$fixed[23,2] + 
  intervals(cbest, which="fixed")$fixed[26,] # chum EBS late 

plot[15,4:6] <- intervals(cbest, which="fixed")$fixed[24,] # chum GOA early
plot[16,4:6] <- intervals(cbest, which="fixed")$fixed[24,2] + 
  intervals(cbest, which="fixed")$fixed[27,] # chum GOA late

plot[17,4:6] <- intervals(cbest, which="fixed")$fixed[25,] # chum South early
plot[18,4:6] <- intervals(cbest, which="fixed")$fixed[25,2] + 
  intervals(cbest, which="fixed")$fixed[28,] # chum South late

plot$era <- reorder(plot$era, -order(plot$era))

plot$species <- reorder(plot$species, rep(c(1,2,3), each=6))

dodge <- position_dodge(width=0.9)

# load the color-blind palette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# make a df for asterisks
label.df <- data.frame(species=rep(c("Pink", "Sockeye", "Chum"), each=3), 
                       region=c("EBS", "GOA", "South"), effect=c(0,0.68,0,0.5,0.5,-0.9,0,0,0.4), 
                       era=as.factor(rep(1,9)), 
                       label=c("", "***", "", "", "", "**", "", "", "**"))

mixed.plot <- ggplot(plot, aes(region, effect, fill=era)) + geom_hline(yintercept = 0, color="dark grey") + 
  geom_bar(position=dodge, stat="identity") + 
  geom_errorbar(aes(ymax=UCI, ymin=LCI), position=dodge, width=0.5) + xlab("") + 
  ylab("SST coefficient") + 
  facet_wrap(~species, scales="free_y") + 
  theme(legend.position = c(0.18, 0.85), legend.title = element_blank()) +
  scale_fill_manual(values=cb[c(2,6,2)], breaks=c("Before 1988/89", "After 1988/89", "")) + 
  geom_text(data = label.df, label = label.df$label, size=6)

mixed.plot <-ggplot(plot, aes(region, effect, fill=era)) + 
  theme_linedraw() +
  geom_bar(position=dodge, stat="identity") + geom_hline(yintercept = 0, color="black", lwd=0.3) + 
  geom_errorbar(aes(ymax=UCI, ymin=LCI), position=dodge, width=0.3, size=0.4) +
  xlab("") + ylab("SST coefficient") + 
  facet_wrap(~species, scales="free_y") + 
  theme(legend.position = c(0.11, 0.1), legend.title = element_blank(), legend.key.size = unit(0.14, "in")) +
  theme(legend.margin = margin(lm,lm,lm,lm,"cm")) +
  scale_fill_manual(values=cb[c(2,6,2)], breaks=c("Before 1988/89", "After 1988/89", "")) + geom_text(data = label.df, label = label.df$label, size=6)
