library(coda)
library(deSolve)
library(adaptMCMC)
library(ggplot2)
library(cowplot)
library(survival)
theme_set(theme_cowplot())

setwd("C:/RData")


# Getting DEB parameters
m = readRDS("FullStarve_SR_adaptMCMC3.Rda")

params = m$samples[which.max(m$log.p),]
params["supply"] = 10 # mg C/day

signif(params, 3)

# The model definition is written in C, so we need to use Rtools and load the model as a dll
rtools <- "C:\\Rtools\\bin"
gcc <- "C:\\Rtools\\gcc-4.6.3\\bin"
path <- strsplit(Sys.getenv("PATH"), ";")[[1]]
new_path <- c(rtools, gcc, path)
new_path <- new_path[!duplicated(tolower(new_path))]
Sys.setenv(PATH = paste(new_path, collapse = ";"))

# compile my model from C definition
dyn.unload("SizeCompModel_Shrink_sim.dll") # unload dll
system("R CMD SHLIB SizeCompModel_Shrink_sim.c")
dyn.load("SizeCompModel_Shrink_sim.dll") # Load dll

#Solves all deterministic DEB skeletons
solve.DEB<-function(params, inits, duration){
  DEB.run <- lsoda(inits, seq(from=0, to=duration, by=1), func = "derivs", dllname = "SizeCompModel_Shrink_sim", 
                   initfunc = "initmod",  nout=4, outnames=c("Survival", "LG", "Survival_C", "LGC"), params, atol=1e-6,  rtol=1e-8, maxsteps=1e8)
  DEB.run
  
}


solve.DEBs<-function(pars, inits, duration){
  parms = as.numeric(pars)
  params = c(iM=parms[1], k=parms[2], M=parms[3], EM=parms[4], Fh=parms[5], muD=parms[6],
             DR=parms[7], fe=parms[8], yRP=parms[9], ph=parms[10], yPE=parms[11], iPM=parms[12],
             eh=parms[13], mP=parms[14], alpha=parms[15], yEF=parms[25], LM=parms[17],
             kd=parms[18], z=parms[19], kk=parms[20], hb=parms[21], theta=parms[22], mR=parms[23], yVE=parms[24], supply=parms[36])
  solve.DEB(params, inits, duration)
}

D_L = function(L, DR){ # estimate a reasonable amount of maturity for a given size
  a = DR/8^2
  return(min(DR, a*L^2))
}

inits = c(F = 0, L=4, e=0.8, D = D_L(4, params["DR"]), RH = 0, P = 2.85e-5, RP = 0, DAM=0, HAZ=0, LC=4, eC=0.8, DC=D_L(4, params["DR"]), RHC=0, DAMC=0, HAZC=0)
inits["LC"] = 0
solve.DEBs(params, inits, 112)

summarize.SW = function(F_size = 4, C_size = 4, supply=10, parms=params){
  initsI = c(F = 0, L=F_size, e=0.8, D = D_L(F_size, parms["DR"]), RH = 0, P = 2.85e-5, RP = 0, DAM=0, HAZ=0, LC=C_size, eC=0.8, DC=D_L(C_size, params["DR"]), RHC=0, DAMC=0, HAZC=0)

  parms["supply"] = supply; inits["F"] = supply

  predsI = solve.DEBs(parms, initsI, 112)
  if(attributes(predsI)$istate[1] != 2){print(c(parms, F_size, C_size))}

  predsI = data.frame(predsI)
  # expected parasite reproduction
  deathsI = c(predsI$Survival[1:112] - predsI$Survival[2:113], predsI$Survival[113])
  # Daily
  Cercs = sum(deathsI * predsI$RP[1:113]/(4e-5 * (1:113) ))*7 # *7 to convert from shedding per day to per week, which is form of data

  # Max gigantism
  Final_Length = sum(deathsI*predsI$LG[1:113])
  
  # Castration phenotype
  egglaying = predsI$RH[2:113] - predsI$RH[1:112]
  Inf_Eggs = max(which(egglaying > 0.015))
  
  # Expected lifespan
  median.lifespan = sum(predsI$Survival)
  
  return(c(Cercs, Final_Length, Inf_Eggs, median.lifespan))
}

summarize.SW(F_size=4,  C_size = 0, supply=10)

F_vals = 4 
C_vals = seq(from=0, to=16, by=0.25)
F.vec = rep(F_vals, times = length(C_vals))
C.vec = rep(C_vals, each = length(F_vals))


# Get a priori model predictions using the fits from Civitello et al. 2020 and the
# Experimental conditions and duration from this experiment
output = mapply(summarize.SW, F_size=F.vec, C_size=C.vec, supply=21.42/7) # Supply rate corresponds to this exp
output2 = data.frame(t(output))
colnames(output2) = c("Cercariae", "Final_Length", "I_Eggs", "Lifespan")
output2[,("F_size")] = F.vec
output2[,"C_size"] = C.vec


### Basic stats for Rachel's paper
setwd("C:/RData")
snails = read.csv("Size_comp_LT.csv")

# Data sheet contains time series for each snail, so these functions extract cumulative totals for each snail
MaxL = aggregate(Length_F ~ Snail, data=snails, FUN=max)
Lifespan = aggregate(Alive ~ Snail, data=snails, FUN=sum)
Cerc_total = aggregate(C_Cercs ~ Snail, data=snails, FUN=max)
# These are needed for the interval-censored accelerated failure time model
Time1 = aggregate(Time1 ~ Snail, data=snails, FUN=mean)
Time2 = aggregate(Time2 ~ Snail, data=snails, FUN=mean)
Event = aggregate(Event ~ Snail, data=snails, FUN=mean)

# Collect all of this processed data into a single data frame
snail_stats = data.frame("ID" = MaxL$Snail, 
                         "Infected" = rep(c("Yes", "No", "Yes", "No", "Yes", "No", "Yes", "No", "Yes", "No", "Yes", "No"), times=c(16, 5, 19, 5, 12, 5, 13, 5, 13, 5, 11, 5)), 
                         "Competitor" = rep(c(0, 2, 4, 8, 13, 16), times=c(16+5, 19+5, 12+5, 13+5, 13+5, 11+5)),  
                         "Final_Length" = MaxL$Length_F, "Lifespan" = Lifespan$Alive, 
                         "Time1" = Time1$Time1, 
                         "Time2" = Time2$Time2, 
                         "Event" = Event$Event, 
                         "Total Cercs" = Cerc_total$C_Cercs, 
                         "Rate Cercs" = Cerc_total$C_Cercs/Lifespan$Alive)

# Analysis for Panel A, effect of competitor size on cercariae per shed?
# Log transformed to deal w/ heteroskedaskticity
mA = lm(log(Rate.Cercs) ~ Competitor, data=subset(snail_stats, Infected == "Yes"))
summary(mA) # Yes

# Analysis for Panel B, effect of competitor size on final length?
mB = lm(Final_Length ~ Competitor, data=subset(snail_stats, Infected == "Yes"))
summary(mB) # Yes

# Analysis for Panel C, effect of competitor size on lifespan?
# Accelerated failure time model with Weibull distribution and interval censored observations
mC = survreg(Surv(time=Time1, time2=Time2, event=Event, type="interval") ~ Competitor, data=subset(snail_stats, Infected == "Yes"))
summary(mC) # Yes


### Collecting treatment means and SEs for ggplotting
SEM = function(x){
  sd(x)/sqrt(length(na.omit(x)))
}

# Cercs
cerc_means = aggregate(Rate.Cercs ~ Competitor, FUN=mean, data=subset(snail_stats, Infected == "Yes"))
cerc_ses = aggregate(Rate.Cercs ~ Competitor, FUN=SEM, data=subset(snail_stats, Infected == "Yes"))
cerc_means = cbind(cerc_means, "SE" = cerc_ses$Rate.Cercs)
cerc_means

# Lengths
L_means = aggregate(Final_Length ~ Competitor, FUN=mean, data=subset(snail_stats, Infected == "Yes"))
L_ses = aggregate(Final_Length ~ Competitor, FUN=SEM, data=subset(snail_stats, Infected == "Yes"))
L_means = cbind(L_means, "SE" = L_ses$Final_Length)

# Survival
S_means = aggregate(Lifespan ~ Competitor, FUN=mean, data=subset(snail_stats, Infected == "Yes"))
S_ses = aggregate(Lifespan ~ Competitor, FUN=SEM, data=subset(snail_stats, Infected == "Yes"))
S_means = cbind(S_means, "SE" = S_ses$Lifespan)
S_means[,c("Lifespan", "SE")] = S_means[,c("Lifespan", "SE")]*7

### ggplotting
pA = ggplot(data=output2, aes(x=C_size, y=Cercariae)) +
  theme(plot.margin = unit(c(0, 0.25, 0.75, 0.5), "cm")) + 
  theme(legend.position="none",
        axis.ticks.length = unit(-1.5, "mm"),
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_text(size=8, margin = margin(r=6))) + ylim(c(0, 800)) +
  geom_line() +
  geom_point(data=cerc_means, inherit.aes=F, aes(x=Competitor, y=Rate.Cercs)) + 
  geom_linerange(data=cerc_means, inherit.aes=F, aes(x=Competitor, ymin=Rate.Cercs - SE, ymax=Rate.Cercs + SE))

pB = ggplot(data=output2, aes(x=C_size, y=Final_Length)) +
  theme(plot.margin = unit(c(0, 0.25, 0.75, 0.5), "cm")) + 
  theme(legend.position="none",
        axis.ticks.length = unit(-1.5, "mm"),
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_text(size=8, margin = margin(r=6))) + ylim(c(0, 16)) +
  geom_line() +
  geom_point(data=L_means, inherit.aes=F, aes(x=Competitor, y=Final_Length)) + 
  geom_linerange(data=L_means, inherit.aes=F, aes(x=Competitor, ymin=Final_Length - SE, ymax=Final_Length + SE))

pC = ggplot(data=output2, aes(x=C_size, y=Lifespan)) +
  theme(plot.margin = unit(c(0, 0.25, 0.75, 0.5), "cm")) + 
  theme(legend.position="none",
        axis.ticks.length = unit(-1.5, "mm"),
        axis.title.x = element_blank(), axis.text.x = element_text(size=8, margin = margin(t=6)),
        axis.title.y = element_blank(), axis.text.y = element_text(size=8, margin = margin(r=6))) + ylim(c(0, 120)) +
  geom_line() +
  geom_point(data=S_means, inherit.aes=F, aes(x=Competitor, y=Lifespan)) + 
  geom_linerange(data=S_means, inherit.aes=F, aes(x=Competitor, ymin=Lifespan - SE, ymax=Lifespan + SE))



Fig1 = plot_grid(pA, pB, pC, align="v", ncol=1, nrow=3, axis="rltb", scale=1) +
  # y-axis labels
  draw_label("Daily cercarial production rate", x=0.03, y=0.85, angle=90, size=10) +
  draw_label("Length at death, mm", x=0.03, y=0.51, angle=90, size=10) +
  draw_label("Host survival, days", x=0.03, y=0.2, angle=90, size=10) +
  # x-axis labels
  draw_label("Competitor Size, mm", x = 0.55, y = 0.02, size=10) +
  # panel labels
  draw_label("A", x = 0.26, y = 0.985, size=10) + 
  draw_label("Competitor size:", x = 0.75, y = 0.95, size=8, hjust=0.5) +
  draw_label(expression(paste("P < ", 10^-12)), x = 0.75, y = 0.93, size=8, hjust=0.5) +
  draw_label(expression(paste(R^2, " = 0.46")), x = 0.75, y = 0.91, size=8, hjust=0.5) +
  draw_label("B", x = 0.26, y = 0.65, size=10) +
  draw_label("Competitor size:", x = 0.75, y = 0.45, size=8, hjust=0.5) +
  draw_label(expression(paste("P < ", 10^-12)), x = 0.75, y = 0.43, size=8, hjust=0.5) +
  draw_label(expression(paste(R^2, " = 0.46")), x = 0.75, y = 0.41, size=8, hjust=0.5) +
  draw_label("C", x = 0.26, y = 0.32, size=10) +
  draw_label("Competitor size:", x = 0.75, y = 0.16, size=8, hjust=0.5) +
  draw_label("P = 0.003", x = 0.75, y = 0.14, size=8, hjust=0.5)

save_plot("Fig1_SizeComp.png", Fig1, ncol=1, nrow=3, base_height=2, base_aspect_ratio = 1.05, dpi=600, units="in")

Fig1
