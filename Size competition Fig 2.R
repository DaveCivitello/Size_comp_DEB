library(coda)
library(deSolve)
library(adaptMCMC)

setwd("C:/RData")

thin.int = 1000
m = readRDS("Full_Fitting_Size_AR1.Rda") #$samples[seq(from=50001, to=200000, by=thin.int),c(1:25)]
params = m$samples[which.max(m$log.p),]
#signif(params, 3)

# compile my model from C definition
dyn.unload("SizeCompModel_Shrink.dll") # unload dll
# system("R CMD SHLIB SizeCompModel_Shrink.c")
dyn.load("SizeCompModel_Shrink.dll") # Load dll

# Lines up events
experiment.events = function(initial.food){
  time.F = c(sort(c((5:19)*7,((5:19)*7 - 3))), 56) # gets the fourth and seventh day of every week
  
  var.F = c(rep("F", times = length(time.F)-1),  "HAZ")
  
  value.F = c(rep(as.numeric(initial.food), times = length(time.F)-1), 0)
  
  method = rep("replace", times = length(time.F))
  out = data.frame(var = var.F, time = time.F, value = value.F, method = method)
  out[order(out$time),]
}

experiment.events(10.76)

#Solves all deterministic DEB skeletons
solve.DEB<-function(params, inits, duration, events){
  DEB.run <- as.data.frame(lsoda(inits, 28:duration, func = "derivs", dllname = "SizeCompModel_Shrink",
                                 initfunc = "initmod",  nout=4, outnames=c("Survival", "LG", "Survival_C", "LGC"), params,  
                                 rtol=1e-6, atol=1e-6, maxsteps=500000, events = list(data = events)))
  DEB.run = subset(DEB.run, select=c("time", "L", "LG", "e", "RP", "RH", "Survival", "LGC", "RHC", "Survival_C"))
  DEB.run[,"Rtotal"] = DEB.run[,"RH"] + DEB.run[,"RHC"]
  result <- rbind(DEB.run)
  result
  
}

solve.DEBs<-function(parms, duration){
  result = matrix(nrow=0, ncol=11)
  comp_size_vector = c(0, 2.3, 4.0, 8.1, 13.4, 16.9)
  comp_D_vector = as.numeric(c(0, parms["DR"]/8, parms["DR"]/4, parms["DR"], parms["DR"], parms["DR"]))
  for ( i in 1:length(comp_size_vector)){
    parms = as.numeric(parms)
    params = c(iM=parms[1], k=parms[2], M=parms[3], EM=parms[4], Fh=parms[5], muD=parms[6],
               DR=parms[7], fe=parms[8], yRP=parms[9], ph=parms[10], yPE=parms[11], iPM=parms[12],
               eh=parms[13], mP=parms[14], alpha=parms[15], yEF=parms[25], LM=parms[17],
               kd=parms[18], z=parms[19], kk=parms[20], hb=parms[21], theta=parms[22], mR=parms[23], yVE=parms[24], startage=28)
    inits = c(F = 10.78, L=3.98, e=0.9, D = as.numeric(params["DR"]/4), RH = 0, P = 2.85e-5, RP = 0, DAM=0, HAZ=0, LC=comp_size_vector[i], eC=0.9, DC=comp_D_vector[i], RHC=0, DAMC=0, HAZC=0)
    events = experiment.events(inits["F"])
    simulation = solve.DEB(params, inits, duration, events)
    simulation[which(simulation$time <= 56), "Survival"] = 1 # Conditions on individuals surviving to diagnosis
    result = rbind(result, simulation)
  }
  result = subset(result, time%%7==0)
  result
}

predsI = solve.DEBs(params, 133)
predsI["Color"] = rep(c("black", "brown", "purple", "blue", "green", "red"), each=16)


# ### Uninfected predictions
solve.DEBs<-function(parms, duration){
  result = matrix(nrow=0, ncol=11)
  comp_size_vector = c(0, 2.3, 4.0, 8.1, 13.4, 16.9)
  comp_D_vector = as.numeric(c(0, parms["DR"]/8, parms["DR"]/4, parms["DR"], parms["DR"], parms["DR"]))
  for ( i in 1:length(comp_size_vector)){
    parms = as.numeric(parms)
    params = c(iM=parms[1], k=parms[2], M=parms[3], EM=parms[4], Fh=parms[5], muD=parms[6],
               DR=parms[7], fe=parms[8], yRP=parms[9], ph=parms[10], yPE=parms[11], iPM=parms[12],
               eh=parms[13], mP=parms[14], alpha=parms[15], yEF=parms[25], LM=parms[17],
               kd=parms[18], z=parms[19], kk=parms[20], hb=parms[21], theta=parms[22], mR=parms[23], yVE=parms[24], startage=28)
    inits = c(F = 10.78, L=3.98, e=0.6, D = as.numeric(params["DR"]/4), RH = 0, P = 0, RP = 0, DAM=0, HAZ=0, LC=comp_size_vector[i], eC=0.6, DC=comp_D_vector[i], RHC=0, DAMC=0, HAZC=0)
    events = experiment.events(inits["F"])
    events = subset(events, var == "F")
    simulation = solve.DEB(params, inits, duration, events)
    #simulation[which(simulation$time == 56), "Survival"] = 1
    result = rbind(result, simulation)
  }
  result = subset(result, time%%7==0)
  result
}

predsU = solve.DEBs(params, 133)
predsU["Color"] = rep(c("black", "brown", "purple", "blue", "green", "red"), each=16)


####### ggplotting ######

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(gridExtra)
# Trying to make ggplots for Figs for starvation experiment paper

# Get infected and uninfected snails
setwd("C:/RData")
snails = read.csv("Size_comp_LT.csv")

head(snails)

# Simple function to calculate Standard Errors
SEM = function(x){
  sd(x)/sqrt(length(x))
}

# Summarize all data
snail_mean_Ls = aggregate(Length_F ~ Week*Competitor_Size*Exposure, data=snails, FUN=mean, drop=F)
snail_SE_Ls = aggregate(Length_F ~  Week*Competitor_Size*Exposure, data=snails, FUN=SEM, drop=F)

compt_mean_Ls = aggregate(Length_C ~ Week*Competitor_Size*Exposure, data=snails, FUN=mean, drop=F)
compt_SE_Ls = aggregate(Length_C ~  Week*Competitor_Size*Exposure, data=snails, FUN=SEM, drop=F)

snail_mean_Es = aggregate(log10(C_Eggs+1) ~   Week*Competitor_Size*Exposure, data=snails, FUN=mean, drop=F)
colnames(snail_mean_Es)[4] = "C_Eggs"
snail_SE_Es = aggregate(log10(C_Eggs+1) ~   Week*Competitor_Size*Exposure, data=snails, FUN=SEM, drop=F)
colnames(snail_SE_Es)[4] = "C_Eggs"

snail_summary = cbind(snail_mean_Ls, "Length_F_SE" = snail_SE_Ls$Length)
head(snail_summary)

snail_mean_Ws = aggregate(log10(1+C_Cercs) ~  Week*Competitor_Size, data=subset(snails, Exposure=="Y"), FUN=mean, drop=F)
colnames(snail_mean_Ws)[3] = "C_Worms"
#snail_mean_Ws
snail_SE_Ws = aggregate(log10(1+C_Cercs) ~  Week*Competitor_Size, data=subset(snails, Exposure=="Y"), FUN=SEM, drop=F)
colnames(snail_SE_Ws)[3] = "C_Worms"

snail_mean_Alive = aggregate(Alive ~  Week*Competitor_Size*Exposure, data=snails, FUN=mean, drop=F)

# competitor sizes
comp_Ls = aggregate(Length_C ~ Competitor_Size, data=subset(snails, Week==0), FUN=mean, drop=F)
comp_SE_Ls = aggregate(Length_C ~  Competitor_Size, data=subset(snails, Week==0), FUN=SEM, drop=F)
comp_lengths = cbind(comp_Ls, "Length_C_SE" = comp_SE_Ls$Length_C)

anova(lm(Length_C ~ Competitor_Size*Exposure, data=subset(snails, Week==0)))

#### build the summarized dataframe
snail_summary = cbind(snail_mean_Ls,"Length_F_SE" = snail_SE_Ls$Length)
snail_summary = cbind(snail_summary,"C_Eggs" = snail_mean_Es$C_Eggs)
snail_summary = cbind(snail_summary,"C_Eggs_SE" = snail_SE_Es$C_Eggs)
snail_summary = cbind(snail_summary, "Survival" = snail_mean_Alive$Alive)

snail_summary.Inf = subset(snail_summary, Exposure == "Y")
snail_summary.Inf = cbind(snail_summary.Inf, "C_Worms" = snail_mean_Ws$C_Worms, "C_Worms_SE" = snail_SE_Ws$C_Worms)
snail_summary.Inf = cbind(snail_summary.Inf, "Color" = c(rep(c("black", "brown", "purple", "blue", "green", "red"), each=16)))

snail_summary.Un = subset(snail_summary, Exposure == "N")
snail_summary.Un = cbind(snail_summary.Un, "Color" = c(rep(c("black", "brown", "purple", "blue", "green", "red"), each=16)))

# # Lining up predictions
snail_summary.Inf = cbind(snail_summary.Inf, subset(predsI, select=c("LG", "RP", "LGC", "Rtotal", "Survival")))
snail_summary.Inf[,13] = log10(1+snail_summary.Inf[,13]/4e-5)
snail_summary.Inf[,15] = log10(1+snail_summary.Inf[,15]/0.015)
colnames(snail_summary.Inf)[16] = "SurvProb"

snail_summary.Inf = cbind(snail_summary.Inf, "Length_C" = c(rep(NA, times=16), compt_mean_Ls$Length_C[81:160]))
snail_summary.Inf = cbind(snail_summary.Inf,  "Length_C_SE" = c(rep(NA, times=16), compt_SE_Ls$Length_C[81:160]))
# 
snail_summary.Un = cbind(snail_summary.Un, subset(predsU, select=c("LG", "RP", "LGC", "Rtotal", "Survival")))
snail_summary.Un[,13] = log10(1+snail_summary.Un[,13]/0.015)
colnames(snail_summary.Un)[14] = "SurvProb"
snail_summary.Un = cbind(snail_summary.Un, "Length_C" = c(rep(NA, times=16), compt_mean_Ls$Length_C[1:80]))
snail_summary.Un = cbind(snail_summary.Un, "Length_C_SE" = c(rep(NA, times=16), compt_SE_Ls$Length_C[1:80]))

# Make plots 
p1 = 
  ggplot(data=subset(snail_summary.Inf, !is.nan(Length_F) & is.finite(Length_F)), aes(x=Week, y=Length_F, group=Color, color=factor(Color))) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  +  theme(axis.text = element_text(size=10)) + 
  theme(legend.position="None", axis.ticks.length = unit(-1.5, "mm"), axis.text.x = element_blank(),
        axis.text.y = element_text(margin = margin(r=6)), axis.title.x= element_blank(), axis.title.y=element_blank()) +
  scale_y_continuous(breaks=c(5, 10, 15), limits=c(0,20)) +
  xlim(c(0, 16)) +   labs(x="", y="") + 
  geom_point(size=1) + 
  geom_linerange(aes(ymin=Length_F-Length_F_SE, ymax=Length_F + Length_F_SE)) +
  scale_color_manual(values=c("#67001f","#4393c3", "#b2182b", "#2166ac", "#f4a582",  "#053061"))

p1a = p1 +
  geom_line(aes(x=Week, y=LG, group=Color, color=Color))


p2 = ggplot(data=subset(snail_summary.Un, !is.nan(Length_F) & is.finite(Length_F)), aes(x=Week, y=Length_F, group=Color, color=factor(Color))) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  +  theme(axis.text = element_text(size=10)) + 
  theme(legend.position="None", axis.ticks.length = unit(-1.5, "mm"), axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.title.x= element_blank(), axis.title.y=element_blank()) +
  scale_y_continuous(breaks=c(0, 5, 10, 15), limits=c(0,20)) +
  xlim(c(0, 16)) +  labs(x="", y="") + 
  geom_point(size=1) + 
  geom_linerange(aes(ymin=Length_F-Length_F_SE, ymax=Length_F+ Length_F_SE)) +
  scale_color_manual(values=c("#67001f","#4393c3", "#b2182b", "#2166ac", "#f4a582",  "#053061"))

p2a = p2 +
  geom_line(aes(x=Week, y=LG, group=Color, color=Color))

p3 =   ggplot(data=subset(snail_summary.Inf, !is.nan(Length_C) & is.finite(Length_C)), aes(x=Week, y=Length_C, group=Color, color=factor(Color))) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  +  theme(axis.text = element_text(size=10)) + 
  theme(legend.position="None", axis.ticks.length = unit(-1.5, "mm"), axis.text.x = element_blank(),
        axis.text.y = element_text(margin = margin(r=6)), axis.title.x= element_blank(), axis.title.y=element_blank()) +
  scale_y_continuous(breaks=c(5, 10, 15), limits=c(0,20)) +
  xlim(c(0, 16)) +   labs(x="", y="") + 
  geom_point(size=1) + 
  geom_linerange(aes(ymin=Length_C-Length_F_SE, ymax=Length_C + Length_C_SE)) +
  scale_color_manual(values=c("#67001f","#4393c3", "#b2182b", "#2166ac", "#f4a582",  "#053061"))

p3a = p3  +
  geom_line(aes(x=Week, y=LGC, group=Color, color=Color))


p4 = ggplot(data=subset(snail_summary.Un, !is.nan(Length_C) & is.finite(Length_C)), aes(x=Week, y=Length_C, group=Color, color=factor(Color))) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  +  theme(axis.text = element_text(size=10)) + 
  theme(legend.position="None", axis.ticks.length = unit(-1.5, "mm"), axis.text.x = element_blank(),
        axis.text.y = element_text(margin = margin(r=6)), axis.title.x= element_blank(), axis.title.y=element_blank()) +
  scale_y_continuous(breaks=c(5, 10, 15), limits=c(0,20)) +
  xlim(c(0, 16)) +   labs(x="", y="") + 
  geom_point(size=1) + 
  geom_linerange(aes(ymin=Length_C-Length_F_SE, ymax=Length_C + Length_C_SE)) +
  scale_color_manual(values=c("#67001f","#4393c3", "#b2182b", "#2166ac", "#f4a582",  "#053061"))

p4a = p4 +
  geom_line(aes(x=Week, y=LGC, group=Color, color=Color))

p5 = ggplot(subset(snail_summary.Inf, Week >= 4), aes(x=Week, y=Survival, group=Color, color=factor(Color))) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))   +  theme(axis.text = element_text(size=10)) +
  scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) +
  theme(legend.position="None", axis.ticks.length = unit(-1.5, "mm"), axis.text.x = element_blank(),
        axis.text.y = element_text(margin = margin(r=6)), axis.title.x= element_blank(), axis.title.y=element_blank()) +
  xlim(c(0, 16)) + labs(x="", y="") +
  geom_point(size=1) +
  scale_color_manual(values=c("#67001f","#4393c3", "#b2182b", "#2166ac", "#f4a582",  "#053061"))

p5a = p5  +
  geom_line(aes(x=Week, y=SurvProb, group=Color, color=Color))


p6 = ggplot(snail_summary.Un, aes(x=Week, y=Survival, group=Color, color=factor(Color))) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))   +  theme(axis.text = element_text(size=10)) +
  scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0, 1)) +
  theme(legend.position="None", axis.ticks.length = unit(-1.5, "mm"), axis.ticks.y = ,
        axis.text.x = element_text(margin = margin(t=6)), axis.text.y = element_blank(), 
        axis.title.x= element_blank(), axis.title.y=element_blank()) +
  xlim(c(0, 16)) + labs(x="", y=NULL) +
  geom_point(size=1) +
  scale_color_manual(values=c("#67001f","#4393c3", "#b2182b", "#2166ac", "#f4a582",  "#053061"))

p6a = p6  +
  geom_line(aes(x=Week, y=SurvProb, group=Color, color=Color))

p7 = ggplot(data=subset(snail_summary.Inf, !is.nan(C_Worms) & is.finite(C_Worms)), aes(x=Week, y=C_Worms, group=Color, color=factor(Color))) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +  theme(axis.text = element_text(size=10)) +
  theme(legend.position="None", axis.ticks.length = unit(-1.5, "mm"), #axis.title = element_text(size=10),
        axis.text.x = element_text(margin = margin(t=6)), axis.text.y = element_text(margin = margin(r=6)),
        axis.title.x= element_blank(), axis.title.y=element_blank()) +
  xlim(c(0, 16)) +  labs(x="", y="") + 
  geom_point(size=1) + 
  geom_linerange(aes(ymin=C_Worms-C_Worms_SE, ymax=C_Worms+ C_Worms_SE)) +
  scale_color_manual(values=c("#67001f","#4393c3", "#b2182b", "#2166ac", "#f4a582",  "#053061"))

p7a = p7 +
  geom_line(aes(x=Week, y=RP, group=Color, color=Color))

spacer = ggplot(data=snail_summary.Inf, aes(x=Week, y=C_Worms)) +
  geom_blank() + theme_void()

Fig1 = plot_grid(spacer, spacer, spacer,spacer,
                 spacer, p1a, p2a, spacer,
                 spacer, p3a, p4a, spacer,
                 spacer, p5a, p6a, spacer,
                 spacer, p7a, spacer,spacer,
                 spacer, spacer, spacer,spacer,
                 align="hv", ncol=4, nrow=6, rel_widths=c(0.15, 1, 1, 0.05), rel_heights=c(0.1, 1, 1, 1, 1, 0.1), axis="rltb", scale=1) +
  # Column labels
  draw_label("Infected", x=0.35, y=0.99, vjust=1) +
  draw_label("Uninfected", x=0.8, y=0.99, vjust=1) +
  # y-axis labels
  draw_label("Host length\n(mm) ± SE", x=0.04, y=0.86, angle=90, size=12) +
  draw_label("Competitor length\n(mm) ± SE", x=0.04, y=0.63, angle=90, size=12) +
  draw_label("Host survival", x=0.04, y=0.4, angle=90, size=12) +
  draw_label("Cumulative cercarial", x=0.02, y=0.15, angle=90, size=12) +
  draw_label(expression(paste("production, ", paste(log[10], "(x+1) ± SE"))), x=0.06, y=0.15, angle=90, size=12) +
  # x-axis labels
  draw_label("Weeks post infection", x = 0.55, y = 0.01, vjust=0 ,size=12) +
  # panel labels
  draw_label("A", x = 0.15, y = 0.96) + 
  draw_label("B", x = 0.6, y = 0.96) +
  draw_label("C", x = 0.15, y = 0.73) +
  draw_label("D", x = 0.6, y = 0.73) +
  draw_label("E", x = 0.15, y = 0.465) +
  draw_label("F", x = 0.6, y = 0.465) +
  draw_label("G", x = 0.15, y = 0.25) +
  # fit statistics
  # draw_label(expression(paste(r[c], " = 0.95")), x=0.45, y=0.78, size=10)+
  # draw_label(expression(paste(r[c], " = 0.97")), x=0.92, y=0.78, size=10)+
  # draw_label(expression(paste(r[c], " = 0.71")), x=0.47, y=0.73, size=10)+
  # draw_label(expression(paste(r[c], " = 0.93")), x=0.9, y=0.73, size=10)+
  # draw_label("AUC = 0.87", x=0.25, y=0.3, size=10)+
  # draw_label("AUC = 0.77", x=0.7, y=0.3, size=10)+
  # draw_label(expression(paste(r[c], " = 0.97")), x=0.47, y=0.06, size=10)+
  # legend
  draw_label("Competitor size, mm", x=0.78, y=0.23, size=11)+
  draw_label(".", x=0.75, y=0.20, size=28, colour="#67001f")+ draw_label("None", x=0.77, y=0.20, size=10, hjust=0)+ # Alt+0149 for bullet point
  draw_label(".", x=0.75, y=0.18, size=28, colour="#b2182b")+ draw_label("2", x=0.77, y=0.18, size=10, hjust=0)+
  draw_label(".", x=0.75, y=0.16, size=28, colour="#f4a582")+ draw_label("4", x=0.77, y=0.16, size=10, hjust=0)+
  draw_label(".", x=0.75, y=0.14, size=28, colour="#4393c3")+ draw_label("8", x=0.77, y=0.14, size=10, hjust=0)+
  draw_label(".", x=0.75, y=0.12, size=28, colour="#2166ac")+ draw_label("12-15", x=0.77, y=0.12, size=10, hjust=0)+
  draw_label(".", x=0.75, y=0.10, size=28, colour="#053061")+ draw_label(">15", x=0.77, y=0.10, size=10, hjust=0)
save_plot("Fig2_SizeComp.png", Fig1, ncol=2, nrow=4, base_height=2, base_aspect_ratio = 1.1, dpi=600, units="in")
