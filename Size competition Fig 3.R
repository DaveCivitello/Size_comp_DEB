library(coda)
library(deSolve)
library(adaptMCMC)

setwd("C:/RData")

thin.int = 10000
m = readRDS("Full_Fitting_SizeComp.Rda")$samples[seq(from=50001, to=250000, by=thin.int),c(1:25)]
m2 = readRDS("Full_Fitting_SizeComp3.Rda")$samples[seq(from=50001, to=250000, by=thin.int),c(1:25)]
m3 = readRDS("Full_Fitting_SizeComp4.Rda")$samples[seq(from=50001, to=250000, by=thin.int),c(1:25)]

m = rbind(m, m2, m3)

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



#Solves all deterministic DEB skeletons
solve.DEB<-function(params, inits, duration, events){
  DEB.run <- as.data.frame(lsoda(inits, 28:duration, func = "derivs", dllname = "SizeCompModel_Shrink",
                                 initfunc = "initmod",  nout=4, outnames=c("Survival", "LG", "Survival_C", "LGC"), parms = params,  
                                 rtol=1e-6, atol=1e-6, maxsteps=5e5, events = list(data = events)))
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
               kR=parms[18], delta0=parms[19], hdelta=parms[20], hb=parms[21], theta=parms[22], mR=parms[23], yVE=parms[24], startage=28)
    inits = c(F = 10.78, L=3.98, e=0.9, D = as.numeric(params["DR"]/4), RH = 0, P = 2.85e-5, RP = 0, DAM=0, HAZ=0, LC=comp_size_vector[i], eC=0.9, DC=comp_D_vector[i], RHC=0, DAMC=0, HAZC=0)
    events = experiment.events(inits["F"])
    simulation = solve.DEB(params, inits, duration, events)
    simulation[which(simulation$time <= 56), "Survival"] = 1 # Conditions on individuals surviving to diagnosis
    result = rbind(result, simulation)
  }
  result = subset(result, time%%7==0)
  result
}


L.matrix = matrix(,nrow=0, ncol=16*6)
LC.matrix = matrix(,nrow=0, ncol=16*6)
C.matrix = matrix(,nrow=0, ncol=16*6)
E.matrix = matrix(,nrow=0, ncol=16*6)
S.matrix = matrix(,nrow=0, ncol=16*6)


for(i in 1:length(m[,1])){
  pars = m[i,]
  results = solve.DEBs(pars, 133)
  L.matrix = rbind(L.matrix, results$LG)
  LC.matrix = rbind(LC.matrix, results$LGC)
  C.matrix = rbind(C.matrix, results$RP)
  E.matrix = rbind(E.matrix, results$RH)
  S.matrix = rbind(S.matrix, results$Survival)
  print(i)
}

Hi.L = apply(L.matrix, 2, quantile, 0.995)
Med.L = apply(L.matrix, 2, quantile, probs=0.5)
Lo.L = apply(L.matrix, 2, quantile, probs=0.005)

Hi.LC = apply(LC.matrix, 2, quantile, 0.995)
Med.LC = apply(LC.matrix, 2, quantile, probs=0.5)
Lo.LC = apply(LC.matrix, 2, quantile, probs=0.005)

Hi.C = apply(C.matrix, 2, quantile, 0.995)
Med.C = apply(C.matrix, 2, quantile, probs=0.5)
Lo.C = apply(C.matrix, 2, quantile, probs=0.005)

Hi.E = apply(E.matrix, 2, quantile, 0.995)
Med.E = apply(E.matrix, 2, quantile, probs=0.5)
Lo.E = apply(E.matrix, 2, quantile, probs=0.005)

Hi.S = apply(S.matrix, 2, quantile, 0.995)
Med.S = apply(S.matrix, 2, quantile, probs=0.5)
Lo.S = apply(S.matrix, 2, quantile, probs=0.005)

output = data.frame("time" = rep((4:19)*7, times=6), "L_H" = Hi.L, "L_L" = Lo.L, "LC_H" = Hi.LC, "LC_L" = Lo.LC, "C_H"=Hi.C, "C_L"=Lo.C, "E_H"=Hi.E, "E_L"=Lo.E, "S_H"=Hi.S, "S_L"=Lo.S, "L"=Med.L, "LC"=Med.LC, "E"=Med.E, "C"=Med.C, "S"=Med.S, "Color"=rep(c("red", "green", "blue", "purple", "brown", "black"), each=16))
head(output)
write.csv(output, "SizeComp_fit_INF.csv")

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
               kR=parms[18], delta0=parms[19], hdelta=parms[20], hb=parms[21], theta=parms[22], mR=parms[23], yVE=parms[24], startage=28)
    inits = c(F = 10.78, L=3.98, e=0.9, D = as.numeric(params["DR"]/4), RH = 0, P = 0, RP = 0, DAM=0, HAZ=0, LC=comp_size_vector[i], eC=0.9, DC=comp_D_vector[i], RHC=0, DAMC=0, HAZC=0)
    events = experiment.events(inits["F"])
    events = subset(events, var == "F")
    simulation = solve.DEB(params, inits, duration, events)
    result = rbind(result, simulation)
  }
  result = subset(result, time%%7==0)
  result
}


L.matrix = matrix(,nrow=0, ncol=16*6)
LC.matrix = matrix(,nrow=0, ncol=16*6)
C.matrix = matrix(,nrow=0, ncol=16*6)
E.matrix = matrix(,nrow=0, ncol=16*6)
S.matrix = matrix(,nrow=0, ncol=16*6)


for(i in 1:length(m[,1])){
  pars = m[i,]
  results = solve.DEBs(pars, 133)
  L.matrix = rbind(L.matrix, results$LG)
  LC.matrix = rbind(LC.matrix, results$LGC)
  C.matrix = rbind(C.matrix, results$RP)
  E.matrix = rbind(E.matrix, results$RH)
  S.matrix = rbind(S.matrix, results$Survival)
  print(i)
}

Hi.L = apply(L.matrix, 2, quantile, 0.995)
Med.L = apply(L.matrix, 2, quantile, probs=0.5)
Lo.L = apply(L.matrix, 2, quantile, probs=0.005)

Hi.LC = apply(LC.matrix, 2, quantile, 0.995)
Med.LC = apply(LC.matrix, 2, quantile, probs=0.5)
Lo.LC = apply(LC.matrix, 2, quantile, probs=0.005)

Hi.C = apply(C.matrix, 2, quantile, 0.995)
Med.C = apply(C.matrix, 2, quantile, probs=0.5)
Lo.C = apply(C.matrix, 2, quantile, probs=0.005)

Hi.E = apply(E.matrix, 2, quantile, 0.995)
Med.E = apply(E.matrix, 2, quantile, probs=0.5)
Lo.E = apply(E.matrix, 2, quantile, probs=0.005)

Hi.S = apply(S.matrix, 2, quantile, 0.995)
Med.S = apply(S.matrix, 2, quantile, probs=0.5)
Lo.S = apply(S.matrix, 2, quantile, probs=0.005)

output2 = data.frame("time" = rep((4:19)*7, times=6), "L_H" = Hi.L, "L_L" = Lo.L, "LC_H" = Hi.LC, "LC_L" = Lo.LC, "C_H"=Hi.C, "C_L"=Lo.C, "E_H"=Hi.E, "E_L"=Lo.E, "S_H"=Hi.S, "S_L"=Lo.S, "L"=Med.L, "LC"=Med.LC, "E"=Med.E, "C"=Med.C, "S"=Med.S, "Color"=rep(c("red", "green", "blue", "purple", "brown", "black"), each=16))
head(output2)
write.csv(output2, "SizeComp_fit_UN.csv")

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
  sd(x)/sqrt(na.omit(length(x)))
}

# Summarize all data
snail_mean_Ls = aggregate(Length_F ~ Week*Competitor_Size*Exposure, data=snails, FUN=mean, drop=F)
snail_SE_Ls = aggregate(Length_F ~  Week*Competitor_Size*Exposure, data=snails, FUN=SEM, drop=F)

# Have to deal with missing competitor lengths for snails with no competitors
compt_mean_Ls = aggregate(Length_C ~ Week*Competitor_Size*Exposure, data=snails, FUN=mean, drop=F)[,"Length_C"]
compt_mean_Ls = c(rep(NA, times=16), compt_mean_Ls[1:80], rep(NA, times=16), compt_mean_Ls[81:160])
compt_SE_Ls = aggregate(Length_C ~  Week*Competitor_Size*Exposure, data=snails, FUN=SEM, drop=F)[,"Length_C"]
compt_SE_Ls = c(rep(NA, times=16), compt_SE_Ls[1:80], rep(NA, times=16), compt_SE_Ls[81:160])

snail_mean_Es = aggregate(log10(C_Eggs+1) ~   Week*Competitor_Size*Exposure, data=snails, FUN=mean, drop=F)
colnames(snail_mean_Es)[4] = "C_Eggs"
snail_SE_Es = aggregate(log10(C_Eggs+1) ~   Week*Competitor_Size*Exposure, data=snails, FUN=SEM, drop=F)
colnames(snail_SE_Es)[4] = "C_Eggs"

snail_summary = cbind(snail_mean_Ls, "Length_F_SE" = snail_SE_Ls$Length, "Length_C" = compt_mean_Ls, "Length_C_SE" = compt_SE_Ls)
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

#### build the summarized dataframe
snail_summary = cbind(snail_summary,"C_Eggs" = snail_mean_Es$C_Eggs)
snail_summary = cbind(snail_summary,"C_Eggs_SE" = snail_SE_Es$C_Eggs)
snail_summary = cbind(snail_summary, "Survival" = snail_mean_Alive$Alive)

snail_summary.Inf = subset(snail_summary, Exposure == "Y")
snail_summary.Inf = cbind(snail_summary.Inf, "C_Worms" = snail_mean_Ws$C_Worms, "C_Worms_SE" = snail_SE_Ws$C_Worms)
snail_summary.Inf = cbind(snail_summary.Inf, "Color" = c(rep(c("black", "brown", "purple", "blue", "green", "red"), each=16)))

snail_summary.Un = subset(snail_summary, Exposure == "N")
snail_summary.Un = cbind(snail_summary.Un, "Color" = c(rep(c("black", "brown", "purple", "blue", "green", "red"), each=16)))

# Lining up predictions
predictions = read.csv("SizeComp_fit_INF.csv")
snail_summary.Inf = cbind(snail_summary.Inf, subset(predictions, select=c("L_H", "L_L", "LC_H", "LC_L", "C_H", "C_L", "E_H", "E_L", "S_H", "S_L", "L", "LC", "C", "E", "S")))
snail_summary.Inf[,c(18:19, 26)] = log10(1+snail_summary.Inf[,c(18:19, 26)]/4e-5)
snail_summary.Inf[,c(20:21, 27)] = log10(1+snail_summary.Inf[,c(20:21, 27)]/0.015)

predictions = read.csv("SizeComp_fit_UN.csv")
snail_summary.Un = cbind(snail_summary.Un, subset(predictions, select=c("L_H", "L_L", "LC_H", "LC_L", "E_H", "E_L", "S_H", "S_L", "L", "LC", "E", "S")))
snail_summary.Un[,c(16:17, 22)] = log10(1+snail_summary.Un[,c(16:17, 22)]/0.015)


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
  scale_color_manual(values=c("#67001f","#4393c3", "#b2182b", "#2166ac", "#f4a582",  "#053061")) +
  geom_ribbon(aes(x=Week, ymin=L_L, ymax=L_H, group=Color, color=NULL, fill=factor(Color)), alpha=0.3)  +
  geom_line(aes(x=Week, y=L, group=Color, color=Color)) +
  scale_fill_manual(values=c("#67001f","#4393c3", "#b2182b", "#2166ac", "#f4a582",  "#053061"), name="fill")


p2 = ggplot(data=subset(snail_summary.Un, !is.nan(Length_F) & is.finite(Length_F)), aes(x=Week, y=Length_F, group=Color, color=factor(Color))) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  +  theme(axis.text = element_text(size=10)) + 
  theme(legend.position="None", axis.ticks.length = unit(-1.5, "mm"), axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.title.x= element_blank(), axis.title.y=element_blank()) +
  scale_y_continuous(breaks=c(0, 5, 10, 15), limits=c(0,20)) +
  xlim(c(0, 16)) +  labs(x="", y="") + 
  geom_point(size=1) + 
  geom_linerange(aes(ymin=Length_F-Length_F_SE, ymax=Length_F + Length_F_SE)) +
  scale_color_manual(values=c("#67001f","#4393c3", "#b2182b", "#2166ac", "#f4a582",  "#053061")) +
  geom_ribbon(aes(x=Week, ymin=L_L, ymax=L_H, group=Color, color=NULL, fill=factor(Color)), alpha=0.3)  +
  geom_line(aes(x=Week, y=L, group=Color, color=Color)) +
  scale_fill_manual(values=c("#67001f","#4393c3", "#b2182b", "#2166ac", "#f4a582",  "#053061"), name="fill")

p3 =   ggplot(data=subset(snail_summary.Inf, !is.nan(Length_C) & is.finite(Length_C)), aes(x=Week, y=Length_C, group=Color, color=factor(Color))) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  +  theme(axis.text = element_text(size=10)) + 
  theme(legend.position="None", axis.ticks.length = unit(-1.5, "mm"), axis.text.x = element_blank(),
        axis.text.y = element_text(margin = margin(r=6)), axis.title.x= element_blank(), axis.title.y=element_blank()) +
  scale_y_continuous(breaks=c(5, 10, 15), limits=c(0,20)) +
  xlim(c(0, 16)) +   labs(x="", y="") + 
  geom_point(size=1) + 
  geom_linerange(aes(ymin=Length_C-Length_C_SE, ymax=Length_C + Length_C_SE)) +
  scale_color_manual(values=c("#4393c3", "#b2182b", "#2166ac", "#f4a582",  "#053061")) +
  geom_ribbon(aes(x=Week, ymin=LC_L, ymax=LC_H, group=Color, color=NULL, fill=factor(Color)), alpha=0.3)  +
  geom_line(aes(x=Week, y=LC, group=Color, color=Color)) +
  scale_fill_manual(values=c("#4393c3", "#b2182b", "#2166ac", "#f4a582",  "#053061"), name="fill")

p4 = ggplot(data=subset(snail_summary.Un, !is.nan(Length_C) & is.finite(Length_C)), aes(x=Week, y=Length_C, group=Color, color=factor(Color))) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  +  theme(axis.text = element_text(size=10)) + 
  theme(legend.position="None", axis.ticks.length = unit(-1.5, "mm"), axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.title.x= element_blank(), axis.title.y=element_blank()) +
  scale_y_continuous(breaks=c(5, 10, 15), limits=c(0,20)) +
  xlim(c(0, 16)) +   labs(x="", y="") + 
  geom_point(size=1) + 
  geom_linerange(aes(ymin=Length_C-Length_C_SE, ymax=Length_C + Length_C_SE)) +
  scale_color_manual(values=c("#4393c3", "#b2182b", "#2166ac", "#f4a582",  "#053061")) +
  geom_ribbon(aes(x=Week, ymin=LC_L, ymax=LC_H, group=Color, color=NULL, fill=factor(Color)), alpha=0.3)  +
  geom_line(aes(x=Week, y=LC, group=Color, color=Color)) +
  scale_fill_manual(values=c("#4393c3", "#b2182b", "#2166ac", "#f4a582",  "#053061"), name="fill")

p5 = ggplot(subset(snail_summary.Inf, Week >= 4), aes(x=Week, y=Survival, group=Color, color=factor(Color))) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))   +  theme(axis.text = element_text(size=10)) +
  scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) +
  theme(legend.position="None", axis.ticks.length = unit(-1.5, "mm"), axis.text.x = element_blank(),
        axis.text.y = element_text(margin = margin(r=6)), axis.title.x= element_blank(), axis.title.y=element_blank()) +
  xlim(c(0, 16)) + labs(x="", y="") +
  geom_point(size=1) +
  scale_color_manual(values=c("#67001f","#4393c3", "#b2182b", "#2166ac", "#f4a582",  "#053061"))+
  geom_ribbon(aes(x=Week, ymin=S_L, ymax=S_H, group=Color, color=NULL, fill=factor(Color)), alpha=0.3)  +
  geom_line(aes(x=Week, y=S, group=Color, color=Color)) +
  scale_fill_manual(values=c("#67001f","#4393c3", "#b2182b", "#2166ac", "#f4a582",  "#053061"), name="fill")

p6 = ggplot(snail_summary.Un, aes(x=Week, y=Survival, group=Color, color=factor(Color))) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))   +  theme(axis.text = element_text(size=10)) +
  scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0, 1)) +
  theme(legend.position="None", axis.ticks.length = unit(-1.5, "mm"), axis.ticks.y = ,
        axis.text.x = element_text(margin = margin(t=6)), axis.text.y = element_blank(), 
        axis.title.x= element_blank(), axis.title.y=element_blank()) +
  xlim(c(0, 16)) + labs(x="", y=NULL) +
  geom_point(size=1) +
  scale_color_manual(values=c("#67001f","#4393c3", "#b2182b", "#2166ac", "#f4a582",  "#053061"))+
  geom_ribbon(aes(x=Week, ymin=S_L, ymax=S_H, group=Color, color=NULL, fill=factor(Color)), alpha=0.3)  +
  geom_line(aes(x=Week, y=S, group=Color, color=Color)) +
  scale_fill_manual(values=c("#67001f","#4393c3", "#b2182b", "#2166ac", "#f4a582",  "#053061"), name="fill")

p7 = ggplot(data=subset(snail_summary.Inf, !is.nan(C_Worms) & is.finite(C_Worms)), aes(x=Week, y=C_Worms, group=Color, color=factor(Color))) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +  theme(axis.text = element_text(size=10)) +
  theme(legend.position="None", axis.ticks.length = unit(-1.5, "mm"), #axis.title = element_text(size=10),
        axis.text.x = element_text(margin = margin(t=6)), axis.text.y = element_text(margin = margin(r=6)),
        axis.title.x= element_blank(), axis.title.y=element_blank()) +
  xlim(c(0, 16)) +  labs(x="", y="") + 
  geom_point(size=1) + 
  geom_linerange(aes(ymin=C_Worms-C_Worms_SE, ymax=C_Worms+ C_Worms_SE)) + 
  scale_color_manual(values=c("#67001f","#4393c3", "#b2182b", "#2166ac", "#f4a582",  "#053061")) +
  geom_ribbon(aes(x=Week, ymin=C_L, ymax=C_H, group=Color, color=NULL, fill=factor(Color)), alpha=0.3)  +
  geom_line(aes(x=Week, y=C, group=Color, color=Color)) +
  scale_fill_manual(values=c("#67001f","#4393c3", "#b2182b", "#2166ac", "#f4a582",  "#053061"), name="fill")

spacer = ggplot(data=snail_summary.Inf, aes(x=Week, y=C_Worms)) +
  geom_blank() + theme_void()

Fig1 = plot_grid(spacer, spacer, spacer,spacer,
                 spacer, p1, p2, spacer,
                 spacer, p3, p4, spacer,
                 spacer, p5, p6, spacer,
                 spacer, p7, spacer,spacer,
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
  draw_label(expression(paste(r[c], " = 0.91")), x=0.45, y=0.78, size=10)+
  draw_label(expression(paste(r[c], " = 0.92")), x=0.92, y=0.78, size=10)+
  draw_label(expression(paste(r[c], " = 0.96")), x=0.47, y=0.73, size=10)+
  draw_label(expression(paste(r[c], " = 0.91")), x=0.9, y=0.73, size=10)+
  draw_label("AUC = 0.89", x=0.25, y=0.3, size=10)+
  draw_label("AUC = 0.64", x=0.7, y=0.3, size=10)+
  draw_label(expression(paste(r[c], " = 0.96")), x=0.47, y=0.06, size=10)+
  # legend
  draw_label("Competitor size, mm", x=0.78, y=0.23, size=11)+
  draw_label("•", x=0.75, y=0.20, size=28, colour="#67001f")+ draw_label("None", x=0.77, y=0.20, size=10, hjust=0)+ # Alt+0149 for bullet point
  draw_label("•", x=0.75, y=0.18, size=28, colour="#b2182b")+ draw_label("2", x=0.77, y=0.18, size=10, hjust=0)+
  draw_label("•", x=0.75, y=0.16, size=28, colour="#f4a582")+ draw_label("4", x=0.77, y=0.16, size=10, hjust=0)+
  draw_label("•", x=0.75, y=0.14, size=28, colour="#4393c3")+ draw_label("8", x=0.77, y=0.14, size=10, hjust=0)+
  draw_label("•", x=0.75, y=0.12, size=28, colour="#2166ac")+ draw_label("12-15", x=0.77, y=0.12, size=10, hjust=0)+
  draw_label("•", x=0.75, y=0.10, size=28, colour="#053061")+ draw_label(">15", x=0.77, y=0.10, size=10, hjust=0)
save_plot("Fig3_SizeComp.png", Fig1, ncol=2, nrow=4, base_height=2, base_aspect_ratio = 1.1, dpi=600, units="in")

