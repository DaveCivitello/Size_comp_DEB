 ### Adaptive MCMC to fit "shrinking and regression" model simultaneously to resource supply and periodic 
### starvation experiments using Biomphalaria glabrata and Schistosoma mansoni

library("adaptMCMC") # needed for MCMC
library("deSolve") # needed to simulate the models
library("mvtnorm")
library("LaplacesDemon")
library("coda")

# The model definition is written in C, so we need to use Rtools and load the model as a dll
rtools <- "C:\\Rtools\\bin"
gcc <- "C:\\Rtools\\gcc-4.6.3\\bin"
path <- strsplit(Sys.getenv("PATH"), ";")[[1]]
new_path <- c(rtools, gcc, path)
new_path <- new_path[!duplicated(tolower(new_path))]
Sys.setenv(PATH = paste(new_path, collapse = ";"))

# compile my model from C definition
setwd("C:/RData")
dyn.unload("IndividualModel_Shrink_P_dilute.dll") # unload dll
system("R CMD SHLIB IndividualModel_Shrink_P_dilute.c")
dyn.load("IndividualModel_Shrink_P_dilute.dll") # Load dll

dyn.unload("SizeCompModel_Shrink.dll") # unload dll
system("R CMD SHLIB SizeCompModel_Shrink.c")
dyn.load("SizeCompModel_Shrink.dll") # Load dll


#### Fixed information ####

# 1 - Data
setwd("C:/RData")

data = read.csv("ResourceSupplyExp.csv")
data2 = read.csv("PeriodicStarveExp.csv")
data3 = read.csv("Size_comp_LT.csv")

data = list(t = data$Date, L = data$Length, Negg = data$C_Eggs, Nworms = data$C_Worms, Alive=data$Last_Alive,
            L2 = data2$Length, Negg2 = data2$C_Eggs, Nworms2 = data2$C_Worms, Alive2=data2$Last_Alive,
            L3 = data3$Length_F, L3C = data3$Length_C, Negg3 = data3$C_Eggs, Nworms3 = data3$C_Cercs, Alive3=data3$F_last_alive)


# 2 - Initial conditions
setinits.Food<-function(F0 = 16.5, L0=0.7, e0=0.9, D0 = 0, RH0 = 0, P0 = 0, RP0 = 0, DAM0=0, HAZ0=0){
  inits<-c(F = F0, L = L0, e = e0, D = D0, RH = RH0, P = P0, RP = RP0, DAM=DAM0, HAZ=HAZ0)
  return(inits)
}

setinits.Starve<-function(F0 = 10.74, L0=0.7, e0=0.9, D0 = 0, RH0 = 0, P0 = 0, RP0 = 0, DAM0=0, HAZ0=0){
  inits<-c(F = F0, L = L0, e = e0, D = D0, RH = RH0, P = P0, RP = RP0, DAM=DAM0, HAZ=HAZ0)
  return(inits)
}

setinits.Size<-function(F0 = 10.76, L0=3.98, e0=0.9, D0 = 0, RH0 = 0, P0 = 2.85e-5, RP0 = 00, DAM0=0, HAZ0=0,
                        LC0=0, eC0=0.9, DC0=0, RHC0=0, DAMC0=0, HAZC0=0){
  inits<-c(F = F0, L = L0, e = e0, D = D0, RH = RH0, P = P0, RP = RP0, DAM=DAM0, HAZ=HAZ0,
           LC = LC0, eC = eC0, DC = DC0, RHC=RHC0, DAMC = DAMC0, HAZC = HAZC0)
  return(inits)
}


# 3 - Functions for feeding events (these match the timing of feeding events in the actual ecxperiments)

# For the size competition experiment
size.events = function(initial.food){
  time.F = c(sort(c((5:19)*7,((5:19)*7 - 3))), 56) # gets the fourth and seventh day of every week
  
  var.F = c(rep("F", times = length(time.F)-1),  "HAZ")
  
  value.F = c(rep(as.numeric(initial.food), times = length(time.F)-1), 0)
  
  method = rep("replace", times = length(time.F))
  out = data.frame(var = var.F, time = time.F, value = value.F, method = method)
  out[order(out$time),]
}

# For the starvation experiment
periodic.starvation.events = function(initial.food, starve.period, infection.date, weeks.duration){
  infection.value = 2.85e-5
  all.fed.dates = infection.date + c(1,4) # Snails were infected on a Thursday, All fed on Friday and Monday, then treatments
  potential.feedings = infection.date + 5 + sort(c((1:weeks.duration)*7 - 6,((1:weeks.duration)*7 - 4),((1:weeks.duration)*7 - 2)))
  #print(0:(weeks.duration/starve.period - 1)*6 + 1)
  potential.feedings.vals = rep(initial.food, times = length(potential.feedings))
  # work out starvation
  if( starve.period == 2){
    potential.feedings.vals[(potential.feedings - (infection.date + 5))  %% 14 <= 7] = 0
  }
  if (starve.period == 3){
    potential.feedings.vals[(potential.feedings - (infection.date + 5)) %% 21 <= 14] = 0
  }
  if (starve.period == 4){
    potential.feedings.vals[(potential.feedings - (infection.date + 5)) %% 28 <= 21] = 0
  }
  event.dates = c(infection.date, all.fed.dates, potential.feedings, infection.date)
  event.values = c(infection.value, initial.food, initial.food, potential.feedings.vals, 0)
  methods = rep("replace", times=length(event.dates))
  vars = c("P", rep("F", times = length(event.dates)-2), "HAZ")
  data.frame(var=vars, time= event.dates, value = event.values, method= methods)
}

# For the resource supply experiment
feeding.events <- function(dates, var="F", value, method="replace", Infected=0){
  # Assemble all data
  events = length(dates) # number of events
  vars = rep(var, times = events)
  values = rep(value, times = events)
  methods = rep(method, times = events)
  
  #build data.frame
  result = data.frame(var = vars, time = dates, value = values, method = methods)
  if(Infected > 0){
    result = rbind(result, data.frame(var="P", time=Infected, value=2.85e-5, method="replace"))
  }
  result = rbind(result, data.frame(var="HAZ", time=63, value=0, method="replace"))
  result
}

dur.S = 133 # This experiment starts on t = 28, ends on t = 133
dur.R = 245 # Starts on t = 0
dur.P = 140 # Starts on t = 0

in.S = setinits.Size()
in.R = setinits.Food()
in.P = setinits.Starve()

Feed.R = feeding.events(dates = sort(c((1:35)*7,((1:35)*7 - 3))), var="F", in.R[1], method="replace", Infected=28)
Feed.S = size.events(10.76)

# 4 - Functions to solve DEBs
solve.DEB.food<-function(params, inits, duration, feeding.events){
  feed.sup <- feeding.events
  feed.sup.U <- subset(feed.sup, var != "P")
  parms = as.numeric(params[1:24])
  params = c(iM=parms[1], k=parms[2], M=parms[3], EM=parms[4], Fh=parms[5], muD=parms[6],
             DR=parms[7], fe=parms[8], yRP=parms[9], ph=parms[10], yPE=parms[11], iPM=parms[12],
             eh=parms[13], mP=parms[14], alpha=parms[15], yEF=parms[16], LM=parms[17],kR=parms[18], 
             delta0=parms[19], hdelta=parms[20], hb=parms[21], theta=parms[22], mR=parms[23], yVE=parms[24], startage=28)
  
  capture.output(Sup.6 <- data.frame(lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                                           initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                                           params[1:25],  rtol=1e-6, atol=1e-6,   
                                           events = list(data = feed.sup))))
  
  capture.output(Sup.6U <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                                 initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                                 params[1:25],  rtol=1e-6, atol=1e-6,   
                                 events = list(data = feed.sup.U)))
  if(attributes(Sup.6U)$istate[1] != 2)(return(Sup.6))
  
  feed.sup[1:70,3] <- 11
  feed.sup.U[1:70,3] <- 11
  capture.output(Sup.5 <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                                initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                                params[1:25],  rtol=1e-6, atol=1e-6,   
                                events = list(data = feed.sup)))
  if(attributes(Sup.5)$istate[1] != 2)(return(Sup.6))
  capture.output(Sup.5U <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                                 initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                                 params[1:25],  rtol=1e-6, atol=1e-6,   
                                 events = list(data = feed.sup.U)))
  if(attributes(Sup.5U)$istate[1] != 2)(return(Sup.6))
  
  feed.sup[1:70,3] <- 5.5
  feed.sup.U[1:70,3] <- 5.5
  capture.output(Sup.4 <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                                initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                                params[1:25],  rtol=1e-6, atol=1e-6,   
                                events = list(data = feed.sup)))  
  if(attributes(Sup.4)$istate[1] != 2)(return(Sup.6))
  capture.output(Sup.4U <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                                 initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                                 params[1:25],  rtol=1e-6, atol=1e-6,   
                                 events = list(data = feed.sup.U)))
  if(attributes(Sup.4U)$istate[1] != 2)(return(Sup.6))
  
  feed.sup[1:70,3] <- 2.75
  feed.sup.U[1:70,3] <- 2.75
  capture.output(Sup.3 <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                                initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                                params[1:25],  rtol=1e-6, atol=1e-6,   
                                events = list(data = feed.sup)))
  if(attributes(Sup.3)$istate[1] != 2)(return(Sup.6))
  capture.output(Sup.3U <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                                 initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                                 params[1:25],  rtol=1e-6, atol=1e-6,   
                                 events = list(data = feed.sup.U)))
  if(attributes(Sup.3U)$istate[1] != 2)(return(Sup.6))
  
  feed.sup[1:70,3] <- 1.375
  feed.sup.U[1:70,3] <- 1.375
  capture.output(Sup.2 <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                                initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                                params[1:25],  rtol=1e-6, atol=1e-6,   
                                events = list(data = feed.sup)))
  if(attributes(Sup.2)$istate[1] != 2)(return(Sup.6))
  capture.output(Sup.2U <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                                 initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                                 params[1:25],  rtol=1e-6, atol=1e-6,   
                                 events = list(data = feed.sup.U)))
  if(attributes(Sup.2U)$istate[1] != 2)(return(Sup.6))
  
  feed.sup[1:70,3] <- 0.6875
  feed.sup.U[1:70,3] <- 0.6875
  capture.output(Sup.1 <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                                initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                                params[1:25],  rtol=1e-6, atol=1e-6,   
                                events = list(data = feed.sup)))
  if(attributes(Sup.1)$istate[1] != 2)(return(Sup.6))
  capture.output(Sup.1U <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                                 initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                                 params[1:25],  rtol=1e-6, atol=1e-6,   
                                 events = list(data = feed.sup.U)))
  if(attributes(Sup.1U)$istate[1] != 2)(return(Sup.6))
  
  result <- rbind(#Infecteds (n=66)
    Sup.1, Sup.1, Sup.1, Sup.1, Sup.1, Sup.1, Sup.1, Sup.1, Sup.1, Sup.1, 
    Sup.2, Sup.2, Sup.2, Sup.2, Sup.2, Sup.2, Sup.2, Sup.2, Sup.2, Sup.2, Sup.2, Sup.2, 
    Sup.3, Sup.3, Sup.3, Sup.3, Sup.3, Sup.3, Sup.3, Sup.3, Sup.3, Sup.3, Sup.3, Sup.3, 
    Sup.4, Sup.4, Sup.4, Sup.4, Sup.4, Sup.4, Sup.4, Sup.4, Sup.4, Sup.4, 
    Sup.5, Sup.5, Sup.5, Sup.5, Sup.5, Sup.5, Sup.5, Sup.5, Sup.5, Sup.5, 
    Sup.6, Sup.6, Sup.6, Sup.6, Sup.6, Sup.6, Sup.6, Sup.6, Sup.6, Sup.6, Sup.6, Sup.6, 
    #Uninfecteds (n=30)
    Sup.1U, Sup.1U, Sup.1U, Sup.1U, Sup.1U, 
    Sup.2U, Sup.2U, Sup.2U, Sup.2U, Sup.2U, 
    Sup.3U, Sup.3U, Sup.3U, Sup.3U, Sup.3U, 
    Sup.4U, Sup.4U, Sup.4U, Sup.4U, Sup.4U, 
    Sup.5U, Sup.5U, Sup.5U, Sup.5U, Sup.5U,
    Sup.6U, Sup.6U, Sup.6U, Sup.6U, Sup.6U)
  
  result
  
}

# This function is customized to my model and data
solve.DEB.starve<-function(params, inits, duration){
  # Collect the params the way C likes them
  parms = as.numeric(params[1:25])
  params = c(iM=parms[1], k=parms[2], M=parms[3], EM=parms[4], Fh=parms[5], muD=parms[6],
             DR=parms[7], fe=parms[8], yRP=parms[9], ph=parms[10], yPE=parms[11], iPM=parms[12],
             eh=parms[13], mP=parms[14], alpha=parms[15], yEF=parms[25], LM=parms[17],kR=parms[18], 
             delta0=parms[19], hdelta=parms[20], hb=parms[21], theta=parms[22], mR=parms[23], yVE=parms[24], startage=14)
  
  # Simulate dynamics for 1-0
  food_1_0 = periodic.starvation.events(2.69, 0, 14, 18)
  food_1_0.U <- subset(food_1_0, var != "P")
  capture.output(out_1_0 <- data.frame(ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                                           initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), method="lsoda",
                                           params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=5e5,
                                           events = list(data = food_1_0))))
  capture.output(out_1_0U <- ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                                 initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), method="lsoda",
                                 params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=5e5,
                                 events = list(data = food_1_0.U)))
  if(attributes(out_1_0U)$istate[1] != 2)(return(out_1_0))
  
  # Simulate dynamics for 2-2
  food_2_2 = periodic.starvation.events(5.37, 2, 14, 18)
  food_2_2.U <- subset(food_2_2, var != "P")
  capture.output(out_2_2 <- ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                                initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), method="lsoda",
                                params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=5e5,
                                events = list(data = food_2_2)))
  if(attributes(out_2_2)$istate[1] != 2)(return(out_1_0))
  capture.output(out_2_2U <- ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                                 initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), method="lsoda",
                                 params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=5e5,
                                 events = list(data = food_2_2.U)))
  if(attributes(out_2_2U)$istate[1] != 2)(return(out_1_0))
  
  # Simulate dynamics for 3-3
  food_3_3 = periodic.starvation.events(8.06, 3, 14, 18)
  food_3_3.U <- subset(food_3_3, var != "P")
  capture.output(out_3_3 <- ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                                initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), method="lsoda",
                                params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=5e5,
                                events = list(data = food_3_3)))
  if(attributes(out_3_3)$istate[1] != 2)(return(out_1_0))
  capture.output(out_3_3U <- ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                                 initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), method="lsoda",
                                 params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=5e5,
                                 events = list(data = food_3_3.U)))
  if(attributes(out_3_3U)$istate[1] != 2)(return(out_1_0))
  
  # Simulate dynamics for 4-4
  food_4_4 = periodic.starvation.events(10.74, 4, 14, 18)
  food_4_4.U <- subset(food_4_4, var != "P")
  capture.output(out_4_4 <- ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                                initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), method="lsoda",
                                params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=5e5,
                                events = list(data = food_4_4)))
  if(attributes(out_4_4)$istate[1] != 2)(return(out_1_0))
  capture.output(out_4_4U <- ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                                 initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), method="lsoda",
                                 params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=5e5,
                                 events = list(data = food_4_4.U)))
  if(attributes(out_4_4U)$istate[1] != 2)(return(out_1_0))
  
  
  result <- rbind(
    out_1_0, out_1_0, out_1_0, out_1_0, out_1_0, # 5
    out_2_2, out_2_2, out_2_2, # 3
    out_3_3, out_3_3, out_3_3, out_3_3, out_3_3, out_3_3, out_3_3, out_3_3, # 8    
    out_4_4, out_4_4, out_4_4, out_4_4,     
    
    out_1_0U, out_1_0U, out_1_0U, out_1_0U, out_1_0U,
    out_2_2U, out_2_2U, out_2_2U, out_2_2U, out_2_2U,
    out_3_3U, out_3_3U, out_3_3U, out_3_3U, out_3_3U,
    out_4_4U, out_4_4U, out_4_4U, out_4_4U, out_4_4U)
  
  result
  
}

solve.DEB.Size<-function(params, inits, duration, feeding.events){
  feed.size <- feeding.events
  
  inits["D"] = params["DR"]/4
  
  initsU = inits; initsU["P"] = 0
  feed.sizeU = subset(feed.size, var == "F")
  
  # Collect the params the way C likes them
  parms = as.numeric(params[1:25])
  params = c(iM=parms[1], k=parms[2], M=parms[3], EM=parms[4], Fh=parms[5], muD=parms[6],
             DR=parms[7], fe=parms[8], yRP=parms[9], ph=parms[10], yPE=parms[11], iPM=parms[12],
             eh=parms[13], mP=parms[14], alpha=parms[15], yEF=parms[25], LM=parms[17],kR=parms[18], 
             delta0=parms[19], hdelta=parms[20], hb=parms[21], theta=parms[22], mR=parms[23], yVE=parms[24])
  
  capture.output(Size.0 <- data.frame(lsoda(inits, 28:duration, func = "derivs", dllname = "SizeCompModel_Shrink", 
                             initfunc = "initmod",  nout=4, outnames=c("Survival", "LG", "Survival_C", "LGC"), 
                             maxsteps=5e5, params[1:25],  rtol=1e-6, atol=1e-6, events = list(data = feed.size))))

  capture.output(Size.0U <- lsoda(initsU, 28:duration, func = "derivs", dllname = "SizeCompModel_Shrink", 
                   initfunc = "initmod",  nout=4, outnames=c("Survival", "LG", "Survival_C", "LGC"), 
                   maxsteps=5e5, params[1:25],  rtol=1e-6, atol=1e-6, events = list(data = feed.sizeU)))
  if(attributes(Size.0U)$istate[1] != 2)(return(Size.0))
  
  inits["LC"] = 2.3; inits["DC"] = params["DR"]/8
  initsU["LC"] = 2.3; initsU["DC"] = params["DR"]/8
  capture.output(Size.2 <- lsoda(inits, 28:duration, func = "derivs", dllname = "SizeCompModel_Shrink", 
                  initfunc = "initmod",  nout=4, outnames=c("Survival", "LG", "Survival_C", "LGC"), 
                  maxsteps=5e5, params[1:25],  rtol=1e-6, atol=1e-6, events = list(data = feed.size)))
  if(attributes(Size.2)$istate[1] != 2)(return(Size.0))
  
  capture.output(Size.2U <- lsoda(initsU, 28:duration, func = "derivs", dllname = "SizeCompModel_Shrink", 
                   initfunc = "initmod",  nout=4, outnames=c("Survival", "LG", "Survival_C", "LGC"), 
                   maxsteps=5e5, params[1:25],  rtol=1e-6, atol=1e-6, events = list(data = feed.sizeU)))
  if(attributes(Size.2U)$istate[1] != 2)(return(Size.0))
  
  inits["LC"] = 4.0; inits["DC"] = params["DR"]/4
  initsU["LC"] = 4.0; initsU["DC"] = params["DR"]/4
  capture.output(Size.4 <- lsoda(inits, 28:duration, func = "derivs", dllname = "SizeCompModel_Shrink", 
                  initfunc = "initmod",  nout=4, outnames=c("Survival", "LG", "Survival_C", "LGC"), 
                  maxsteps=5e5, params[1:25],  rtol=1e-6, atol=1e-6, events = list(data = feed.size)))
  if(attributes(Size.4)$istate[1] != 2)(return(Size.0))
  
  capture.output(Size.4U <- lsoda(initsU, 28:duration, func = "derivs", dllname = "SizeCompModel_Shrink", 
                   initfunc = "initmod",  nout=4, outnames=c("Survival", "LG", "Survival_C", "LGC"), 
                   maxsteps=5e5, params[1:25],  rtol=1e-6, atol=1e-6, events = list(data = feed.sizeU)))
  if(attributes(Size.4U)$istate[1] != 2)(return(Size.0))
  
  inits["LC"] = 8.1; inits["DC"] = params["DR"]
  initsU["LC"] = 8.1; initsU["DC"] = params["DR"]
  capture.output(Size.8 <- lsoda(inits, 28:duration, func = "derivs", dllname = "SizeCompModel_Shrink", 
                  initfunc = "initmod",  nout=4, outnames=c("Survival", "LG", "Survival_C", "LGC"), 
                  maxsteps=5e5, params[1:25],  rtol=1e-6, atol=1e-6, events = list(data = feed.size)))
  if(attributes(Size.8)$istate[1] != 2)(return(Size.0))
  
  capture.output(Size.8U <- lsoda(initsU, 28:duration, func = "derivs", dllname = "SizeCompModel_Shrink", 
                   initfunc = "initmod",  nout=4, outnames=c("Survival", "LG", "Survival_C", "LGC"), 
                   maxsteps=5e5, params[1:25],  rtol=1e-6, atol=1e-6, events = list(data = feed.sizeU)))
  if(attributes(Size.8U)$istate[1] != 2)(return(Size.0))
  
  inits["LC"] = 13.4; inits["DC"] = params["DR"]
  initsU["LC"] = 13.4; initsU["DC"] = params["DR"]
  capture.output(Size.13 <- lsoda(inits, 28:duration, func = "derivs", dllname = "SizeCompModel_Shrink", 
                   initfunc = "initmod",  nout=4, outnames=c("Survival", "LG", "Survival_C", "LGC"), 
                   maxsteps=5e5, params[1:25],  rtol=1e-6, atol=1e-6, events = list(data = feed.size)))
  if(attributes(Size.13)$istate[1] != 2)(return(Size.0))
  
  capture.output(Size.13U <- lsoda(initsU, 28:duration, func = "derivs", dllname = "SizeCompModel_Shrink", 
                    initfunc = "initmod",  nout=4, outnames=c("Survival", "LG", "Survival_C", "LGC"), 
                    maxsteps=5e5, params[1:25],  rtol=1e-6, atol=1e-6, events = list(data = feed.sizeU)))
  if(attributes(Size.13U)$istate[1] != 2)(return(Size.0))
  
  inits["LC"] = 16.9; inits["DC"] = params["DR"]
  initsU["LC"] = 16.9; initsU["DC"] = params["DR"]
  capture.output(Size.16 <- lsoda(inits, 28:duration, func = "derivs", dllname = "SizeCompModel_Shrink", 
                   initfunc = "initmod",  nout=4, outnames=c("Survival", "LG", "Survival_C", "LGC"), 
                   maxsteps=5e5, params[1:25],  rtol=1e-6, atol=1e-6, events = list(data = feed.size)))
  if(attributes(Size.16)$istate[1] != 2)(return(Size.0))
  
  capture.output(Size.16U <- lsoda(initsU, 28:duration, func = "derivs", dllname = "SizeCompModel_Shrink", 
                    initfunc = "initmod",  nout=4, outnames=c("Survival", "LG", "Survival_C", "LGC"), 
                    maxsteps=5e5, params[1:25],  rtol=1e-6, atol=1e-6, events = list(data = feed.sizeU)))
  if(attributes(Size.16U)$istate[1] != 2)(return(Size.0))
  
  
  result <- rbind(#Infecteds (n=84) Uninfecteds (n=30)
    Size.0, Size.0, Size.0, Size.0, Size.0, Size.0, Size.0, Size.0, Size.0, Size.0, Size.0, Size.0, Size.0, Size.0, Size.0, Size.0,
    Size.2, Size.2, Size.2, Size.2, Size.2, Size.2, Size.2, Size.2, Size.2, Size.2, Size.2, Size.2, Size.2, Size.2, Size.2, Size.2, Size.2, Size.2, Size.2, 
    Size.4, Size.4, Size.4, Size.4, Size.4, Size.4, Size.4, Size.4, Size.4, Size.4, Size.4, Size.4,
    Size.8, Size.8, Size.8, Size.8, Size.8, Size.8, Size.8, Size.8, Size.8, Size.8, Size.8, Size.8, Size.8, 
    Size.13, Size.13, Size.13, Size.13, Size.13, Size.13, Size.13, Size.13, Size.13, Size.13, Size.13, Size.13, Size.13,
    Size.16, Size.16, Size.16, Size.16, Size.16, Size.16, Size.16, Size.16, Size.16, Size.16, Size.16,
    
    Size.0U, Size.0U, Size.0U, Size.0U, Size.0U,
    Size.2U, Size.2U, Size.2U, Size.2U, Size.2U, 
    Size.4U, Size.4U, Size.4U, Size.4U, Size.4U, 
    Size.8U, Size.8U, Size.8U, Size.8U, Size.8U, 
    Size.13U, Size.13U, Size.13U, Size.13U, Size.13U, 
    Size.16U, Size.16U, Size.16U, Size.16U, Size.16U)
  
  result[,"Rtotal"] = result[,"RH"] + result[,"RHC"]
  
  result
  
}


extract.data<-function(data, w.t=7, start.age=0){
  ww<- which(data$time%%w.t==0 & data$time>=start.age)
  data[ww,]
}

## This function takes the simulator, params, inits, etc, and runs a
## forward simulation, and then extracts a subset of the points from
## the simulator (determined by w.t), without noise
make.states<-function(params, inits.R, inits.P, inits.S, duration.R, duration.P, duration.S, feeding.events.R, feeding.events.S, w.t=7){
  result.R = solve.DEB.food(params, inits.R, duration.R, feeding.events.R)
  result.R = extract.data(result.R, start.age=28)
  Survival.R = result.R$Survival
  
  # This fix for survival data is currently specific to the sample size (96) and duration (32) of the experiment
  Survival.R[-((1:96)*32)] =   Survival.R[-((1:96)*32)] - Survival.R[-(1+(0:95)*32)]
  
  result.P = solve.DEB.starve(params, inits.P, duration.P)
  result.P = extract.data(result.P, start.age=14)
  Survival.P = result.P$Survival
  # This fix for survival data is currently specific to the sample size (40) and duration (19) of the experiment
  Survival.P[-((1:40)*19)] =   Survival.P[-((1:40)*19)] - Survival.P[-(1+(0:39)*19)]
  
  result.S = solve.DEB.Size(params, inits.S, duration.S, feeding.events.S)
  result.S = extract.data(result.S, start.age=28)
  Survival.S = result.S$Survival
  Survival.S[(0:83)*16 + 5] = 1 # This correctly overwrites the survival probability to condition on survival to day 56 for infecteds
  # This fix for survival data is currently specific to the sample size (114) and duration (15) of the experiment
  Survival.S[-((1:114)*16)] =   Survival.S[-((1:114)*16)] - Survival.S[-(1+(0:113)*16)]
  
  return(list(time=result.R$time, L=result.R$LG, RH=result.R$RH, RP=result.R$RP, SurvR=Survival.R, 
              L2=result.P$LG, W2=result.P$RP, E2=result.P$RH, SurvP=Survival.P, 
              L3F = result.S$LG, L3C = result.S$LGC, W3 = result.S$RP, E3=result.S$Rtotal, SurvS=Survival.S))
}

# 5 - Prior likelihood
prior.likelihood = function(x){
  prior.lik = with(as.list(x),
                   sum(dbeta(c(yPE, yEF, yEF2, yRP, yVE, mP, eh, k, fe), 1, 1, log=T)) + 
                     sum(dunif(c(sd.L, sd.E, sd.W), min=0, max=10, log=T)) +
                     sum(dunif(c(ph, alpha, iPM, EM, DR, Fh, muD, kR, delta0, hdelta, theta, mR, hb), min=0, max=1000000, log=T)) +
                     dnorm(iM, mean=0.0183, sd=0.0016, log=T) + dnorm(M, mean=0.004, sd=0.00047, log=T) + dnorm(LM, mean=35, sd=2, log=T)
  )
  return(prior.lik)
}


# 6 - data likelihood
#### autocorrelation functions for fitting ##############
dmvnorm2 = function(observations, predictions, sd, rho){
  obs = as.vector(na.omit(observations))
  l = length(obs)
  preds = predictions[1:l]
  vc = diag(l)
  vc <- sd^2/(1 - rho^2) * rho^abs(row(vc)-col(vc))
  
  dmvnorm(x = obs, mean=preds, sigma=vc, log=T)
  
}

dmvnorm_rep = function(observations, predictions, sd, rho, ts.length){
  LL_total = 0 # accumulate Log-like values
  n = length(observations)/ts.length # sample size
  for (i in 1:n){
    LL_individual = dmvnorm2(observations[1:ts.length + ((i-1)*ts.length)], predictions[1:ts.length + ((i-1)*ts.length)], sd, rho)
    LL_total = LL_total + LL_individual
  }
  LL_total
}


full.likelihood<-function(x){
  # simulate data
  sim.data = make.states(x, in.R, in.P, in.S, dur.R, dur.P, dur.S, Feed.R, Feed.S, w.t=7)
  
  # data likelihood
  e.c <- 1
  
  ## observation model
  gammaH <- 0.015 # C content of eggs
  gammaP <- 4e-5 # C content of cercs
  
  ## convert predictions into correct count units
  l.temp<-sim.data$L
  n.temp<-sim.data$RH/gammaH
  w.temp<-sim.data$RP/gammaP
  
  
  l2.temp<-sim.data$L2
  n2.temp<-sim.data$E2/gammaH
  w2.temp<-sim.data$W2/gammaP
  
  l3f.temp<-sim.data$L3F
  l3c.temp<-sim.data$L3C
  n3.temp<-sim.data$E3/gammaH
  w3.temp<-sim.data$W3/gammaP
  
  SR.temp<-sim.data$SurvR
  SP.temp<-sim.data$SurvP
  SS.temp<-sim.data$SurvS
  
  sd.L<-as.numeric(x["sd.L"])
  sd.L2<-as.numeric(x["sd.L2"])
  sd.L3<-as.numeric(x["sd.L3"])
  
  sd.E<-as.numeric(x["sd.E"])
  sd.E2<-as.numeric(x["sd.E2"])
  sd.E3<-as.numeric(x["sd.E3"])
  
  sd.W<-as.numeric(x["sd.W"])
  sd.W2<-as.numeric(x["sd.W2"])
  sd.W3<-as.numeric(x["sd.W3"])

  # Avoids simulations that fell short
  NObs = length(data$L)
  NObs2 = length(data$L2)
  if(length(n.temp) != NObs){print("Simulation 1 too short"); return(-1e6)}
  if(length(n2.temp) != NObs2){print("Simulation 2 too short"); return(-1e6)}    
  
  if(anyNA(l.temp)){print("NaNs in l.temp"); return(-1e6)}
  if(anyNA(l2.temp)){print("NaNs in l2.temp"); return(-1e6)}
  if(anyNA(l3f.temp)){print("NaNs in l3f.temp"); return(-1e6)}
  if(anyNA(l3c.temp)){print("NaNs in l3c.temp"); return(-1e6)}
  if(anyNA(SR.temp)){print("NaNs in SR.temp");return(-1e6)}
  if(anyNA(SP.temp)){print("NaNs in SP.temp");return(-1e6)}
  if(anyNA(SS.temp)){print("NaNs in SS.temp"); return(-1e6)}

  ## likelihoods from food gradient
  llik.L<- sum(dnorm(log(data$L), mean=log(l.temp), sd=sd.L, log=TRUE), na.rm=T)
  llik.Negg<- sum(dnorm(log(data$Negg+e.c), mean=log(n.temp+e.c), sd=sd.E, log=TRUE), na.rm=T)
  llik.Nworms<- sum(dnorm(log(data$Nworms+e.c), mean=log(w.temp+e.c), sd=sd.W, log=TRUE), na.rm=T)

  SR =SR.temp[which(data$Alive == 1)]
  llik.Survival <- sum(log(SR))
  
  ## likelihoods from periodic starvation
  llik.L2<- sum(dnorm(log(data$L2), mean=log(l2.temp), sd=sd.L2, log=TRUE), na.rm=T)
  llik.Negg2<- sum(dnorm(log(data$Negg2+e.c), mean=log(n2.temp+e.c), sd=sd.E2, log=TRUE), na.rm=T)
  llik.Nworms2<- sum(dnorm(log(data$Nworms2+e.c), mean=log(w2.temp+e.c), sd=sd.W2, log=TRUE), na.rm=T)

  SP =SP.temp[which(data$Alive2 == 1)]
  llik.Survival2 <- sum(log(SP)) 
  
  ## likelihoods from size comp
  llik.L3f<- sum(dnorm(log(data$L3), mean=log(l3f.temp), sd=sd.L3, log=TRUE), na.rm=T)
  llik.L3c<- sum(dnorm(log(data$L3C[c(273:1344, 1425:1824)]), mean=log(l3c.temp[c(273:1344, 1425:1824)]), sd=sd.L3, log=TRUE), na.rm=T)
  llik.Nworms3<- sum(dnorm(log(data$Nworms3[1:1344]+e.c), mean=log(w3.temp[1:1344] + e.c), sd=sd.W3, log=TRUE), na.rm=T)
 
  SS =SS.temp[which(data$Alive3 == 1)]
  llik.Survival3 <- sum(log(SS)) 
  
  llik<-(llik.L + llik.Negg + llik.Nworms + llik.Survival + 
           llik.L2 + llik.Negg2 + llik.Nworms2 + llik.Survival2 +
           llik.L3f + llik.L3c + llik.Nworms3 + llik.Survival3)

  if(is.na(llik)|!is.finite(llik)){
    print("Infinite NLL")
    return(-1e6)}
  
  lprior = prior.likelihood(x)
  return(llik + lprior)
}


### Tuning ###

setwd("C:/RData")
samps = readRDS("Sizecomp_kickstart.Rda")
pars = samps$samples[which.max(samps$log.p),]
pars = as.numeric(pars[1:34])
names(pars) = c("iM", "k", "M", "EM", "Fh", "muD", "DR", "fe", "yRP",
                "ph", "yPE", "iPM", "eh", "mP", "alpha", "yEF", "LM",
                "kR", "delta0", "hdelta", "hb", "theta", "mR", "yVE", "yEF2",
                "sd.L", "sd.E", "sd.W", "sd.L2", "sd.E2", "sd.W2", "sd.L3", "sd.E3", "sd.W3")

# Can use the last tuned proposal matrix as the starting point
# variances = samps$cov.jump

# Provides a starting point for uncorrelated parameter proposals
variances = diag(34)
diag(variances) = pars*1e-8

full.likelihood(pars)


### running the mcmc ###
model_fit = MCMC(full.likelihood, init=pars, scale=as.matrix(variances), adapt=50000, acc.rate = 0.3, n=250000)
saveRDS(model_fit, file="Full_fitting_SizeComp.Rda")
