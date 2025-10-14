###############################################################################################
##### Define the differential equation, x is T (unit:K), y is v (unit:K/yr)
###############################################################################################

diffeq <- function(time, state, parameters) {
  
  T <- state[1]
  v <- state[2]
  
  c_heat <- 2.08*10^8/(365*24*60*60)
  B <-parameters[,"B"]
  C_store<- parameters[,"C_store"]
  T_optim <- parameters[,"T_optim"]
  T_bal <- parameters[,"T_bal"]
  kT <- parameters[,"kT"]
  c_max <- parameters[,"c_max"]
  c_other <- parameters[,"c_other"]
  
  c<-c_other+c_max*exp(-kT*(T-T_optim)^2)
  
  dTdt <- v
  dvdt <- (c-B)*v/c_heat - (T-T_bal)/(c_heat*C_store)
  
  return(list(c(dTdt, dvdt)))
  
}

################################################################################################################################
##### Define the differential equation of corresponding simple harmonic oscillation, x is T (unit:K), y is v (unit:K/yr)
################################################################################################################################

diffeq_sh <- function(time, state, parameters) {
  
  T <- state[1]
  v <- state[2]
  
  c_heat <- 2.08*10^8/(365*24*60*60)
  B <-parameters[,"B"]
  C_store<- parameters[,"C_store"]
  T_bal <- parameters[,"T_bal"]

  dTdt <- v
  dvdt <-  - (T-T_bal)/(c_heat*C_store)
  
  return(list(c(dTdt, dvdt)))
  
}

###############################################################################################
##### Get the stable cycle from the solution
###############################################################################################

get_stable_cycle<-function(solution,diff_amplitude_threshold=10^(-2),diff_period_threshold=25){
  
  # Find out trough and peak
  solution_Tmin_Tmax<-solution[which(diff(sign(solution$v)) != 0) + 1,]
  amplitude_series<-abs(diff(solution_Tmin_Tmax$T))
  
  # Set a threshold for stable cycles
  diff_amplitude<-abs(diff(amplitude_series))
  stable_indices <- which(diff_amplitude <= diff_amplitude_threshold)
  
  if (length(stable_indices) > 0 & length(amplitude_series)>2) {
    first_stable_index <- stable_indices[1] + 1  
    stable_solution_Tmin_Tmax <- solution_Tmin_Tmax[first_stable_index:nrow(solution_Tmin_Tmax), ]
    stable_solution<-solution[solution$time>=stable_solution_Tmin_Tmax[1,"time"],]
    
    #check if period is stable too
    period_series<-diff(stable_solution_Tmin_Tmax$time, lag = 2)
    diff_period<-abs(diff(period_series))
    
    if(!all(diff_period <= diff_period_threshold)){
      message("System does not reach stability within the given time.")
      stable_solution<-NA
    }
    
  } else {
    message("System does not reach stability within the given time.")
    stable_solution <- NA  
  }
  
  return(stable_solution)
  
}

###############################################################################################
##### Get analysis of the stable cycle of the solution
###############################################################################################

analysis_of_stable_cycle<-function(stable_solution,T1,T2){
  
  ########### Get amplitude and period
  
  stable_solution_Tmin_Tmax<-stable_solution[which(diff(sign(stable_solution$v)) != 0) + 1,]
  
  period<-mean(diff(stable_solution_Tmin_Tmax$time, lag = 2))
  amplitude<-mean(abs(diff(stable_solution_Tmin_Tmax$T)))
  
  T_min<-min(stable_solution_Tmin_Tmax$T)
  T_max<-max(stable_solution_Tmin_Tmax$T)
  T_centre<-(T_max+T_min)/2
  
  v_min<-min(stable_solution$v)
  v_max<-max(stable_solution$v)
  
  ########### Analysis of different parts
  
  stable_solution_Tmin<-stable_solution[which(diff(sign(stable_solution$v)) > 0) + 1,]
  stable_solution_Tmax<-stable_solution[which(diff(sign(stable_solution$v)) < 0) + 1,]
  stable_solution_T1_warming<-stable_solution[which(diff(sign(stable_solution$T-T1)) > 0) + 1,]
  stable_solution_T1_cooling<-stable_solution[which(diff(sign(stable_solution$T-T1)) < 0) + 1,]
  stable_solution_T2_warming<-stable_solution[which(diff(sign(stable_solution$T-T2)) > 0) + 1,]
  stable_solution_T2_cooling<-stable_solution[which(diff(sign(stable_solution$T-T2)) < 0) + 1,]
  
  # get the time point
  t_Tmin<-stable_solution_Tmin$time
  t_Tmax<-stable_solution_Tmax$time
  t_T1_warming<-stable_solution_T1_warming$time
  t_T1_cooling<-stable_solution_T1_cooling$time
  t_T2_warming<-stable_solution_T2_warming$time
  t_T2_cooling<-stable_solution_T2_cooling$time
  
  # make sure t_Tmin is the first that occurs
  t_Tmax<-t_Tmax[which(t_Tmax>t_Tmin[1])]
  t_T1_warming<-t_T1_warming[which(t_T1_warming>t_Tmin[1])]
  t_T1_cooling<-t_T1_cooling[which(t_T1_cooling>t_Tmin[1])]
  t_T2_warming<-t_T2_warming[which(t_T2_warming>t_Tmin[1])]
  t_T2_cooling<-t_T2_cooling[which(t_T2_cooling>t_Tmin[1])]
  
  # get the duration for each part
  t1<-t_T1_warming[1]-t_Tmin[1] #period_slow_warming
  t2<-t_Tmax[1]-t_T1_warming[1] #period_fast_warming
  t3<-t_T2_cooling[1]-t_Tmax[1] #period_slow_cooling
  t4<-t_Tmin[2]-t_T2_cooling[1] #period_fast_cooling
  
  property<-data.frame(period,amplitude,T_min,T_max,T_centre,v_min,v_max,
                       t1,t2,t3,t4)
  
  output<-list(property=property,
               t_Tmin=t_Tmin,t_Tmax=t_Tmax,
               t_T1_warming=t_T1_warming,t_T1_cooling=t_T1_cooling,
               t_T2_warming=t_T2_warming,t_T2_cooling=t_T2_cooling)
  
  return(output)
  
}

###############################################################################################
##### Get the c(T)
###############################################################################################

get_c<-function(T,parameters){
  
  T_optim <- parameters[,"T_optim"]
  kT <- parameters[,"kT"]
  c_max <- parameters[,"c_max"]
  c_other <- parameters[,"c_other"]
  
  c<-c_other+c_max*exp(-kT*(T-T_optim)^2)
  
  return(c)
  
}

###############################################################################################
##### Sensitivity test
###############################################################################################

get_sensitivity<-function(param_combinations ,initial_state, times){
  
  output<-data.frame()
  
  for(i in 1:nrow(param_combinations)){
    
    print(i)
    
    # Solve the differential equations
    solution <- ode(y = initial_state, times = times, func = diffeq, parms = param_combinations[i,])
    solution<-as.data.frame(solution)
    
    # Get only the stable part
    stable_solution<-get_stable_cycle(solution)

    if(!all(is.na(stable_solution))){
      
      # Analysis of the stable part
      property<-analysis_of_stable_cycle(stable_solution,T1=param_combinations[i,"T1"],T2=param_combinations[i,"T2"])[["property"]]
      
      # Output
      each_output<-cbind.data.frame(param_combinations[i,],property)
      output<-rbind.data.frame(output,each_output)
    }
    
  }
  
  return(output)

}

###############################################################################################
##### Plot the solution
###############################################################################################

plot_solution_figure<-function(solution_plot,param_solution){
  
  # Plot into figures
  x_text_position<-max(solution_plot$time)
  
  p1<-ggplot()+theme_bw()+
    geom_point(data=solution_plot,aes(time,T),size=0.1)+
    geom_line(data=solution_plot,aes(time,T),size=0.1)+
    geom_rect(aes(ymin=as.numeric(param_solution[,"T1"]),ymax=as.numeric(param_solution[,"T2"]),
                  xmin=-Inf,xmax=Inf ), fill="gray",alpha=0.4)+
    geom_hline(size=0.5,yintercept=as.numeric(param_solution[,"T_optim"]),linetype="dashed")+
    geom_text(aes(y=as.numeric(param_solution[,"T_optim"]),x=0.1*x_text_position),label=expression(T[optim]))+
    geom_hline(size=0.5,yintercept=as.numeric(param_solution[,"T1"]),linetype="dashed")+
    geom_text(aes(y=as.numeric(param_solution[,"T1"]),x=0.1*x_text_position),label=expression(T[1]))+
    geom_hline(size=0.5,yintercept=as.numeric(param_solution[,"T2"]),linetype="dashed")+
    geom_text(aes(y=as.numeric(param_solution[,"T2"]),x=0.1*x_text_position),label=expression(T[2]))+
    geom_hline(size=0.8,yintercept=as.numeric(param_solution[,"T_eq"]),linetype="dashed",col="dodgerblue4")+
    geom_text(aes(y=as.numeric(param_solution[,"T_eq"]),x=x_text_position),label=expression(T[eq]),col="dodgerblue4")+
    geom_hline(size=0.5,yintercept=as.numeric(param_solution[,"T_centre"]),linetype="dashed",col="red")+
    geom_text(aes(y=as.numeric(param_solution[,"T_centre"]),x=0.8*x_text_position),label=expression(T[centre]),col="red")
  
  # plot phase diagram
  y_text_position<-min(solution_plot$v)
  
  p2<-ggplot() + theme_bw()+ 
    geom_point(data=solution_plot, aes(x = T, y = v),size=0.05)+
    geom_path(data=solution_plot, aes(x = T, y = v),size=0.05)+
    geom_rect(aes(xmin=as.numeric(param_solution[,"T1"]),xmax=as.numeric(param_solution[,"T2"]),
                  ymin=-Inf,ymax=Inf ), fill="gray",alpha=0.4)+
    geom_vline(size=0.5,xintercept=as.numeric(param_solution[,"T_optim"]),linetype="dashed")+
    geom_text(aes(x=as.numeric(param_solution[,"T_optim"]),y=0.8*y_text_position),label=expression(T[optim]))+
    geom_vline(size=0.5,xintercept=as.numeric(param_solution[,"T1"]),linetype="dashed")+
    geom_text(aes(x=as.numeric(param_solution[,"T1"]),y=0.8*y_text_position),label=expression(T[1]))+
    geom_vline(size=0.5,xintercept=as.numeric(param_solution[,"T2"]),linetype="dashed")+
    geom_text(aes(x=as.numeric(param_solution[,"T2"]),y=0.8*y_text_position),label=expression(T[2]))+
    geom_vline(size=0.8,xintercept=as.numeric(param_solution[,"T_eq"]),linetype="dashed",col="dodgerblue4")+
    geom_text(aes(x=as.numeric(param_solution[,"T_eq"]),y=-0.8*y_text_position),label=expression(T[eq]),col="dodgerblue4")+
    geom_vline(size=0.5,xintercept=as.numeric(param_solution[,"T_centre"]),linetype="dashed",col="red")+
    geom_text(aes(x=as.numeric(param_solution[,"T_centre"]),y=y_text_position),label=expression(T[centre]),col="red")
  
  return(list(p1,p2))
  
}

###############################################################################################
##### Get the slow and fast manifold
###############################################################################################

# Get slow manifold
get_slow_manifold<-function(T,parameters){
  
  c_heat <- 2.08*10^8/(365*24*60*60)
  B <-parameters[,"B"]
  C_store<- parameters[,"C_store"]
  T_optim <- parameters[,"T_optim"]
  T_bal <- parameters[,"T_bal"]
  kT <- parameters[,"kT"]
  c_max <- parameters[,"c_max"]
  c_other <- parameters[,"c_other"]
  
  c<-c_other+c_max*exp(-kT*(T-T_optim)^2)
  
  slow_manifold<-(T - T_bal) / (C_store * (c - B))
  
  return(slow_manifold)
  
}

# Get fast manifold
get_fast_manifold<-function(T,parameters){
  
  c_heat <- 2.08*10^8/(365*24*60*60)
  B <-parameters[,"B"]
  C_store<- parameters[,"C_store"]
  T_optim <- parameters[,"T_optim"]
  T_bal <- parameters[,"T_bal"]
  kT <- parameters[,"kT"]
  c_max <- parameters[,"c_max"]
  c_other <- parameters[,"c_other"]
  T1 <- parameters[,"T1"]
  T2 <- parameters[,"T2"]
  
  diff_upper<-erf(sqrt(kT)*(T-T_optim))-erf(sqrt(kT)*(T1-T_optim))
  fast_manifold_upper<- ( (c_other-B)/c_heat )*(T-T1) + (c_max/c_heat)*sqrt(pi/kT)*(1/2)*diff_upper
  
  diff_lower<-erf(sqrt(kT)*(T-T_optim))-erf(sqrt(kT)*(T2-T_optim))
  fast_manifold_lower<- ( (c_other-B)/c_heat )*(T-T2) + (c_max/c_heat)*sqrt(pi/kT)*(1/2)*diff_lower
  
  fast_manifold<-cbind.data.frame(fast_manifold_upper,fast_manifold_lower)
  
  return(fast_manifold)
  
}

###############################################################################################
##### Estimate the peak and trough
###############################################################################################

# T_max≈T_1+c_max/(2(B-c_other)) √(π/k_T ) (Φ(T_max )-Φ(T_1 ))
solve_root_Tmax<-function(x,parameters){
  
  B <-parameters[,"B"]
  T_optim <- parameters[,"T_optim"]
  kT <- parameters[,"kT"]
  c_max <- parameters[,"c_max"]
  c_other <- parameters[,"c_other"]
  T1 <- parameters[,"T1"]
  T2 <- parameters[,"T2"]
  
  phi_x <- erf(sqrt(kT)*(x-T_optim))
  phi_T1 <- erf(sqrt(kT)*(T1-T_optim))
  
  output<-x- T1 - ( c_max/(2*(B-c_other)) ) * sqrt(pi/kT) * (phi_x - phi_T1)
  
  return(output)
  
}

# T_min≈T_2+c_max/(2(B-c_other)) √(π/k_T ) (Φ(T_min )-Φ(T_2 ))
solve_root_Tmin<-function(x,parameters){
  
  B <-parameters[,"B"]
  T_optim <- parameters[,"T_optim"]
  kT <- parameters[,"kT"]
  c_max <- parameters[,"c_max"]
  c_other <- parameters[,"c_other"]
  T1 <- parameters[,"T1"]
  T2 <- parameters[,"T2"]
  
  phi_x <- erf(sqrt(kT)*(x-T_optim))
  phi_T2 <- erf(sqrt(kT)*(T2-T_optim))
  
  output<-x- T2 - ( c_max/(2*(B-c_other)) ) * sqrt(pi/kT) * (phi_x - phi_T2)
  
  return(output)
  
}

# Get Tmax 
estimate_Tmax<-function(parameters){
  
  B <-parameters[,"B"]
  T_optim <- parameters[,"T_optim"]
  kT <- parameters[,"kT"]
  c_max <- parameters[,"c_max"]
  c_other <- parameters[,"c_other"]
  T1 <- parameters[,"T1"]
  T2 <- parameters[,"T2"]
  
  lim<-( c_max/(2*(B-c_other)) ) * sqrt(pi/kT) * 2
  output<-uniroot(solve_root_Tmax, interval = c(T1-lim, T1 + lim), parameters = parameters)$root
  
  return(output)
}

# Get Tmin
estimate_Tmin<-function(parameters){
  
  B <-parameters[,"B"]
  T_optim <- parameters[,"T_optim"]
  kT <- parameters[,"kT"]
  c_max <- parameters[,"c_max"]
  c_other <- parameters[,"c_other"]
  T1 <- parameters[,"T1"]
  T2 <- parameters[,"T2"]
  
  lim<-( c_max/(2*(B-c_other)) ) * sqrt(pi/kT) * 2
  output<-uniroot(solve_root_Tmin, interval = c(T2-lim, T2+lim), parameters = parameters)$root
  
  return(output)
}

###############################################################################################
##### Estimate the period
###############################################################################################

# Let the time used at slow warming phase (t1) and slow cooling phase (t3) approximate the period
estimate_period<-function(parameters,T_min,T_max,bin){
  
  T1 <- parameters[,"T1"]
  T2 <- parameters[,"T2"]
  
  I<-data.frame(T_range=seq(T_min,T_max,by=bin))
  I[,"slow_manifold"]<-get_slow_manifold(I$T_range,parameters)

  I1<-I[which(I$T_range>=T_min & I$T_range<=T1),]
  t1<-sum(bin/abs(I1$slow_manifold))

  I3<-I[which(I$T_range>=T2 & I$T_range<=T_max),]
  t3<-sum(bin/abs(I3$slow_manifold))
  
  output<-data.frame(t1,t3,t=t1+t3)
  
  return(output)
  
}


