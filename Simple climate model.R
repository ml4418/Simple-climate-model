rm(list=ls())

if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
if(!require(egg)){ install.packages("egg");library(egg)}
if(!require(zoo)){ install.packages("zoo");library(zoo)}
if(!require(deSolve)){ install.packages("deSolve");library(deSolve)}
if(!require(pracma)){ install.packages("pracma");library(pracma)}
if(!require(cowplot)){ install.packages("cowplot");library(cowplot)}
if(!require(magick)){ install.packages("magick");library(magick)}

wd<-"D:/PhD Project/Simple climate model/Data and codes/"
source(paste(wd,"Simple climate model_Functions.R",sep=""), encoding = 'UTF-8')

################ some basic settings

c_heat <- 2.08*10^8/(365*24*60*60)
initial_state <- c(T = 10, v = 1/2000) # Initial state
times <- seq(0, 5000, by = 0.2) # Time points expressed in years


#######################################################################################################################
###############################
###############################     Principle of simple climate model
###############################
#######################################################################################################################

## Greenland temperature and LOVECLIM global mean temperature
T_Greenland<-read.csv(paste(wd,"Input data/Greenland temperature.csv",sep=""),row.names = 1)
LOVECLIM_T_mean_globe<-read.csv(paste(wd,"Input data/LOVECLIM global mean temperature.csv",sep=""),row.names = 1)
LOVECLIM_T_mean_globe$T_anomaly<-LOVECLIM_T_mean_globe$T-LOVECLIM_T_mean_globe[1,"T"]

DO<-read.csv(paste(wd,"Input data/official DO start dates.csv",sep=""),row.names = 1)

DO_lab<-seq(1:20)
DO_data<-cbind.data.frame(DO,DO_lab)
DO<-DO[5:12,];DO_data<-DO_data[5:12,]

xbreak <- 10000*seq(3,5,by=0.2)

p1<-ggplot()+
  geom_point(data = T_Greenland, aes(x = age, y = T),size=0.1)+
  geom_line(data = T_Greenland, aes(x = age, y = T),size=0.2)+
  scale_x_reverse(breaks=xbreak,limits=c(50000,30000),labels=xbreak/1000)+
  labs(y = "Greenland", x = NULL)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size = 16),axis.title = element_text(size = 16))+
  geom_text(data=DO_data,aes(x=DO, label=DO_lab, y=Inf),size=6, hjust=0.4, vjust=-0.6)+  
  coord_cartesian(clip = "off")+
  theme(plot.margin = margin(20, 2, 2, 2, unit = "pt"))

p2<-ggplot()+
  geom_point(data = LOVECLIM_T_mean_globe, aes(x = age, y = T_anomaly),size=0.1)+
  geom_line(data = LOVECLIM_T_mean_globe, aes(x = age, y = T_anomaly),size=0.2)+
  scale_x_reverse(breaks=xbreak,limits=c(50000,30000),labels=xbreak/1000)+
  labs(y = "Global mean", x = "Age (kyr BP)")+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size = 16),axis.title = element_text(size = 16))+
  coord_cartesian(clip = "off")

p<-ggarrange(p1,p2)

ggsave(file=paste(wd,"Input data/fast slow pattern.jpeg",sep=""),p,width=7,height=5,dpi=700)


## c(T)
param_example<-data.frame(C_store=20^2,
                          T_optim=10,
                          T_bal=10,
                          kT=10^2,
                          c_max=4,
                          c_other=0.4,
                          B=1.51)

plot_c_T<-data.frame(T=10+seq(-0.5,0.5,by=0.005))
plot_c_T$c<-get_c(plot_c_T$T,param_example)

p<-ggplot(data=plot_c_T,aes(T,c))+theme_classic()+
  geom_line(size=0.7)+labs(y="c(T)")+ylim(-1,6)+xlim(10-0.5,10+1.3)+
  theme(axis.text = element_blank(),axis.ticks = element_blank())+
  geom_vline(xintercept=as.numeric(param_example$T_optim),linetype="dashed")+
  geom_text(size=4,aes(x=as.numeric(param_example$T_optim),y=-1),label=expression(T[optim]))+
  geom_text(size=4,aes(x=quantile(T, 0.1),y=as.numeric(param_example$c_other)),label=expression(c[other]))+
  geom_text(size=4,aes(x=as.numeric(param_example$T_optim),y=as.numeric(param_example$c_other+param_example$c_max)),label=expression(c[other]+c[max]))+
  geom_hline(yintercept=0,col="gray",size=0.3)+
  geom_text(size=4,aes(x=10+0.8,y=0+0.25),label=expression("case:"~c[other]>=B))+
  geom_hline(yintercept=2.5,col="gray",size=0.3)+
  geom_text(size=4,aes(x=10+0.8,y=2.5+0.25),label=expression("case:"~c[other]~"< B <"~c[other]+c[max]))+
  geom_hline(yintercept=5,col="gray",size=0.3)+
  geom_text(size=4,aes(x=10+0.8,y=5+0.25),label=expression("case:"~c[other]+c[max]<=B))


ggsave(file=paste(wd,"Input data/c(T) shape.jpeg",sep=""),p,width=5,height=2.5,dpi=700)


#######################################################################################################################
###############################
###############################     Three cases of the relationship between c(T) and B
###############################
#######################################################################################################################

param_case<-data.frame(C_store=20^2,
                       T_optim=10,
                       T_bal=10,
                       kT=10^2,
                       c_max=rep(NA,3),
                       c_other=rep(NA,3),
                       B=1.51)

param_case$T1<- param_case$T_optim-sqrt(log((param_case$B-param_case$c_other)/param_case$c_max)/(-param_case$kT))
param_case$T2<- param_case$T_optim+sqrt(log((param_case$B-param_case$c_other)/param_case$c_max)/(-param_case$kT))
(param_case$T2-param_case$T1)/2
param_case$T_eq<-param_case$T_bal
param_case$T_centre<-param_case$T_bal

param_case[1,c("c_max","c_other")]<-c(1,0.5) #c_other+c_max ≤ B
param_case[2,c("c_max","c_other")]<-c(1,1.52) #c_other ≥ B
param_case[3,c("c_max","c_other")]<-c(4,0.4) #c_other < B < c_other+c_max

# Solve the differential equations
solution_case1 <- ode(y = initial_state, times = times, func = diffeq, parms = param_case[1,])
solution_case1<-as.data.frame(solution_case1)

solution_case2 <- ode(y = initial_state, times = times, func = diffeq, parms = param_case[2,])
solution_case2<-as.data.frame(solution_case2)

solution_case3 <- ode(y = initial_state, times = times, func = diffeq, parms = param_case[3,])
solution_case3<-as.data.frame(solution_case3)

# Plot into figures
p_left_1<-plot_solution_figure(solution_case1,param_case[1,])[[1]]+ggtitle(expression(italic("(a)")~c[other]+c[max]<=B)) +xlim(0,5000)+theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank()) 
p_left_2<-plot_solution_figure(solution_case2,param_case[2,])[[1]]+ggtitle(expression(italic("(b)")~c[other]>=B))+xlim(0,5000)+theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank())
p_left_3<-plot_solution_figure(solution_case3,param_case[3,])[[1]]+ggtitle(expression(italic("(c)")~c[other]~"< B <"~c[other]+c[max]))+xlim(0,5000)+theme(axis.text=element_blank(),axis.ticks=element_blank())

p_right_1<-plot_solution_figure(solution_case1,param_case[1,])[[2]]+theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank())
p_right_2<-plot_solution_figure(solution_case2,param_case[2,])[[2]]+theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank())
p_right_3<-plot_solution_figure(solution_case3,param_case[3,])[[2]]+theme(axis.text=element_blank(),axis.ticks=element_blank())

p<-ggarrange(p_left_1,p_right_1,
             p_left_2,p_right_2,
             p_left_3,p_right_3,
             ncol=2)

ggsave(file=paste(wd,"Output data/Plots/Figure 2_Three cases of the relationship between c(T) and B.jpeg",sep=""),p,width=8,height=9,dpi=700)




#######################################################################################################################
######################
######################     Different situations of T_eq in the case c_other < B < c_other+c_max
######################
#######################################################################################################################

param_situation<-data.frame(C_store=20^2,
                            T_optim=10,
                            T_bal=10+seq(-0.2,0.2,by=0.1),
                            kT=10^2,
                            c_max=4,
                            c_other=0.4,
                            B=1.51)

param_situation$T1<- param_situation$T_optim-sqrt(log((param_situation$B-param_situation$c_other)/param_situation$c_max)/(-param_situation$kT))
param_situation$T2<- param_situation$T_optim+sqrt(log((param_situation$B-param_situation$c_other)/param_situation$c_max)/(-param_situation$kT))
(param_situation$T2-param_situation$T1)/2
param_situation$T_eq<-param_situation$T_bal
param_situation$T_centre<-NA


# Solve the differential equations
solution_situation1 <- ode(y = initial_state, times = times, func = diffeq, parms = param_situation[1,])
solution_situation1<-as.data.frame(solution_situation1)

solution_situation2 <- ode(y = initial_state, times = times, func = diffeq, parms = param_situation[2,])
solution_situation2<-as.data.frame(solution_situation2)

solution_situation3 <- ode(y = initial_state, times = times, func = diffeq, parms = param_situation[3,])
solution_situation3<-as.data.frame(solution_situation3)

solution_situation4 <- ode(y = initial_state, times = times, func = diffeq, parms = param_situation[4,])
solution_situation4<-as.data.frame(solution_situation4)

solution_situation5 <- ode(y = initial_state, times = times, func = diffeq, parms = param_situation[5,])
solution_situation5<-as.data.frame(solution_situation5)


# Get T_centre
stable_solution2<-get_stable_cycle(solution_situation2)
param_situation[2,"T_centre"]<-analysis_of_stable_cycle(stable_solution2,T1=param_situation[2,"T1"],T2=param_situation[2,"T2"])[["property"]][,"T_centre"]

stable_solution3<-get_stable_cycle(solution_situation3)
param_situation[3,"T_centre"]<-analysis_of_stable_cycle(stable_solution3,T1=param_situation[3,"T1"],T2=param_situation[3,"T2"])[["property"]][,"T_centre"]

stable_solution4<-get_stable_cycle(solution_situation4)
param_situation[4,"T_centre"]<-analysis_of_stable_cycle(stable_solution4,T1=param_situation[4,"T1"],T2=param_situation[4,"T2"])[["property"]][,"T_centre"]


# Plot into figures
T_lim<-c(10-0.5,10+0.5);v_lim<-c(-0.07,0.07)

p_left_1<-plot_solution_figure(solution_situation1,param_situation[1,])[[1]]+ggtitle(expression(italic("(a)")~T[eq]<T[1]))  +xlim(0,5000)+scale_y_continuous(breaks=10+seq(-0.4,0.4,by=0.2),limits=c(T_lim[1],T_lim[2]))+theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank())
p_left_2<-plot_solution_figure(solution_situation2,param_situation[2,])[[1]]+ggtitle(expression(italic("(b)")~T[1]~"≤"~T[eq]~"<"~T[optim])) +xlim(0,5000)+scale_y_continuous(breaks=10+seq(-0.4,0.4,by=0.2),limits=c(T_lim[1],T_lim[2]))+theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank())
p_left_3<-plot_solution_figure(solution_situation3,param_situation[3,])[[1]]+ggtitle(expression(italic("(c)")~T[eq]~"="~T[optim])) +xlim(0,5000)+scale_y_continuous(breaks=10+seq(-0.4,0.4,by=0.2),limits=c(T_lim[1],T_lim[2]))+theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank())
p_left_4<-plot_solution_figure(solution_situation4,param_situation[4,])[[1]]+ggtitle(expression(italic("(d)")~T[optim]~"<"~T[eq]~"≤"~T[2]))+xlim(0,5000)+scale_y_continuous(breaks=10+seq(-0.4,0.4,by=0.2),limits=c(T_lim[1],T_lim[2]))+theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank())
p_left_5<-plot_solution_figure(solution_situation5,param_situation[5,])[[1]]+ggtitle(expression(italic("(e)")~T[eq]>T[2])) +xlim(0,5000)+scale_y_continuous(breaks=10+seq(-0.4,0.4,by=0.2),limits=c(T_lim[1],T_lim[2]))+theme(axis.text=element_blank(),axis.ticks=element_blank())

p_right_1<-plot_solution_figure(solution_situation1,param_situation[1,])[[2]]+scale_x_continuous(breaks=10+seq(-0.4,0.4,by=0.2),limits=c(T_lim[1],T_lim[2]))+ylim(v_lim[1],v_lim[2])+theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank())
p_right_2<-plot_solution_figure(solution_situation2,param_situation[2,])[[2]]+scale_x_continuous(breaks=10+seq(-0.4,0.4,by=0.2),limits=c(T_lim[1],T_lim[2]))+ylim(v_lim[1],v_lim[2])+theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank())
p_right_3<-plot_solution_figure(solution_situation3,param_situation[3,])[[2]]+scale_x_continuous(breaks=10+seq(-0.4,0.4,by=0.2),limits=c(T_lim[1],T_lim[2]))+ylim(v_lim[1],v_lim[2])+theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank())
p_right_4<-plot_solution_figure(solution_situation4,param_situation[4,])[[2]]+scale_x_continuous(breaks=10+seq(-0.4,0.4,by=0.2),limits=c(T_lim[1],T_lim[2]))+ylim(v_lim[1],v_lim[2])+theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank())
p_right_5<-plot_solution_figure(solution_situation5,param_situation[5,])[[2]]+scale_x_continuous(breaks=10+seq(-0.4,0.4,by=0.2),limits=c(T_lim[1],T_lim[2]))+ylim(v_lim[1],v_lim[2])+theme(axis.text=element_blank(),axis.ticks=element_blank())

p<-ggarrange(p_left_1,p_right_1,
             p_left_2,p_right_2,
             p_left_3,p_right_3,
             p_left_4,p_right_4,
             p_left_5,p_right_5,
             ncol=2)

ggsave(file=paste(wd,"Output data/Plots/Figure 3_Different situations of T_eq.jpeg",sep=""),p,width=8,height=10.5,dpi=700)



#######################################################################################################################
######################
######################     GIF for movement of T_eq in the case c_other < B < c_other+c_max
######################
#######################################################################################################################

param_situation<-data.frame(C_store=20^2,
                            T_optim=10,
                            T_bal=10+seq(-0.2,0.2,by=0.01),
                            kT=10^2,
                            c_max=4,
                            c_other=0.4,
                            B=1.51)

param_situation$T1<- param_situation$T_optim-sqrt(log((param_situation$B-param_situation$c_other)/param_situation$c_max)/(-param_situation$kT))
param_situation$T2<- param_situation$T_optim+sqrt(log((param_situation$B-param_situation$c_other)/param_situation$c_max)/(-param_situation$kT))
(param_situation$T2-param_situation$T1)/2
param_situation$T_eq<-param_situation$T_bal
param_situation$T_centre<-NA


# Solve the differential equations
T_lim<-c(10-0.5,10+0.5);v_lim<-c(-0.07,0.07)

for(i in 1:nrow(param_situation)){
  
  solution_situation_each <- ode(y = initial_state, times = times, func = diffeq, parms = param_situation[i,])
  solution_situation_each<-as.data.frame(solution_situation_each)
  list_figure<-plot_solution_figure(solution_situation_each,param_situation[i,])
  
  p_left<-list_figure[[1]]+xlim(0,5000)+ylim(T_lim[1],T_lim[2])+theme(axis.text=element_blank(),axis.ticks=element_blank())
  p_right<-list_figure[[2]]+xlim(T_lim[1],T_lim[2])+ylim(v_lim[1],v_lim[2])+theme(axis.text=element_blank(),axis.ticks=element_blank())
  p_each<-ggarrange(p_left,p_right,ncol=2)
  
  ggsave(file=paste(wd,"Output data/Plots/GIF source file for T_eq/GIF",i,".jpeg",sep=""),p_each,width=7,height=3.3,dpi=700)
  

}


# Set the path to your images
img_files<-paste(wd, "Output data/Plots/GIF source file for T_eq/GIF",c(1:nrow(param_situation)), ".jpeg", sep="")
img_list <- lapply(img_files, image_read)

# Duplicate the last frame to create a pause (e.g., 10 extra frames)
pause_frames <- rep(img_list[[length(img_list)]], times = 10)  # Adjust times for longer pause
img_list_extended <- c(img_list, pause_frames)

# Create animation
animation <- image_animate(image_join(img_list_extended), fps = 1)  # Adjust fps (frames per second)
image_write(animation, paste(wd, "Output data/Plots/Supplementary GIF_Movement of T_eq.gif", sep=""))



#######################################################################################################################
###############################
###############################     Sensitivity test for T1<=T_eq<=T2 in the case c_other < B < c_other+c_max
###############################
#######################################################################################################################

# Define ranges for sensitivity analysis
param_ranges <- list(
  C_store = seq(10,50,by=10)^2,
  T_optim = 10+ seq(-1,1,by=0.1),
  T_bal = 10+ seq(-1,1,by=0.1),
  kT = seq(2,10,by=2)^2,
  c_max = seq(2, 4, by=0.5),
  c_other = seq(0, 1, by=0.2),
  B=1.51
)

param_combinations <- expand.grid(param_ranges)
nrow(param_combinations)

# Rule out combinations, only use combinations where cother < B < cother+cmax and T1<=Tbal<=T2 
param_combinations<-param_combinations[which(param_combinations$c_other<param_combinations$B & param_combinations$c_other+param_combinations$c_max>param_combinations$B),]
param_combinations$T1<- param_combinations$T_optim-sqrt(log((param_combinations$B-param_combinations$c_other)/param_combinations$c_max)/(-param_combinations$kT))
param_combinations$T2<- param_combinations$T_optim+sqrt(log((param_combinations$B-param_combinations$c_other)/param_combinations$c_max)/(-param_combinations$kT))
param_combinations<-param_combinations[which(param_combinations$T_bal>=param_combinations$T1 & param_combinations$T_bal<=param_combinations$T2),]
nrow(param_combinations)

# Check the distribution of L_allow and amplitude threshold
param_combinations$L_allow<-param_combinations$T2-param_combinations$T1
summary(param_combinations$L_allow);hist(param_combinations$L_allow)
summary(param_combinations$T_bal-param_combinations$T_optim);hist(param_combinations$T_bal-param_combinations$T_optim)
summary(abs(param_combinations$T_bal-param_combinations$T_optim));hist(abs(param_combinations$T_bal-param_combinations$T_optim))

# Do sensitivity test, use a longer time series to capture properties
sensitivity_data<-get_sensitivity(param_combinations ,initial_state, times =seq(0, 10000, by = 1))
write.csv(sensitivity_data,paste(wd,"Output data/sensitivity_data.csv",sep=""))

################################################################################
#######################      Get features
################################################################################

sensitivity_data <- read.csv(paste(wd,"Output data/sensitivity_data.csv",sep=""),row.names = 1)
sensitivity_data<-na.omit(sensitivity_data)


#get the period of corresponding simple harmonic oscillation 
sensitivity_data$period_sh<-2*pi*sqrt(c_heat*sensitivity_data$C_store)

#get the estimated amplitude threshold and estimated period
sensitivity_data$est_amplitude_threshold<-sensitivity_data$T1-sensitivity_data$T2+2*(sensitivity_data$c_max/(sensitivity_data$B-sensitivity_data$c_other))*sqrt(pi/sensitivity_data$kT)

for(i in 1:nrow(sensitivity_data)){
  
  if(i%%1000==0){print(i)}
  
  sensitivity_data[i,c("T_min_est")]<-estimate_Tmin(parameters=sensitivity_data[i,])
  sensitivity_data[i,c("T_max_est")]<-estimate_Tmax(parameters=sensitivity_data[i,])
  
  sensitivity_data[i,c("t1_est","t3_est","slow_phase_est")]<-estimate_period(parameters=sensitivity_data[i,],
                                                        T_min=sensitivity_data[i,"T_min_est"],
                                                        T_max=sensitivity_data[i,"T_max_est"],
                                                        bin=0.005)
}

#make sure that the root is not T2 in estimate_Tmin(), and not T1 in estimate_Tmax()
any(sensitivity_data$T_min_est==sensitivity_data$T2)
any(sensitivity_data$T_max_est==sensitivity_data$T1)

write.csv(sensitivity_data,paste(wd,"Output data/sensitivity_data_with_feature.csv",sep=""))

################################################################################
#######################      Plot relationship
################################################################################

sensitivity_data <- read.csv(paste(wd,"Output data/sensitivity_data_with_feature.csv",sep=""),row.names = 1)
summary(sensitivity_data)

############################# General features

p1<-ggplot(data=sensitivity_data,aes(T_bal-T_optim,T_centre-T_optim,col=sqrt(C_store)))+
  theme_bw()+geom_point(size=0.5)+geom_abline(slope=1,intercept=0)+
  scale_color_gradientn(colours  = rev(c("purple","blue","dodgerblue2","green","khaki","gold","orange","red")),na.value = "grey90")+
  labs(x=expression(T[eq]-T[optim]~"(K)"),y=expression(T[centre]-T[optim]~"(K)"))+
  theme(legend.position="none")+
  scale_x_continuous(breaks=seq(-0.5,0.5,by=0.25))+
  scale_y_continuous(breaks=seq(-0.5,0.5,by=0.25))


p2<-ggplot(data=sensitivity_data,aes(est_amplitude_threshold,amplitude,col=sqrt(C_store)))+
  theme_bw()+geom_point(size=0.5)+ 
  scale_color_gradientn(colours  = rev(c("purple","blue","dodgerblue2","green","khaki","gold","orange","red")),na.value = "grey90")+
  labs(x=expression("Estimated threshold of peak to trough amplitude (K)"),y=expression("Peak to trough amplitude (K)"))+
  geom_abline(slope=1,intercept=0)+
  theme(legend.position="none")+
  scale_x_continuous(breaks=c(0,3,6,9,12))+
  scale_y_continuous(breaks=c(0,3,6,9,12))

p3<-ggplot(data=sensitivity_data,aes(slow_phase_est,period,col=sqrt(C_store)))+
  theme_bw()+geom_point(size=0.5)+geom_abline(slope=1,intercept=0)+
  scale_color_gradientn(colours  = rev(c("purple","blue","dodgerblue2","green","khaki","gold","orange","red")),na.value = "grey90")+
  labs(x="Estimated slow phase time (yr)",y="Period (yr)",col=expression(sqrt(C[store])))+
  theme(legend.position = "bottom")

p<-ggarrange(p1,p2,p3,ncol=1,labels=c("(a)","(b)","(c)"),
             label.args = list(gp = grid::gpar(font = 4,cex=1)))

ggsave(file=paste(wd,"Output data/Plots/Figure 4_Centre amplitude period.jpeg",sep=""),p,width=4.5,height=9,dpi=700)


############################# Slow warming and slow cooling

p_t1<-ggplot(data=sensitivity_data,aes(t1_est,t1,col=sqrt(C_store)))+
  theme_bw()+geom_point(size=0.5)+geom_abline(slope=1,intercept=0)+
  scale_color_gradientn(colours  = rev(c("purple","blue","dodgerblue2","green","khaki","gold","orange","red")),na.value = "grey90")+
  labs(x="Estimated slow warming phase time (yr)",y="Slow warming phase time (yr)",col=expression(sqrt(C[store])))+
  theme(legend.position = "none")

p_t3<-ggplot(data=sensitivity_data,aes(t3_est,t3,col=sqrt(C_store)))+
  theme_bw()+geom_point(size=0.5)+geom_abline(slope=1,intercept=0)+
  scale_color_gradientn(colours  = rev(c("purple","blue","dodgerblue2","green","khaki","gold","orange","red")),na.value = "grey90")+
  labs(x="Estimated slow cooling phase time (yr)",y="Slow cooling phase time (yr)",col=expression(sqrt(C[store])))+
  theme(legend.position = "right")


p_slow<-ggarrange(p_t1,p_t3,ncol=2,labels=c("(a)","(b)"),label.args = list(gp = grid::gpar(font = 4,cex=1)))

ggsave(file=paste(wd,"Output data/Plots/Supplementary Figure 1_Slow warming and slow cooling.jpeg",sep=""),p_slow,width=8,height=3.5,dpi=700)

############################# Amplification effect of feedback

p_amplify1<-ggplot(data=sensitivity_data,aes(period/period_sh))+theme_bw()+
  geom_histogram(aes(y = after_stat(count / sum(count))),bins=50,fill="grey",col="black")+
  scale_x_continuous(breaks=c(1,3,5,7))+geom_vline(xintercept = 1)+
  labs(y="Frequency",x=expression("The ratio of period to period"[sh]))

p_amplify2<-ggplot(data=sensitivity_data,aes((t1+t3)/period,  period/period_sh,col=sqrt(C_store)))+
  theme_bw()+geom_point(size=0.5)+
  scale_color_gradientn(colours  = rev(c("purple","blue","dodgerblue2","green","khaki","gold","orange","red")),na.value = "grey90")+
  labs(x="The proportion of slow phase in period",y=expression("The ratio of period to period"[sh]),col=expression(sqrt(C[store])))+
  theme(legend.position = "right")+
  scale_x_continuous(breaks=seq(0.5,1,by=0.1))+
  scale_y_continuous(breaks=seq(1,8,by=1))

ggsave(file=paste(wd,"Output data/Plots/Supplementary Figure 2_Amplification of the period.jpeg",sep=""),p_amplify2,width=7,height=5,dpi=700)

#######################################################################################################################
###########################
###########################     Fast-slow dynamics of the limit cycle
###########################
#######################################################################################################################

################################ fast and slow manifold

param_example<-data.frame(C_store=20^2,
                          T_optim=10,
                          T_bal=10,
                          kT=10^2,
                          c_max=4,
                          c_other=0.4,
                          B=1.51)

param_example$T1<- param_example$T_optim-sqrt(log((param_example$B-param_example$c_other)/param_example$c_max)/(-param_example$kT))
param_example$T2<- param_example$T_optim+sqrt(log((param_example$B-param_example$c_other)/param_example$c_max)/(-param_example$kT))
(param_example$T2-param_example$T1)/2
param_example$T_eq<-param_example$T_bal
param_example$T_centre<-param_example$T_bal

# Calculate
solution <- ode(y = initial_state, times = times, func = diffeq, parms = param_example)
solution<-as.data.frame(solution)

# Define T range for plotting
T_vals <- c(seq(param_example$T1 - 2, param_example$T1 - 0.01, length.out = 10^3),
            seq(param_example$T1 - 0.01, param_example$T1 + 0.01, length.out = 10^3),
            seq(param_example$T1 + 0.01, param_example$T2 - 0.01, length.out = 10^3),
            seq(param_example$T2 - 0.01, param_example$T2 + 0.01, length.out = 10^3),
            seq(param_example$T2 + 0.01, param_example$T2 + 2, length.out = 10^3))

manifold<-data.frame(T_vals=T_vals)
manifold[,"slow_manifold"]<-get_slow_manifold(T_vals,param_example)
manifold[,c("fast_manifold_upper","fast_manifold_lower")]<-get_fast_manifold(T_vals,param_example)

# Plot
T_lim<-c(10-0.7,10+0.7);v_lim<-c(-0.08,0.08)

p_manifold<-ggplot()+theme_bw()+
  geom_point(data=solution,aes(T,v),size=0.05)+
  geom_path(data=solution,aes(T,v),size=0.05)+
  geom_line(data=manifold,aes(T_vals,slow_manifold),col="dodgerblue3",size=0.3)+
  geom_line(data=manifold,aes(T_vals,fast_manifold_upper),col="red",size=0.3)+
  geom_line(data=manifold,aes(T_vals,fast_manifold_lower),col="red",size=0.3)+
  scale_x_continuous(breaks=10+seq(-0.4,0.4,by=0.2),limits=c(T_lim[1],T_lim[2]))+ylim(v_lim[1],v_lim[2])+
  theme(axis.text=element_blank(),axis.ticks=element_blank())+
  geom_rect(data=param_example,aes(xmin=as.numeric(param_example[,"T1"]),xmax=as.numeric(param_example[,"T2"]),
                                   ymin=-Inf,ymax=Inf ), fill="gray",alpha=0.4)


################################ distortion

### solution

param_distortion<-data.frame(C_store=20^2,
                             T_optim=10,
                             T_bal=10,
                             kT=10^2,
                             c_max=4,
                             c_other=0.4,
                             B=1.51)

param_distortion$T1<- param_distortion$T_optim-sqrt(log((param_distortion$B-param_distortion$c_other)/param_distortion$c_max)/(-param_distortion$kT))
param_distortion$T2<- param_distortion$T_optim+sqrt(log((param_distortion$B-param_distortion$c_other)/param_distortion$c_max)/(-param_distortion$kT))
(param_distortion$T2-param_distortion$T1)/2
param_distortion$T_eq<-param_distortion$T_bal

solution <- ode(y = initial_state, times = seq(0, 5000, by = 0.1), func = diffeq, parms = param_distortion)
solution<-as.data.frame(solution)
stable_solution<-get_stable_cycle(solution)
analysis<-analysis_of_stable_cycle(stable_solution,T1=param_distortion$T1,T2=param_distortion$T2)

### solution of corresponding simple harmonic oscillation with same amplitude

initial_state_sh<-c(T = analysis[["property"]]$T_min, v = 0) # Initial state

solution_sh <- ode(y = initial_state_sh, times = seq(0, 5000, by = 0.1), func = diffeq_sh, parms = param_distortion)
solution_sh<-as.data.frame(solution_sh)
stable_solution_sh<-get_stable_cycle(solution_sh)
analysis_sh<-analysis_of_stable_cycle(stable_solution_sh,T1=param_distortion$T1,T2=param_distortion$T2)

###  use only only cycle of the stable part

solution_plot<-solution[solution$time>=analysis[["t_Tmin"]][1]&solution$time<=analysis[["t_Tmin"]][2],]
solution_sh_plot<-solution_sh[solution_sh$time>=analysis_sh[["t_Tmin"]][1]&solution_sh$time<=analysis_sh[["t_Tmin"]][2],]

### plot the principle of distortion

time_diff<-analysis[["t_Tmax"]][1]-analysis_sh[["t_Tmax"]][1]

p_time<-ggplot()+theme_bw()+xlim(analysis[["t_Tmin"]][1]-200,analysis[["t_Tmin"]][2]+200)+
  geom_point(data=solution_plot,aes(time,T),size=0.1)+
  geom_line(data=solution_plot,aes(time,T),size=0.1)+
  geom_point(data=solution_sh_plot,aes(time+time_diff,T),size=0.1,col="gray")+
  geom_line(data=solution_sh_plot,aes(time+time_diff,T),size=0.1,col="gray")+
  geom_rect(aes(ymin=as.numeric(param_distortion[,"T1"]),ymax=as.numeric(param_distortion[,"T2"]),
                xmin=-Inf,xmax=Inf ), fill="gray",alpha=0.4)+
  geom_text(aes(x=analysis[["t_Tmin"]][1],y=solution[solution$time==analysis[["t_Tmin"]][1],"T"]),label=expression(T[min]),col="red",size=5)+
  geom_text(aes(x=analysis[["t_T1_warming"]][1],y=solution[solution$time==analysis[["t_T1_warming"]][1],"T"]),label=expression(T[1]),col="red",size=5)+
  geom_text(aes(x=analysis[["t_T2_warming"]][1],y=solution[solution$time==analysis[["t_T2_warming"]][1],"T"]),label=expression(T[2]),col="red",size=5)+
  geom_text(aes(x=analysis[["t_Tmax"]][1],y=solution[solution$time==analysis[["t_Tmax"]][1],"T"]),label=expression(T[max]),col="red",size=5)+
  geom_text(aes(x=analysis[["t_T2_cooling"]][1],y=solution[solution$time==analysis[["t_T2_cooling"]][1],"T"]),label=expression(T[2]),col="red",size=5)+
  geom_text(aes(x=analysis[["t_T1_cooling"]][1],y=solution[solution$time==analysis[["t_T1_cooling"]][1],"T"]),label=expression(T[1]),col="red",size=5)+
  geom_text(aes(x=analysis[["t_Tmin"]][2],y=solution[solution$time==analysis[["t_Tmin"]][2],"T"]),label=expression(T[min]),col="red",size=5)+
  theme(axis.text=element_blank(),axis.ticks=element_blank())

y_text_position<-min(solution_plot$v)

p_phase<-ggplot()+theme_bw()+xlim(analysis[["property"]][["T_min"]]-0.1,analysis[["property"]][["T_max"]]+0.1)+
  geom_point(data=solution_plot,aes(T,v),size=0.1)+
  geom_path(data=solution_plot,aes(T,v),size=0.1)+
  geom_point(data=solution_sh_plot,aes(T,v),size=0.1,col="gray")+
  geom_path(data=solution_sh_plot,aes(T,v),size=0.1,col="gray")+
  geom_rect(aes(xmin=as.numeric(param_distortion[,"T1"]),xmax=as.numeric(param_distortion[,"T2"]),
                ymin=-Inf,ymax=Inf ), fill="gray",alpha=0.4)+
  geom_text(aes(x=analysis[["property"]][["T_min"]],y=0),label=expression(T[min]),col="red",size=5)+
  geom_text(aes(x=analysis[["property"]][["T_max"]],y=0),label=expression(T[max]),col="red",size=5)+
  geom_text(aes(x=solution[solution$time==analysis[["t_T1_warming"]][1],"T"],y=solution[solution$time==analysis[["t_T1_warming"]][1],"v"]),label=expression(T[1]),col="red",size=5)+
  geom_text(aes(x=solution[solution$time==analysis[["t_T2_warming"]][1],"T"],y=solution[solution$time==analysis[["t_T2_warming"]][1],"v"]),label=expression(T[2]),col="red",size=5)+
  geom_text(aes(x=solution[solution$time==analysis[["t_T1_cooling"]][1],"T"],y=solution[solution$time==analysis[["t_T1_cooling"]][1],"v"]),label=expression(T[1]),col="red",size=5)+
  geom_text(aes(x=solution[solution$time==analysis[["t_T2_cooling"]][1],"T"],y=solution[solution$time==analysis[["t_T2_cooling"]][1],"v"]),label=expression(T[2]),col="red",size=5)+
  theme(axis.text=element_blank(),axis.ticks=element_blank())




########################### plot together

p<-ggarrange(p_manifold,p_time,p_phase,p_amplify1,ncol=2,labels=c("(a)","(b)","(c)","(d)"),
             label.args = list(gp = grid::gpar(font = 4,cex=1)))


ggsave(file=paste(wd,"Output data/Plots/Figure 5_Distortion and amplification.jpeg",sep=""),p,width=8,height=6,dpi=700)




