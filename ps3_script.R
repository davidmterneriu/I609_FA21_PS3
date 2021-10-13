#PS3 Work 

rm(list=ls())
library(readxl)
library(tidyverse)
library(ggplot2)
library(latex2exp)


SIR_fun=function(s0,i0,r0,beta,gamma,dt,Tend){
  
  beta=beta%>%as.numeric()
  gamma=gamma%>%as.numeric()
  
  iterations=Tend/dt
  
  s_pop=s0
  i_pop=i0
  r_pop=r0
  pop_df<-c(0,s_pop, i_pop,r_pop)
  
  
  for(i in 1:iterations){
    
    ds=-beta*s_pop*i_pop
    di=beta*s_pop*i_pop-gamma*i_pop
    dr=gamma*i_pop
    
    s_new=s_pop+ds*dt
    i_new=i_pop+di*dt
    r_new=r_pop+dr*dt

    new_df=c(i*dt,s_new,i_new,r_new)
    pop_df=rbind(pop_df,new_df)
    
    s_pop<-s_new
    i_pop<-i_new
    r_pop<-r_new
  }
  pop_df=as.data.frame(pop_df)
  rownames(pop_df)<-c()
  colnames(pop_df)<-c("Time","Susceptible","Infected","Removed")
  
  
  
  return(pop_df)
  
}

################################################################################
#Q1
################################################################################

par_set1=c(beta=3,gamma=1.6)


q1_df=SIR_fun(s0=0.87,i0=0.13,r0=0,
        beta=par_set1[1],gamma=par_set1[2],
        dt=0.001,Tend = 20)

plot1=q1_df%>%gather(key="Series",value="share",-1)%>%
  mutate(Series=factor(Series,levels = c("Susceptible","Infected","Removed")))%>%
  ggplot(aes(x=Time,y=share,color=Series))+
  geom_line(size=1)+
  ggthemes::scale_color_fivethirtyeight()+
  theme_bw()+
  scale_y_continuous(breaks = seq(0,1,by=0.1))+
  labs(y="Population Share",title = "SIR Model Simulation",
       subtitle = TeX("($\\beta =3$ and $\\gamma =1.6$)"))

################################################################################
#Q2
################################################################################

q1_df%>%filter(Time==5.2)%>%
  select(Infected)%>%
  unlist()%>%
  as.numeric()

################################################################################
#Q3
################################################################################

q1_df%>%filter(Time==max(Time))%>%
  select(Removed)%>%
  unlist()%>%
  as.numeric()

################################################################################
#Q4
################################################################################

infect_df=q1_df%>%select(Time,Infected)
infect_fun1=approxfun(x=infect_df$Time,y=infect_df$Infected)

infect_over1=function(test_time){
  over_level=0.03
  res=infect_fun1(test_time)-over_level
  return(res)
}

time_over=uniroot(f=infect_over1,lower = 0,upper = 20,
                  tol=10^(-9))
#Time that the plague is over
time_over$root

#Sniff check
infect_fun1(time_over$root)

################################################################################
#Q5
################################################################################

max_infect1=optimize(f=infect_fun1,lower=0,upper = 20,maximum = T)

#infection max time
max_infect1$maximum

#infection max
max_infect1$objective

time_seq=seq(0,time_over$root,length=1000)

q5_data=data.frame(Time=time_seq,Infected=infect_fun1(time_seq))

ggplot(data=q5_data,aes(x=Time,y=Infected))+
  geom_line(color="red",size=1)+
  scale_y_continuous(breaks=seq(0.03,max_infect1$maximum,by=0.025))+
  geom_hline(yintercept = 0.03,linetype="dotted")+
  geom_segment(x=max_infect1$maximum,y=0,xend=max_infect1$maximum,yend=max_infect1$objective,
               linetype="dashed")+
  theme_bw()+
  labs(title = "Infection Curve",
       subtitle = "From the epidemic's beginning to end.")


################################################################################
#Q6
################################################################################

suspect_df=q1_df%>%select(Time,Susceptible)
suspect_fun1=approxfun(x=suspect_df$Time,y=suspect_df$Susceptible)  

#Susceptible share when infections peak 
suspect_fun1(max_infect1$maximum)


################################################################################
#Q8
################################################################################


par_set2=c(beta=1.5,gamma=1.6)


q8_df=SIR_fun(s0=0.87,i0=0.13,r0=0,
              beta=par_set2[1],gamma=par_set2[2],
              dt=0.001,Tend = 20)

plot2=q8_df%>%gather(key="Series",value="share",-1)%>%
  mutate(Series=factor(Series,levels = c("Susceptible","Infected","Removed")))%>%
  ggplot(aes(x=Time,y=share,color=Series))+
  geom_line(size=1)+
  ggthemes::scale_color_fivethirtyeight()+
  theme_bw()+
  scale_y_continuous(breaks = seq(0,1,by=0.1))+
  labs(y="Population Share",title = "SIR Model Simulation",
       subtitle = TeX("($\\beta =1.5$ and $\\gamma =1.6$)"))

################################################################################
#Q9
################################################################################


q8_df%>%filter(Time==max(Time))%>%
  select(Removed)%>%
  unlist()%>%
  as.numeric()

################################################################################
#Bonus
################################################################################

s_max_test=function(beta,gamma,test_val){
  ratio=beta/gamma
  LHS=log(test_val)
  RHS=ratio*(test_val-1)
  res=abs(LHS-RHS)
  return(res)
}



s0_1=optim(f=s_max_test,par=0.5,beta=3,gamma=1.6,lower = 0,upper=1,method = "Brent")$par
s0_1



1/1.6

1/(2*1.6)


#Option 1

q_opt1_df=SIR_fun(s0=0.87,i0=0.13,r0=0,
              beta=1.5,gamma=1.6,
              dt=0.001,Tend = 20)

#Option 2

q_opt2_df=SIR_fun(s0=0.87,i0=0.13,r0=0,
                  beta=3,gamma=3.2,
                  dt=0.001,Tend = 20)

rbind(select(q_opt1_df,Time,Infected)%>%mutate(Policy="Lockdown: 1/2 x Beta"),
      select(q_opt2_df,Time,Infected)%>%mutate(Policy="Medicine: 2 x Gamma"))%>%
  ggplot(aes(x=Time,y=Infected,color=Policy))+
  geom_line(size=1)+
  theme_bw()+
  ggthemes::scale_color_few()+
  labs(title="Infection population over time: Lockdowns vs Medicine")+
  theme(legend.position = "top")

q_opt1_df%>%filter(Time==max(Time))

q_opt2_df%>%filter(Time==max(Time))





trap_int=function(a,b,n,fun1){
  xgrid=seq(a,b,length=n)
  step=(b-a)/n
  res=c(rep(0,n-1))
  for(k in 1:(n-1)){
    res[k]<-fun1(a+k*step)
  }
  res=sum(res)
  res=step*(fun1(a)/2+res+fun1(b)/2)
  return(res)
}

infect_opt1=approxfun(x=q_opt1_df$Time,y=q_opt1_df$Infected)
infect_opt2=approxfun(x=q_opt2_df$Time,y=q_opt2_df$Infected)

area_curve=data.frame(Policy=c("Lockdown: 1/2 x Beta","Medicine: 2 x Gamma"),
                      Lost=c(trap_int(a=0,b=20,n=1000,fun1=infect_opt1),
                             trap_int(a=0,b=20,n=1000,fun1=infect_opt2)))

ggplot(data=area_curve,
       aes(x=Policy,y=Lost,fill=Policy))+
  geom_bar(stat="identity")+
  geom_text(aes(label=round(Lost,2)),vjust=-0.25)+
  theme_bw()+
  ggthemes::scale_fill_few()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "top")+
  labs(y="Total working hours",title="Lost working hours: Lockdowns vs Medicine")




infect_over1=function(x){
  return(infect_opt1(x)-0.03)
}


infect_over2=function(x){
  return(infect_opt2(x)-0.03)
}

uniroot(infect_over1,lower=0,upper=20)$root
uniroot(infect_over2,lower=0,upper=20)$root

uniroot(infect_over1,lower=0,upper=20)$root/uniroot(infect_over2,lower=0,upper=20)$root
