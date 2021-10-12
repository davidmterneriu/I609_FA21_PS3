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


