rm(list=ls())
df<-read.csv("https://raw.githubusercontent.com/jihoo-kim/Coronavirus-Dataset/master/patient.csv")
library(cowplot)
library(ggplot2)

df$confirmed_date<-as.Date(df$confirmed_date,format="%Y-%m-%d")
df$deceased_date<-as.Date(df$deceased_date,format="%Y-%m-%d")


ggplot()+
  geom_bar(data=df,aes(x=confirmed_date), fill="blue")
  #geom_bar(data=df,aes(x=deceased_date),fill="red")

ggplot()+
  #geom_bar(data=df,aes(x=confirmed_date))+
  geom_bar(data=df,aes(x=deceased_date),fill="red")

#calculate age:
df$age<-2020-df$birth_year
ages<-c(0,10,20,30,40,50,60,70,80,90,100)
df$ageband<-cut(df$age, breaks=ages)
hist(df$age)


#table(df$ageband)
inc=unname(as.data.frame(table(df$ageband))[2])
mort=unname(as.data.frame(table(df[!is.na(df$deceased_date),]$ageband))[2])
age_distr=data.frame('ages'=ages[1:10],'inc'=inc,'mort'=mort)

write.csv(age_distr,file="Korea_ages.csv", row.names = FALSE)

library(dplyr)
cases<-df%>%group_by(confirmed_date)%>%summarise(cases = n())
mort<-df[!is.na(df$deceased_date),]%>%group_by(deceased_date)%>%summarise(mort = n())
mort$confirmed_date<-mort$deceased_date

df_cases=merge(cases, mort, by="confirmed_date",all = TRUE)
df_cases=df_cases[,-3]
colnames(df_cases[1])="date"

write.csv(df_cases,file="Korea_cases.csv", row.names = FALSE)
