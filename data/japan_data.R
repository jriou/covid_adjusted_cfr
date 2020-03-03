## 25.02.2020

# japan_20200225.csv is a copy-paste of the case list of the website https://jagjapan.maps.arcgis.com/apps/opsdashboard/index.html#/55c22ee976bc42338cb454765a6edf6b
# here, the data is formatted into a df with one case per line.

setwd("C:/Users/counotte/OneDrive - Universitaet Bern/NCOV")

df<-read.csv(file="japan_20200225.csv", stringsAsFactors = FALSE)

df[seq(3,nrow(df),by=19),]$info
df[seq(2,nrow(df),by=19),]$info
df[seq(1,nrow(df),by=19),]$info


df2<-data.frame('id'=df[seq(1,nrow(df),by=19),]$info,'age'=df[seq(2,nrow(df),by=19),]$info,
                   'sex'=df[seq(3,nrow(df),by=19),]$info,'date'=as.Date(df[seq(4,nrow(df),by=19),]$info, format="%m/%d/%Y"))


library(ggplot2)
# check age distr
ggplot(df2,aes(x=age, fill=sex))+geom_bar()

ggplot(df2,aes(x=date))+geom_bar()
                