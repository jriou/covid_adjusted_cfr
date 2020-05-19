library("readxl")
library("dplyr")
current_wd=getwd()
# setwd("C:/Users/ahauser/Dropbox/covid_adjusted_cfr")
# Contact matrix: data/contact_matrix_prems2017/MUestimates_all_locations_1.xlsx
# and data/contact_matrix_prems2017/MUestimates_all_locations_2.xlsx
# from: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005697#sec020
# Age distribution: data/age_distribution.csv, from (used to transform contact matrix): http://data.un.org/Data.aspx?d=POP&f=tableCode%3A22

#function to symmetrize contact matrix, used in get_sym_contact_matrix10/20/5_85
symmetrize = function(mat,dist){
  d=dim(mat)[1]
  m=mat*matrix(rep(dist,d),nrow=d,byrow=FALSE)
  return((m+t(m))/(2*matrix(rep(dist,d),nrow=d,byrow=FALSE)))
}

#10-year classes
age_class=function(x,min_age,age_range,max_age){
  age_lim=seq(min_age,max_age,age_range)
  return(sapply(as.list(x),function(x) sum(age_lim<=x)))
}
age_class222=function(x,age_range,max_age){
  age_lim=seq(0,max_age,age_range)
  return(sapply(as.list(x),function(x) sum(age_lim<=x)))
}

#------------------------------------------------------------------------------------------------------
#10-year and 20-year classes: symmetric matrices
get_sym_contact_matrix10=function(country,year_dist){
  #load age distribution, 5-year group (as contact matrix is split in 5-year groups)
  excel_age_dist=ifelse(country=="China","data/contact_matrix/age_distribution_china.csv","data/contact_matrix/age_distribution.csv")
  
  age_dist = as.data.frame(read.csv(excel_age_dist)) %>% 
    filter(Country.or.Area==country) %>%
    filter(Year == year_dist) %>%
    filter(Age %in% as.character(0:120)) %>%
    filter(Sex=="Both Sexes") %>%
    mutate(age_group=age_class(as.numeric(levels(Age))[Age],5,5,80)) %>% 
    group_by(age_group) %>%
    summarise(n=sum(Value)) %>%
    pull(n)
  age_dist5 = age_dist/sum(age_dist)
  
  #load age distribution, 20-year group
  age_dist = as.data.frame(read.csv(excel_age_dist)) %>%
    filter(Country.or.Area==country) %>%
    filter(Year == year_dist) %>%
    filter(Age %in% as.character(0:120)) %>%
    filter(Sex=="Both Sexes") %>%
    mutate(age_group=age_class(as.numeric(levels(Age))[Age],10,10,80)) %>% 
    group_by(age_group) %>%
    summarise(n=sum(Value)) %>%
    pull(n)
  age_dist10 = age_dist/sum(age_dist)
  
  #load contact matrix
  country_excel1<-excel_sheets("data/contact_matrix/contact_matrix_prems2017/MUestimates_all_locations_1.xlsx")
  country_excel2<-excel_sheets("data/contact_matrix/contact_matrix_prems2017/MUestimates_all_locations_2.xlsx")
  if(country %in% country_excel1){
    contact <- as.data.frame(((read_excel("data/contact_matrix/contact_matrix_prems2017/MUestimates_all_locations_1.xlsx",sheet=country,col_names=TRUE))))
  }else if(country %in% country_excel2){
    contact <- as.data.frame(((read_excel("data/contact_matrix/contact_matrix_prems2017/MUestimates_all_locations_2.xlsx",sheet=country,col_names=FALSE))))
  }else{
    stop("--------------- No contact matrix for this country.")
  }
  colnames(contact)=paste("X",seq(0,75,5),sep = "")
  
  #divide the >=75 group into 75-79 and >=80 groups
  tmp = contact[,16]*age_dist5[16]/sum(age_dist5[16:17])
  tmp2 = contact[,16]*age_dist5[17]/sum(age_dist5[16:17])
  contact[,16]=tmp
  contact=cbind(contact,tmp2)
  tmp=c(contact[,16]*age_dist5[1:16]/age_dist5[17],contact[16,16])
  contact=rbind(contact,tmp)
  #symmetrize matrix: https://cran.r-project.org/web/packages/socialmixr/vignettes/introduction.html
  contact<-symmetrize(contact,age_dist5)
  rownames(contact)=NULL
  colnames(contact)=paste("X",seq(0,80,5),sep = "")
  
  #transform ino 9 age groups
  contact2= contact%>%
    transmute(X0+X5,X10+X15,X20+X25,X30+X35,X40+X45,X50+X55,X60+X65,X70+X75,X80)  %>%
    mutate(age2=c(0,0,10,10,20,20,30,30,40,40,50,50,60,60,70,70,80),
           age_dist=age_dist5)%>%
    group_by(age2) %>%
    mutate(w=age_dist/sum(age_dist)) %>%
    summarise(sum(`X0 + X5`*w),sum(`X10 + X15`*w),sum(`X20 + X25`*w),sum(`X30 + X35`*w),
              sum(`X40 + X45`*w),sum(`X50 + X55`*w),sum(`X60 + X65`*w),sum(`X70 + X75`*w),sum(`X80`*w))  %>%
    select(-age2) %>%
    as.matrix(.)
  rownames(contact2)=NULL
  colnames(contact2)=NULL
  
  print(paste("c(",paste(t(contact2),sep="' '", collapse=", "),")",sep=" "))
  return(contact2)
}
get_sym_contact_matrix20=function(country){
  #load age distribution, 5-year group (as contact matrix is split in 5-year groups)
  excel_age_dist=ifelse(country=="China","data/contact_matrix/age_distribution_china.csv","data/contact_matrix/age_distribution.csv")
  
  age_dist = as.data.frame(read.csv(excel_age_dist)) %>%
    filter(Country.or.Area==country) %>%
    filter(Year == max(as.numeric(levels(Year))[Year], na.rm=TRUE)) %>%
    filter(Source.Year == max(as.numeric(Source.Year), na.rm=TRUE)) %>%
    filter(Age %in% as.character(0:120)) %>%
    filter(Sex=="Both Sexes") %>%
    mutate(age_group=age_class(as.numeric(levels(Age))[Age],5,5,80)) %>% 
    group_by(age_group) %>%
    summarise(n=sum(Value)) %>%
    pull(n)
  age_dist5 = age_dist/sum(age_dist)
  
  #load age distribution, 20-year group
  age_dist = as.data.frame(read.csv(excel_age_dist)) %>%
    filter(Country.or.Area==country) %>%
    filter(Year == max(as.numeric(levels(Year))[Year], na.rm=TRUE)) %>%
    filter(Source.Year == max(as.numeric(Source.Year), na.rm=TRUE)) %>%
    filter(Age %in% as.character(0:120)) %>%
    filter(Sex=="Both Sexes") %>%
    mutate(age_group=age_class(as.numeric(levels(Age))[Age],20,20,80)) %>% 
    group_by(age_group) %>%
    summarise(n=sum(Value)) %>%
    pull(n)
  age_dist20 = age_dist/sum(age_dist)
  
  #load contact matrix
  country_excel1<-excel_sheets("data/contact_matrix/contact_matrix_prems2017/MUestimates_all_locations_1.xlsx")
  country_excel2<-excel_sheets("data/contact_matrix/contact_matrix_prems2017/MUestimates_all_locations_2.xlsx")
  if(country %in% country_excel1){
    contact <- as.data.frame(((read_excel("data/contact_matrix/contact_matrix_prems2017/MUestimates_all_locations_1.xlsx",sheet=country,col_names=TRUE))))
  }else if(country %in% country_excel2){
    contact <- as.data.frame(((read_excel("data/contact_matrix/contact_matrix_prems2017/MUestimates_all_locations_2.xlsx",sheet=country,col_names=FALSE))))
  }else{
    stop("--------------- No contact matrix for this country.")
  }
  colnames(contact)=paste("X",seq(0,75,5),sep = "")
  print("contact matrix before transformation")
  print(contact)
  
  #divide the >=75 group into 75-79 and >=80 groups
  tmp = contact[,16]*age_dist5[16]/sum(age_dist5[16:17])
  tmp2 = contact[,16]*age_dist5[17]/sum(age_dist5[16:17])
  contact[,16]=tmp
  contact=cbind(contact,tmp2)
  tmp=c(contact[,16]*age_dist5[1:16]/age_dist5[17],contact[16,16])
  contact=rbind(contact,tmp)
  #symmetrize matrix: https://cran.r-project.org/web/packages/socialmixr/vignettes/introduction.html
  contact<-symmetrize(contact,age_dist5)
  rownames(contact)=NULL
  colnames(contact)=paste("X",seq(0,80,5),sep = "")
  
  #transform ino 5 age groups
  contact2= contact%>%
    transmute(g1=X0+X5+X10+X15,g2=X20+X25+X30+X35,g3=X40+X45+X50+X55,g4=X60+X65+X70+X75,g5=X80)  %>%
    mutate(age2=c(rep(1:4,each=4),5),
           age_dist=age_dist5) %>%
    group_by(age2) %>%
    mutate(w=age_dist/sum(age_dist)) %>% summarise(sum(g1 *w),sum(g2 *w),sum(g3 *w),sum(g4 *w),sum(g5 *w)) %>%
    select(-age2) %>%
    as.matrix(.)
  rownames(contact2)=NULL
  colnames(contact2)=NULL
  
  print(age_dist20)
  print(paste("c(",paste(t(contact2),sep="' '", collapse=", "),")",sep=" "))
  return(contact2)
}
#Matrix with 9 age groups: 0-4,5-14,  ,>75 for Austria
get_sym_contact_matrix5_75=function(country,year_dist){
  #load age distribution, 5-year group (as contact matrix is split in 5-year groups)
  excel_age_dist=ifelse(country=="China","data/contact_matrix/age_distribution_china.csv","data/contact_matrix/age_distribution.csv")
  
  age_dist = as.data.frame(read.csv(excel_age_dist)) %>%
    filter(Country.or.Area==country) %>%
    filter(Year == year_dist) %>%
    filter(Age %in% as.character(0:120)) %>%
    filter(Sex=="Both Sexes") %>%
    mutate(age_group=age_class(as.numeric(levels(Age))[Age],5,5,75)) %>% 
    group_by(age_group) %>%
    summarise(n=sum(Value)) %>%
    pull(n)
  age_dist5 = age_dist/sum(age_dist)
  
  #load age distribution, 10-year group
  age_dist = as.data.frame(read.csv(excel_age_dist)) %>%
    filter(Country.or.Area==country) %>%
    filter(Year == year_dist) %>%
    filter(Age %in% as.character(0:120)) %>%
    filter(Sex=="Both Sexes") %>%
    mutate(age_group=age_class(as.numeric(levels(Age))[Age],5,10,75)) %>% 
    group_by(age_group) %>%
    summarise(n=sum(Value)) %>%
    pull(n)
  age_dist10 = age_dist/sum(age_dist)
  
  #load contact matrix
  country_excel1<-excel_sheets("data/contact_matrix/contact_matrix_prems2017/MUestimates_all_locations_1.xlsx")
  country_excel2<-excel_sheets("data/contact_matrix/contact_matrix_prems2017/MUestimates_all_locations_2.xlsx")
  if(country %in% country_excel1){
    contact <- as.data.frame(((read_excel("data/contact_matrix/contact_matrix_prems2017/MUestimates_all_locations_1.xlsx",sheet=country,col_names=TRUE))))
  }else if(country %in% country_excel2){
    contact <- as.data.frame(((read_excel("data/contact_matrix/contact_matrix_prems2017/MUestimates_all_locations_2.xlsx",sheet=country,col_names=FALSE))))
  }else{
    stop("--------------- No contact matrix for this country.")
  }
  colnames(contact)=paste("X",seq(0,75,5),sep = "")
  
  #transform into 9 age groups
  contact2= contact%>%
    transmute(X0,X5+X10,X15+X20,X25+X30,X35+X40,X45+X50,X55+X60,X65+X70,X75)  %>%
    mutate(age2=c(0,5,5,15,15,25,25,35,35,45,45,55,55,65,65,75),
           age_dist=age_dist5)%>%
    group_by(age2) %>%
    mutate(w=age_dist/sum(age_dist)) %>%
    summarise(sum(`X0`*w),sum(`X5 + X10`*w),sum(`X15 + X20`*w),sum(`X25 + X30`*w),
              sum(`X35 + X40`*w),sum(`X45 + X50`*w),sum(`X55 + X60`*w),sum(`X65 + X70`*w),sum(`X75`*w))  %>%
    select(-age2) %>%
    as.matrix(.)
  rownames(contact2)=NULL
  colnames(contact2)=NULL
  
  print(paste("c(",paste(t(contact2),sep="' '", collapse=", "),")",sep=" "))
  return(contact2)
}

# m<-get_sym_contact_matrix20("Switzerland")
# m<-get_sym_contact_matrix10("Italy")
# m<-get_sym_contact_matrix10("Sweden")
# m<-get_sym_contact_matrix10("China")

#------------------------------------------------------------------------------------------------------
#3 age classes, 0-19, 20-59, >=60
age_class3=function(x){
  age_lim=c(0,20,60)
  return(sapply(as.list(x),function(x) sum(age_lim<=x)))
}
get_contact_matrix3=function(country){
  #Age distribution, 5-year group
  excel_age_dist=ifelse(country=="China","data/contact_matrix/age_distribution_china.csv","data/contact_matrix/age_distribution.csv")
  age_dist = as.data.frame(read.csv(excel_age_dist)) %>%
    filter(Country.or.Area==country) %>%
    filter(Year == max(as.numeric(levels(Year))[Year], na.rm=TRUE)) %>%
    filter(Source.Year == max(as.numeric(Source.Year), na.rm=TRUE)) %>%
    filter(Age %in% as.character(0:120)) %>%
    filter(Sex=="Both Sexes") %>%
    mutate(age_group=age_class(as.numeric(levels(Age))[Age],5,5,75)) %>% 
    group_by(age_group) %>%
    summarise(n=sum(Value)) %>%
    pull(n)
  age_dist5 = age_dist/sum(age_dist)
  
  #Age distribution, 10-year group
  age_dist = as.data.frame(read.csv(excel_age_dist)) %>%
    filter(Country.or.Area==country) %>%
    filter(Year == max(as.numeric(levels(Year))[Year], na.rm=TRUE)) %>%
    filter(Source.Year == max(as.numeric(Source.Year), na.rm=TRUE)) %>%
    filter(Age %in% as.character(0:120)) %>%
    filter(Sex=="Both Sexes") %>%
    mutate(age_group=age_class3(as.numeric(levels(Age))[Age])) %>% 
    group_by(age_group) %>%
    summarise(n=sum(Value)) %>%
    pull(n)
  age_dist3 = age_dist/sum(age_dist)
  
  
  
  #Contact matrix
  country_excel1<-excel_sheets("data/contact_matrix/contact_matrix_prems2017/MUestimates_all_locations_1.xlsx")
  country_excel2<-excel_sheets("data/contact_matrix/contact_matrix_prems2017/MUestimates_all_locations_2.xlsx")
  if(country %in% country_excel1){
    contact <- as.data.frame(((read_excel("data/contact_matrix/contact_matrix_prems2017/MUestimates_all_locations_1.xlsx",sheet=country,col_names=TRUE))))
  }else if(country %in% country_excel2){
    contact <- as.data.frame(((read_excel("data/contact_matrix/contact_matrix_prems2017/MUestimates_all_locations_2.xlsx",sheet=country,col_names=FALSE))))
  }else{
    stop("--------------- No contact matrix for this country.")
  }
  colnames(contact)=paste("X",seq(0,75,5),sep = "")
  
  contact2= contact%>%
    transmute(g1=X0+X5+X10+X15,g2=X20+X25+X30+X35+X40+X45+X50+X55,g3=X60+X65+X70+X75)  %>%
    mutate(age2=c(1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3),
           age_dist=age_dist5) %>%
    group_by(age2) %>%
    mutate(w=age_dist/sum(age_dist)) %>% summarise(sum(g1 *w),sum(g2 *w),sum(g3 *w)) %>%
    select(-age2) %>%
    as.matrix(.)
  rownames(contact2)=NULL
  colnames(contact2)=NULL
  
  print("age distribution")
  print(paste("c(",paste(t(age_dist3),sep="' '", collapse=", "),")",sep=" "))
  print("contact matrix")
  print(paste("c(",paste(t(contact2),sep="' '", collapse=", "),")",sep=" "))
  return(contact2)
}

# m<-get_contact_matrix3("Switzerland")

setwd(current_wd)
#####################################################################################################################
#Do not delete this, it might be used later on
#Age distribution, 5-year group
# age_dist = as.data.frame(read.csv("data/age_distribution.csv")) %>%
#   filter(Country.or.Area==country) %>%
#   filter(Sex=="Both Sexes")
# age_dist_t1=age_dist %>%
#   filter(Age %in% paste(seq(0,95,5),"-", seq(4,99,5)))


# get_contact_matrix_old=function(country){
#   #Age distribution, 5-year group
#   excel_age_dist=ifelse(country=="China","data/contact_matrix/age_distribution_china.csv","data/contact_matrix/age_distribution.csv")
#   
#   age_dist = as.data.frame(read.csv(excel_age_dist)) %>%
#     filter(Country.or.Area==country) %>%
#     filter(Year == max(as.numeric(levels(Year))[Year], na.rm=TRUE)) %>%
#     filter(Source.Year == max(as.numeric(Source.Year), na.rm=TRUE)) %>%
#     filter(Age %in% as.character(0:120)) %>%
#     filter(Sex=="Both Sexes") %>%
#     mutate(age_group=age_class(as.numeric(levels(Age))[Age],5,5,75)) %>% 
#     group_by(age_group) %>%
#     summarise(n=sum(Value)) %>%
#     pull(n)
#   age_dist5 = age_dist/sum(age_dist)
#   
#   #Age distribution, 10-year group
#   age_dist = as.data.frame(read.csv(excel_age_dist)) %>%
#     filter(Country.or.Area==country) %>%
#     filter(Year == max(as.numeric(levels(Year))[Year], na.rm=TRUE)) %>%
#     filter(Source.Year == max(as.numeric(Source.Year), na.rm=TRUE)) %>%
#     filter(Age %in% as.character(0:120)) %>%
#     filter(Sex=="Both Sexes") %>%
#     mutate(age_group=age_class(as.numeric(levels(Age))[Age],10,10,80)) %>% 
#     group_by(age_group) %>%
#     summarise(n=sum(Value)) %>%
#     pull(n)
#   age_dist10 = age_dist/sum(age_dist)
#   
#   
#   
#   #Contact matrix
#   country_excel1<-excel_sheets("data/contact_matrix/contact_matrix_prems2017/MUestimates_all_locations_1.xlsx")
#   country_excel2<-excel_sheets("data/contact_matrix/contact_matrix_prems2017/MUestimates_all_locations_2.xlsx")
#   if(country %in% country_excel1){
#     contact <- as.data.frame(((read_excel("data/contact_matrix/contact_matrix_prems2017/MUestimates_all_locations_1.xlsx",sheet=country,col_names=TRUE))))
#   }else if(country %in% country_excel2){
#     contact <- as.data.frame(((read_excel("data/contact_matrix/contact_matrix_prems2017/MUestimates_all_locations_2.xlsx",sheet=country,col_names=FALSE))))
#   }else{
#     stop("--------------- No contact matrix for this country.")
#   }
#   colnames(contact)=paste("X",seq(0,75,5),sep = "")
#   print("contact matrix before transformation")
#   print(contact)
#   
#   contact2= contact%>%
#     transmute(X0+X5,X10+X15,X20+X25,X30+X35,X40+X45,X50+X55,X60+X65,X70+X75)  %>%
#     mutate(age2=c(0,0,10,10,20,20,30,30,40,40,50,50,60,60,70,70),
#            age_dist=age_dist5)%>%
#     group_by(age2) %>%
#     mutate(w=age_dist/sum(age_dist)) %>%
#     summarise(sum(`X0 + X5`*w),sum(`X10 + X15`*w),sum(`X20 + X25`*w),sum(`X30 + X35`*w),
#               sum(`X40 + X45`*w),sum(`X50 + X55`*w),sum(`X60 + X65`*w),sum(`X70 + X75`*w))  %>%
#     select(-age2) %>%
#     as.matrix(.)
#   rownames(contact2)=NULL
#   colnames(contact2)=NULL
#   
#   tmp = contact2[,8]*age_dist10[8]/sum(age_dist10[8:9])
#   tmp2 = contact2[,8]*age_dist10[9]/sum(age_dist10[8:9])
#   contact2[,8]=tmp
#   contact2=cbind(contact2,tmp2)
#   tmp=c(contact2[,9]*age_dist10[1:8]/age_dist10[9],contact2[8,8])
#   contact2=rbind(contact2,tmp)
#   rownames(contact2)=NULL
#   colnames(contact2)=NULL
#   
#   print(paste("c(",paste(t(contact2),sep="' '", collapse=", "),")",sep=" "))
#   return(contact2)
# }