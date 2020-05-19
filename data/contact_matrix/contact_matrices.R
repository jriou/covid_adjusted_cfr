#Matrix of contact between age classes

#https://cran.r-project.org/web/packages/socialmixr/vignettes/introduction.html
library(socialmixr)
library(dplyr)
library(tidyr)
library(ggplot2)

list_surveys()

#China (Shanghai)
shanghai_survey <- get_survey("https://doi.org/10.5281/zenodo.3366396")
m <- contact_matrix(survey = shanghai_survey, age.limits=c(0,10,20,30,40,50,60,70,80),symmetric = TRUE)
data.frame(contacts=as.numeric(m$matrix)) %>%
  tbl_df() %>%
  mutate(age1=rep(1:9,9),age2=rep(1:9,each=9)) %>%
  ggplot() +
  geom_tile(aes(x=age2,y=age1,fill=contacts))
contact = structure(m$matrix,dim=c(9,9))
apply(contact,1,sum)
sum(contact)
paste(t(contact),sep="' '", collapse=", ")


#China (Shanghai)
shanghai_survey <- get_survey("https://doi.org/10.5281/zenodo.3366396")
m <- contact_matrix(survey = shanghai_survey, age.limits=seq(0,75,5),symmetric = TRUE)
data.frame(contacts=as.numeric(m$matrix)) %>%
  tbl_df() %>%
  mutate(age1=rep(1:16,16),age2=rep(1:16,each=16)) %>%
  ggplot() +
  geom_tile(aes(x=age2,y=age1,fill=contacts))
contact = structure(m$matrix,dim=c(16,16))
apply(contact,1,sum)


#Italy, 9 age classes
italy_survey <- get_survey("https://doi.org/10.5281/zenodo.1043437")
m <- contact_matrix(survey = italy_survey, age.limits=c(0,10,20,30,40,50,60,70,80),symmetric = TRUE)
data.frame(contacts=as.numeric(m$matrix)) %>%
  tbl_df() %>%
  mutate(age1=rep(1:9,9),age2=rep(1:9,each=9)) %>%
  ggplot() +
  geom_tile(aes(x=age2,y=age1,fill=contacts))
contact = structure(m$matrix,dim=c(9,9))
apply(contact,1,mean)
sum(contact)
paste(t(contact),sep="' '", collapse=", ")

#Italy, 16 age classes
italy_survey <- get_survey("https://doi.org/10.5281/zenodo.1043437")
m <- contact_matrix(survey = italy_survey, age.limits=seq(0,75,5),symmetric = TRUE)
data.frame(contacts=as.numeric(m$matrix)) %>%
  tbl_df() %>%
  mutate(age1=rep(1:16,16),age2=rep(1:16,each=16)) %>%
  ggplot() +
  geom_tile(aes(x=age2,y=age1,fill=contacts))

#Italy ()
italy_survey <- get_survey("https://doi.org/10.5281/zenodo.1043437")
m <- contact_matrix(survey = italy_survey, age.limits=c(0,10,20,30,40,50,60,70,80),symmetric = TRUE)
m <- contact_matrix(survey = italy_survey, age.limits=c(0,5,20,65),symmetric = TRUE)
contact = structure(m$matrix,dim=c(9,9))
apply(contact,1,sum)
sum(contact)
paste(t(contact),sep="' '", collapse=", ")
contact_it <-c(5.13567073170732, 1.17274819632136, 0.982359525171638, 2.21715890088845, 1.29666356906914, 0.828866413937242, 0.528700773224482, 0.232116187961884, 0.0975205061876398, 1.01399087153423, 10.420788530466, 1.5084165224448, 1.46323525034693, 2.30050630727188, 1.0455742822567, 0.396916593664865, 0.276112578159939, 0.0867321859134207, 0.787940961549209, 1.39931415327149, 4.91448118586089, 2.39551550152373, 2.08291844616138, 1.67353143324194, 0.652483430981848, 0.263165822550241, 0.107498717856296, 1.53454251726848, 1.17129688889679, 2.06708280469829, 3.91165644171779, 2.74588910732349, 1.66499320847473, 1.02145416818956, 0.371633336270256, 0.112670158106901, 0.857264438638371, 1.7590640625625, 1.71686658407219, 2.62294018855816, 3.45916114790287, 1.87635185962704, 0.862205884832066, 0.523958801433231, 0.205791955532149, 0.646645383952458, 0.943424739130445, 1.62776721065554, 1.87677409215498, 2.21415705015835, 2.5920177383592, 1.10525460534109, 0.472961105423521, 0.282448363507455, 0.504954014454259, 0.438441714821823, 0.77694120330432, 1.40954408148402, 1.24556204828388, 1.35307720400585, 1.70385674931129, 0.812686154912104, 0.270111273681845, 0.305701280434649, 0.420580126969344, 0.432113761275257, 0.707170907986224, 1.04376196943771, 0.798427737704416, 1.12065725135372, 1.33035714285714, 0.322575366839763, 0.237578345845701, 0.24437789962337, 0.326505855457376, 0.396586297530862, 0.758318763302674, 0.881999483055259, 0.688988121391528, 0.596692087603768, 0.292682926829268)
sum(contact_it)

contacts = matrix(contact_it,nrow=9,byrow=T) %>%
  data.frame() %>%
  tbl_df()
names(contacts) = c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+")
contacts$i = c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+")
h4 = gather(contacts,"j","c",1:9) %>%
  ggplot() +
  geom_tile(aes(x=i,y=j,fill=c)) +
  scale_fill_gradient(low="darkblue",high="darkgoldenrod1") +
  labs(x="Age group",y="Age group",fill="Contacts")
h4



#############################################################################
#Korea
#Age distribution by 5 years
library("readxl")
age_dist = read_excel("data/age_structure.xlsx") %>%
  filter(country=="Korea, Republic of") %>%
  gather("age","n",3:23) %>%
  mutate(n=as.numeric(n),age2=c(seq(0,70,5),rep(75,6))) %>%
  group_by(age2) %>%
  summarise(n=sum(n)) %>%
  pull(n)
age_dist_korea = age_dist/sum(age_dist)
#Age distribution by 10 years
age_dist = read_excel("data/age_structure.xlsx") %>%
  filter(country=="Korea, Republic of") %>%
  gather("age","n",3:23) %>%
  mutate(n=as.numeric(n),age2=c(0,0,10,10,20,20,30,30,40,40,50,50,60,60,70,70,80,80,80,80,80)) %>%
  group_by(age2) %>%
  summarise(n=sum(n)) %>%
  pull(n)
age_dist_korea10 = age_dist/sum(age_dist)


contact_korea=as.data.frame(((read.csv("C:/Users/ahauser/Documents/GitHub/COVID_age/data/contact_matrix_prem_et_al.csv"))[,-1]))

contact_korea2= contact_korea%>%
  transmute(X0+X5,X10+X15,X20+X25,X30+X35,X40+X45,X50+X55,X60+X65,X70+X75)  %>%
  mutate(age2=c(0,0,10,10,20,20,30,30,40,40,50,50,60,60,70,70),
         age_dist_korea=age_dist_korea)%>%
  group_by(age2) %>%
  mutate(w=age_dist_korea/sum(age_dist_korea)) %>%
  summarise(sum(`X0 + X5`*w),sum(`X10 + X15`*w),sum(`X20 + X25`*w),sum(`X30 + X35`*w),
            sum(`X40 + X45`*w),sum(`X50 + X55`*w),sum(`X60 + X65`*w),sum(`X70 + X75`*w))  %>%
  select(-age2) %>%
  as.matrix(.)
rownames(contact_korea2)=NULL
colnames(contact_korea2)=NULL
contact_korea2

tmp = contact_korea2[,8]*age_dist_korea10[8]/sum(age_dist_korea10[8:9])
tmp2 = contact_korea2[,8]*age_dist_korea10[9]/sum(age_dist_korea10[8:9])
contact_korea2[,8]=tmp
contact_korea2=cbind(contact_korea2,tmp2)
tmp=c(contact_korea2[,9]*age_dist_korea10[1:8]/age_dist_korea10[9],contact_korea2[8,8])
contact_korea2=rbind(contact_korea2,tmp)
rownames(contact_korea2)=NULL
colnames(contact_korea2)=NULL
contact_korea2
paste(t(contact_korea2),sep="' '", collapse=", ")
contact_korea=c(3.521460976678, 0.807606635584979, 0.58017221932137, 1.45700400214389, 0.935832062414259, 0.404107489354786, 0.205352168108353, 0.0486140785410915, 0.0363597064909721, 1.13965890032514, 13.3263555794776, 1.37335196151099, 1.31217092938662, 2.01921605371974, 0.727844341372875, 0.157232043396022, 0.0541469767195827, 0.0404979018420231, 0.524125317641722, 1.77138293677653, 5.37342282955803, 2.6119530994001, 2.12498433097382, 1.49061176889747, 0.25082030298294, 0.0475140570215849, 0.0355369724027551, 1.73548563162998, 1.37863424962868, 2.33028799573755, 5.01693062997484, 3.00626105376023, 1.51976579426414, 0.479233855478861, 0.0777661367602361, 0.0581632727060529, 1.06279703666294, 2.49460796812014, 1.96780746460473, 3.08801652748333, 4.3284390359967, 1.70659785329017, 0.357728965065922, 0.0984010009252109, 0.0735966127391337, 1.01852479190416, 2.2415190271648, 2.27873878669667, 2.47625738362958, 2.98171697267116, 3.09444915679377, 0.606931880875396, 0.117623609620315, 0.0879736909666683, 0.6337638954766, 0.609156467872957, 0.857177966163359, 1.51627761423349, 1.18954954052422, 1.10638728069481, 1.23464881046495, 0.159781597736603, 0.119504723131815, 0.452002800912675, 0.708904909045733, 0.347302041274133, 0.710632862584685, 0.964118814556132, 0.721873750011454, 0.624025387196449, 0.426281556119421, 0.318826824001527, 0.0432872868178882, 0.0530447730927251, 0.0501000809302181, 0.0959421918690952, 0.158160905518649, 0.187119413095164, 0.200487698765865, 0.426281556119421, 0.426281556119421)

#############################################################################
#Spain
setwd("C:/Users/ahauser/Documents/GitHub/COVID_age/")
library("readxl")

#Age distribution, 1-year group
#http://data.un.org/Data.aspx?d=POP&f=tableCode%3A22
age_dist = as.data.frame(read.csv("data/age_dist_spain.csv")) %>%
  filter(Age %in% as.character(0:100)) %>%
  filter(Sex=="Both Sexes") %>%
  pull(Value)
age_dist_spain1 = age_dist/sum(age_dist)

#Age distribution, 10-year group
age_dist = as.data.frame(read.csv("data/age_dist_spain.csv")) %>%
  filter(Age %in% as.character(0:100)) %>%
  filter(Sex=="Both Sexes") %>%
  mutate(age_group=c(rep(1:8,each=10),rep(9,20)))%>%
  group_by(age_group)%>%
  summarise(n=sum(Value))%>%
  pull(n)
age_dist_spain10 = age_dist/sum(age_dist)


#Contact matrix, symmetric
#https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002673
contact_spain=as.data.frame(((read_excel("data/contact_matrix_spain.xls",sheet="Spain",col_names=FALSE))))
contact_spain = contact_spain/age_dist_spain1 #adjust for population size
mat=contact_spain%>%
  transmute(age1=rowSums(.[1:10]),
            age2=rowSums(.[11:20]),
            age3=rowSums(.[21:30]),
            age4=rowSums(.[31:40]),
            age5=rowSums(.[41:50]),
            age6=rowSums(.[51:60]),
            age7=rowSums(.[61:70]),
            age8=rowSums(.[71:80]),
            age9=rowSums(.[81:100]))%>%
  mutate(age_group=c(rep(1:8,each=10),rep(9,20)),
         w=age_dist_spain1)%>%
  group_by(age_group)%>%
  summarise(age1=sum(age1*w)/sum(w),
            age2=sum(age2*w)/sum(w),
            age3=sum(age3*w)/sum(w),
            age4=sum(age4*w)/sum(w),
            age5=sum(age5*w)/sum(w),
            age6=sum(age6*w)/sum(w),
            age7=sum(age7*w)/sum(w),
            age8=sum(age8*w)/sum(w),
            age9=sum(age9*w)/sum(w))%>%
            select(-age_group) %>%
            as.matrix(.)
contact_spain=mat/min(mat)
paste(t(contact_spain),sep="' '", collapse=", ")
contact_spain=c(49.2269063867702, 13.7408907975891, 9.01127205786125, 13.9439087568032, 9.35472847081987, 3.69082917670957, 2.33029827097489, 1.87517583377893, 1.05665063349587, 13.2101405795172, 48.1916365195594, 11.8836388255581, 13.5091202160572, 12.0310032302973, 7.59031475638765, 2.48400285189405, 1.81298159381375, 1, 8.55139520544766, 11.7302635358618, 30.4026349879839, 22.9255971192851, 20.0298719784758, 14.0577361513123, 7.25554787580836, 2.78591726425901, 1.3476662550505, 9.82234685847999, 9.89840568685278, 17.0176863065936, 21.8048759561126, 16.4288705595471, 11.2894944364076, 6.66203415499145, 3.75012858936071, 1.22525889922359, 5.4749545362594, 7.32417349294473, 12.3531190849477, 13.6498005232549, 11.3913519522699, 7.53907516437613, 4.31201812797397, 2.72220210311943, 1.35781295941957, 2.45983549971467, 5.26198195285647, 9.87294303380208, 10.681343444786, 8.58520879931757, 6.37126444547917, 4.22633615791024, 2.53482550207287, 1.42464670350286, 2.06187045956757, 2.28617463827454, 6.76502178385478, 8.368083100646, 6.51900194556572, 5.61089034495126, 4.8314056938546, 4.19654037723667, 1.62110046844043, 2.31769301477266, 2.3308523213618, 3.6285344670085, 6.5800522156975, 5.74890257285618, 4.70088754941817, 5.86213015742467, 5.57559531432068, 3.28317188581957, 1.66389491252605, 1.63795484494177, 2.23627883623868, 2.73899732131152, 3.653297961978, 3.36604202555953, 2.88505796904813, 4.1828682142039, 2.97657103609473)

#############################################################################
#Spain
#Age distribution, 5-year group
#http://data.un.org/Data.aspx?d=POP&f=tableCode%3A22
age_dist = as.data.frame(read.csv("data/age_dist_spain.csv")) %>%
  filter(Age %in% as.character(0:100)) %>%
  filter(Sex=="Both Sexes") %>%
  mutate(age_group=c(rep(seq(0,70,5),each=5),rep(75,25))) %>%
  group_by(age_group) %>%
  summarise(n=sum(Value)) %>%
  pull(n)
age_dist_spain5 = age_dist/sum(age_dist)

#Age distribution, 10-year group
age_dist = as.data.frame(read.csv("data/age_dist_spain.csv")) %>%
  filter(Age %in% as.character(0:100)) %>%
  filter(Sex=="Both Sexes") %>%
  mutate(age_group=c(rep(1:8,each=10),rep(9,20)))%>%
  group_by(age_group)%>%
  summarise(n=sum(Value))%>%
  pull(n)
age_dist_spain10 = age_dist/sum(age_dist)


contact_spain=as.data.frame(read_excel("C:/Users/ahauser/Documents/GitHub/COVID_age/data/age_contact_prem_2.xlsx",sheet="Spain",col_names = paste("X",seq(0,75,5),sep="")))

contact_spain2= contact_spain%>%
  transmute(X0+X5,X10+X15,X20+X25,X30+X35,X40+X45,X50+X55,X60+X65,X70+X75)  %>%
  mutate(age2=c(0,0,10,10,20,20,30,30,40,40,50,50,60,60,70,70),
         age_dist_spain5=age_dist_spain5)%>%
  group_by(age2) %>%
  mutate(w=age_dist_spain5/sum(age_dist_spain5)) %>%
  summarise(sum(`X0 + X5`*w),sum(`X10 + X15`*w),sum(`X20 + X25`*w),sum(`X30 + X35`*w),
            sum(`X40 + X45`*w),sum(`X50 + X55`*w),sum(`X60 + X65`*w),sum(`X70 + X75`*w))  %>%
  select(-age2) %>%
  as.matrix(.)
rownames(contact_spain2)=NULL
colnames(contact_spain2)=NULL
contact_spain2

tmp = contact_spain2[,8]*age_dist_spain10[8]/sum(age_dist_spain10[8:9])
tmp2 = contact_spain2[,8]*age_dist_spain10[9]/sum(age_dist_spain10[8:9])
contact_spain2[,8]=tmp
contact_spain2=cbind(contact_spain2,tmp2)
tmp=c(contact_spain2[,9]*age_dist_spain10[1:8]/age_dist_spain10[9],contact_spain2[8,8])
contact_spain2=rbind(contact_spain2,tmp)
rownames(contact_spain2)=NULL
colnames(contact_spain2)=NULL
contact_spain2
paste(t(contact_spain2),sep="' '", collapse=", ")
contact_spain=c(5.35148784013064, 0.792403698555513, 0.597888237674624, 1.81141261889314, 1.15880446573503, 0.426179823985014, 0.272566579444233, 0.0667965764245607, 0.0524292496334014, 0.962410237985733, 6.96132818205631, 1.05406336733019, 1.13101315162061, 1.62168368448698, 0.541175903065237, 0.131491288694527, 0.0540877655519872, 0.0424539866267057, 0.452019159060351, 1.34005151271159, 5.2217886689935, 3.01397324511203, 2.22431927441854, 1.43787568086573, 0.197174914122326, 0.049091677623549, 0.038532511077935, 1.86647113340635, 1.27397280273288, 2.63210519557764, 5.87106096844647, 3.36191516449852, 1.4530452744515, 0.523291006732789, 0.10038463556598, 0.0787928274047232, 1.10821166176017, 1.96582937687123, 2.05373971751821, 3.3943240797736, 4.28881981759823, 1.48129197804708, 0.305549086413716, 0.108848089031009, 0.0854358701807222, 0.894549948516997, 1.57636342878353, 1.99319025424205, 2.30845452446609, 2.57838327243766, 2.46844694468242, 0.521659898101079, 0.123368775343464, 0.0968332909509898, 0.616712123753189, 0.467315717289364, 0.722196642987352, 1.53068487644559, 0.963926067358351, 0.932445331806381, 1.58924678766151, 0.22742940762708, 0.178511442123438, 0.49296743101773, 0.761168134312711, 0.392685418560048, 0.922634583854662, 1.19274137966985, 0.916395267151485, 0.868894131016553, 0.5953456757922, 0.46729231833502, 0.0825597022962617, 0.0695377130823058, 0.063939746734617, 0.176136931824655, 0.229871638980727, 0.228790005278394, 0.317695213030174, 0.5953456757922, 0.5953456757922)
matrix(contact_spain,nrow=9,byrow=TRUE)







source("data/lombardy/data_management_lombardy.R")
incidence_cases
incidence_deaths
agedistr_cases
agedistr_deaths
pop_t
age_dist


source("data/germany/data_management_bavaria.R")
incidence_cases
incidence_deaths
agedistr_cases
agedistr_deaths
pop_t
age_dist

source("data/germany/data_management_badenw.R")
incidence_cases
incidence_deaths
agedistr_cases
agedistr_deaths
pop_t
age_dist

source("data/austria/data_management_austria.R")
incidence_cases
incidence_deaths
agedistr_cases
agedistr_deaths
pop_t
age_dist

source("data/spain/data_management_spain.R")
incidence_cases
incidence_deaths
agedistr_cases
agedistr_deaths
pop_t
age_dist

source("data/switzerland/data_management_switzerland.R")
incidence_cases
incidence_deaths
agedistr_cases
agedistr_deaths
pop_t
age_dist

