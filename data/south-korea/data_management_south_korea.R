# Setup ----
source("setup.R")

# Controls ----
day_start = as.Date("2020-02-15")
day_data = as.Date("2020-02-16")
day_tmax = as.Date("2020-03-05")
day_quarantine = as.Date("2020-02-20")

# Population data ----


# Load surveillance data ----
data_south_korea = read.csv("data/south-korea/patient.csv",stringsAsFactors = FALSE) %>%
  tbl_df() %>%
  mutate(confirmed_date=ymd(confirmed_date),
         released_date=ymd(released_date),
         deceased_date=ymd(deceased_date),
         age=2020-birth_year
  )
  
# Case & death incidence ----
t1 = data_south_korea %>%
  group_by(confirmed_date) %>%
  summarise(cases=n()) %>%
  rename(date=confirmed_date)
t1 = data_south_korea %>%
  group_by(deceased_date) %>%
  summarise(deaths=n())%>%
  rename(date=deceased_date) %>%
  right_join(t1) %>%
  gather("type","value",2:3)
ggplot(t1) +
  geom_col(aes(x=date,y=value),stat=identity) +
  facet_wrap(~type,scales="free_y",nrow=2)
# No case by disease onset yet!

# Age distributions ----
t2 = data_south_korea %>%
  mutate(age_group=cut(age,breaks=c(seq(0,80,by=10),1+max(data_south_korea$age,na.rm=T)),right = FALSE),
         death=ifelse(state=="deceased",1,0)) %>%
  filter(!is.na(age_group)) %>%
  group_by(age_group) %>%
  summarise(cases=n(),deaths=sum(death)) %>%
  gather("type","value",2:3)
ggplot(t2) +
  geom_col(aes(x=age_group,y=value),stat=identity) +
  facet_wrap(~type,scales="free_y",nrow=2)

# Format for Stan ----



t0 = 0
tmax = as.numeric(day_tmax-day_start)
S = tmax+time_lag
ts = 1:S
t_data = as.numeric(day_data-day_start)
D = tmax-t_data+1
tswitch = as.numeric(day_tswitch-day_start)




