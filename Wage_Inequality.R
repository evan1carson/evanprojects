#****************************************************
#*set up R
rm(list=ls())
pacman::p_load(nnet,mgcv,quantreg,systemfit,AER, car, gmodels, haven, jtools, 
               pastecs, plm, psych, stargazer,summarytools,tidyverse, rdrobust,
               openxlsx,rio,skedastic,ggmap,sf,devtools,forecast,pracma,ggthemes)
setwd("~/Desktop/data task")
options("scipen"=999, digits=3)
cps <- read_csv("Noto_data_task.csv")
#*assumptions: clghsg_all is already logged


cps <- cps %>% 
  mutate(year1 = seq(1:46))

#********************
#*column (1)
#********************
#subset to 1963-1987
cpssub <- subset(cps, year<=1987 & year>=1963)
reg1 <- lm(clphsg_all~eu_lnclg+year1,data=cpssub)
summ(reg1,digits=3)
#********************
#*column (2)
#********************
reg2 <- lm(clphsg_all~eu_lnclg+year1,data=cps)
summ(reg2,digits=3)

#********************
#*column (3)
#********************
sub3 <-  cps %>%  
  mutate(post=as.numeric(cps$year>1992)) %>% 
  mutate(brk = post*year1) %>% 
  mutate(int = post*eu_lnclg)
reg3 <- lm(clphsg_all~eu_lnclg+post+brk+year1,data=sub3)
summ(reg3,digits=2)


#********************
#*column (4)
#********************
sub4 <- cps %>% 
  mutate(time2= year1^2) %>% 
  mutate(yearsq=time2/100)
reg4 <- lm(clphsg_all~eu_lnclg+year1+yearsq,data=sub4)
summ(reg4,digits=3)


#********************
#*column (5)
#********************
sub5 <- sub4 %>% 
  mutate(time3= year1^3) %>% 
  mutate(year3=time3/1000)
reg5 <- lm(clphsg_all~eu_lnclg+year1+yearsq+year3,data=sub5)
summ(reg5,digits=3)

#regressions table
elasticity <-function(regression) {1/summary(regression)$coefficients[2]}
elastcoefs <- c(elasticity(reg1),elasticity(reg2),elasticity(reg3),
                elasticity(reg4),elasticity(reg5))
elastcoefs <- round(elastcoefs,digits=3)
pval <- function(regression){summary(regression)$coefficients[2,4]}
elastp <- c(pval(reg1),pval(reg2),pval(reg3),
            pval(reg4),pval(reg5))
elastp <- round(elastp,digits=3)

stargazer(reg1,reg2,reg3,reg4,reg5,type="text",
          order=c("eu_lnclg","year1","yearsq","year3","brk"),
          covariate.labels=c("CLG/HS relative supply",
                             "Time",
                             "Time2/100",
                             "Time3/1000",
                             "Time x post-1992"),
          omit="post",align=TRUE,dep.var.labels.include = FALSE,dep.var.caption="",
          title="REGRESSION MODELS FOR THE COLLEGE/HIGH SCHOOL LOG WAGE GAP, 1963–2005",
          keep.stat=c("n","rsq", "aic","bic"), 
          add.lines=list(c("Elasticity Estimates:",elastcoefs), 
                         c("p-Values of Elasticity",elastp)),
          no.space=FALSE, df=FALSE,
          report="vcs",notes="Standard errors in parentheses",
          digit.separator = "",digits=3)

#****************************************
#graph detrended supply/wages (panel A)
#*****************************************
cps$detrendsupply<- detrend(cps$eu_lnclg)
cps$detrendwage <- detrend(cps$clphsg_all)
plot <- gather(cps,
               key="type",
               value="logpoints",
               25:26)

ggplot(plot,mapping=aes(x=year,y=logpoints,group=type,color=type))+
  theme_economist()+
  scale_colour_economist(name="",labels=c("Detrended Relative Supply","Detrended Wage Differential"))+
  geom_line()+
  geom_point()+
  geom_hline(yintercept=0)+
  ylab("Log Points")+
  ggtitle("Detrended Collge/High School Wage Differential and Relative Supply, 1963-2008")

#****************************************
#Katz-Murphey Graph (panel B)
#*****************************************

#use PCE deflator to obtain real minimum wage, use 2008 dollars
pce <- read_csv("PCEPI.csv")
pce$date <- as.Date(pce$DATE)
pce$year <- floor_date(pce$date, "year")
pce <- pce %>%
  group_by(year) %>%
  summarize(mean = mean(PCEPI))
pce08 <- pce$mean[[which(pce$year=="2008-01-01")]]
pce$inflation <- (pce08-pce$mean)/pce$mean
cps$rminwage <- cps$minimum*(1+pce$inflation)
cps <- cps %>% 
  mutate(lnrealwage=log(rminwage))
cpssub <- subset(cps,year<=1987)

#Run Katz-Murphy 1963-1987 Model, plug coefficients intro function
km <- lm(clphsg_all~eu_lnclg+year1+mur2554+lnrealwage,data=cpssub)

summ(km,digits=3)
z=function(supply,year,mur,rwage){(-0.811*supply+0.033*year+0.009*mur-0.064*rwage-0.314 )}

#use full data in km model
k <- z(cps$eu_lnclg,cps$year1,cps$mur2554,cps$lnrealwage)

plot1 <- data.frame(year=cps$year,predicted=k,actual=cps$clphsg_all)
plot1a <- gather(plot1,
                 key="level",
                 value="wagedif",
                 2:3)
ggplot(plot1a,mapping=aes(x=year,y=wagedif,group=level,color=level))+
  theme_economist()+
  scale_colour_economist(name="",labels=c("Observed CLG/HS Gap","Katz-Murphy Predicted Wage Gap:1963-1987 Trend"))+
  geom_line()+
  geom_point()+
  geom_vline(xintercept=1987)+
  geom_vline(xintercept=1992)+
  ylab("Log Wage Gap")+
  ggtitle("Katz-Murphy Prediction Model for the College/High School Wage Gap")

######################
#iterate
######################
#create dummies for each year
df = data.frame()
for (i in 1:46){
  output = as.numeric(cps$year1>i)
  df=rbind(df,output)
}
colnames(df) <- sprintf("dum%s",seq(1:46))
#create interaction term
trndbrk <- df*cps$year1
colnames(trndbrk) <- sprintf("brk%s",seq(1:46))
#cps1 <- cbind(cps,df,trndbrk)

#iterate regressions
my_lms <- lapply(1:46, function(x) lm(cps$clphsg_all ~ df[,x] +trndbrk[,x]+ cps$eu_lnclg+ cps$year1))

rsquared <- vector()
for (i in 1:46) {
  rsquared[[i]] = summary(my_lms[[i]])$r.squared
}

max(rsquared)

