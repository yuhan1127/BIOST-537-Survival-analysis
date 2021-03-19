## ----setup, include=FALSE--------------------------------------------------------------------
data<-read.csv("bmt-2.csv")[,-1]
library(survival)
library(table1)


## ----1, echo=TRUE----------------------------------------------------------------------------
km.data<-survfit(Surv(tdfs/30,deltadfs)~1,data=data)
plot(km.data,main="the estimated distribution of disease-free survival time",xlab = "disease-free survival time (in months)",ylab="Survival Probability",xaxt="n")
axis(side=1, at=seq(0,88,by=12))



## ----2, echo=TRUE----------------------------------------------------------------------------
data$male<-factor(data$male,levels=c(0,1),labels=c("Female","Male"))
data$donormale<-factor(data$donormale,levels=c(0,1),labels=c("Female","Male"))
data$cmv<-factor(data$cmv,levels=c(0,1),labels=c("Negative","Positive"))
data$donorcmv<-factor(data$donorcmv,levels=c(0,1),labels=c("Negative","Positive"))
data$mtx<-factor(data$mtx,levels=c(0,1),labels=c("No","Yes"))
data$fab<-factor(data$fab,levels=c(1,0),labels=c("FAB grade 4 or 5 and AML","Other"))

data$disgroup<-factor(data$disgroup,levels=c(1,2,3),labels = c("ALL","AML low risk","AML high risk"))  
label(data$male)<-"Patient sex"
label(data$age)<-"Patient age (years)"
label(data$donormale)<-"Donor sex"
label(data$donorage)<-"Donor age (years)"
label(data$cmv)<-"Patient CMV status"
label(data$donorcmv)<-"Donor CMV status"
label(data$waittime)<-"Waiting time until transplant (days)"
label(data$mtx)<-"Prophylactic use of methotrexate "


table1(~age+male+cmv+donorage+donormale+donorcmv+waittime+mtx|disgroup,data=data)
table1(~age+male+cmv+donorage+donormale+donorcmv+waittime+mtx|fab,data=data)




## ----3, echo=TRUE----------------------------------------------------------------------------
s.data<-with(data,Surv(tdfs,deltadfs))

#univariate screening
coxph(s.data ~ age, data=data)
coxph(s.data ~ male, data=data)
coxph(s.data ~ donorage, data=data)
coxph(s.data ~ donormale, data=data)
coxph(s.data ~ cmv, data=data)
coxph(s.data ~ donorcmv, data=data)
coxph(s.data ~ waittime, data=data)
coxph(s.data ~ disgroup,data=data)
coxph(s.data ~ fab, data=data)
coxph(s.data ~ mtx, data=data)



## ----4, echo=TRUE----------------------------------------------------------------------------

## include time-varying covariate derived from deltadfs
data.tvc<-tmerge(data1=data,
                   data2=data,
                   id=id,
                   deltadfs=event(tdfs,deltadfs),
                   deltar=event(tdfs,deltar),
                   postagvhd=tdc(ta))
head(data.tvc)
s.data.dfs<-with(data.tvc,Surv(tstart,tstop,deltadfs))
s.data.r<-with(data.tvc,Surv(tstart,tstop,deltar))

#aGVHD occurrence and disease free survival
coxph(s.data.dfs~postagvhd+fab,data=data.tvc)
#aGVHD occurrence and relapse
coxph(s.data.r~postagvhd+fab,data=data.tvc)


## ----5, echo=TRUE-------------------------------------------------------------------------------
data_aGVHD<-subset(data,deltaa==1)
s.data.aGVHD<-with(data_aGVHD,Surv(tdfs,deltadfs))

#univariate screening
coxph(s.data.aGVHD ~ age, data=data_aGVHD)
coxph(s.data.aGVHD ~ male, data=data_aGVHD)
coxph(s.data.aGVHD ~ donorage, data=data_aGVHD)
coxph(s.data.aGVHD ~ donormale, data=data_aGVHD)
coxph(s.data.aGVHD ~ cmv, data=data_aGVHD)
coxph(s.data.aGVHD ~ donorcmv, data=data_aGVHD)
coxph(s.data.aGVHD ~ waittime, data=data_aGVHD)
coxph(s.data.aGVHD ~ disgroup,data=data_aGVHD)
coxph(s.data.aGVHD ~ fab, data=data_aGVHD)
coxph(s.data.aGVHD ~ mtx, data=data_aGVHD)



## ----6, echo=TRUE-------------------------------------------------------------------------------
# stratified proportional hazard model adjusting for confounders
# confounders: cmv, patient age, fab
s.data<-with(data,Surv(ta/30,deltaa))
model5<-coxph(s.data~mtx+strata(age)+strata(cmv)+strata(fab),data=data)
model5
## fix the confounders??
# provide an estimate survival function of time from transplant until onset of aGVHD separately for patient without mtx use, aged 28 (??), negative status and fab=other, and smilar patients with mtx use.
# plot(survfit(model5,newdata=data.frame(mtx="No",age=28,cmv="Negative",fab="Other"),conf.int=FALSE),col=4,xlab="Time from transplant until onset of aGVHD (in months)",ylim=c(0.7,1))   
# lines(survfit(model5,newdata=data.frame(mtx="Yes",age=28,cmv="Negative",fab="Other"),conf.int=FALSE),col=2) 

model5_unadjust<-coxph(s.data~mtx,data=data)
plot(survfit(model5_unadjust,newdata=data.frame(mtx="Yes"),conf.int=FALSE),col=4,lwd=1.5,
     xlab="Time from transplant until onset of aGVHD (in months)",ylim=c(0.7,1))   
lines(survfit(model5_unadjust,newdata=data.frame(mtx="No"),conf.int=FALSE),col=2,lwd=1.5,)
legend("topright",legend=c("Yes","No"),col=c(4,2),lwd=c(1.5,1.5))


## ----7, echo=TRUE-------------------------------------------------------------------------------
# Similar to Q4
## include time-varying covariate derived from deltap
data.tvc<-tmerge(data1=data,
                   data2=data,
                   id=id,
                   deltadfs=event(tdfs,deltadfs),
                   deltar=event(tdfs,deltar),
                   postP=tdc(tp))
head(data.tvc)
s.data.dfs<-with(data.tvc,Surv(tstart,tstop,deltadfs))
s.data.r<-with(data.tvc,Surv(tstart,tstop,deltar))

#aGVHD occurrence and disease free survival
coxph(s.data.dfs~postP+fab,data=data.tvc)
#aGVHD occurrence and relapse
coxph(s.data.r~postP+fab,data=data.tvc)

