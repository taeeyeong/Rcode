###################################################################################################################
################################# 1) time dependent covariates CoxPH #################################################################
###################################################################################################################

setwd("/Users/taeyeong/Desktop/BIBS/코로나임상/대구데이터/약물반응")
library(aods3)
library(stringr)

dataset <- read.csv('Daegu_Dataset_V3.1.csv', header = T, stringsAsFactors = F)
dim(dataset) #55484 315

dataset[which(str_detect(dataset$Daily.New.KSW.CSS,"경증")),"CSS4stage"] <- 1
dataset[which(str_detect(dataset$Daily.New.KSW.CSS,"중등도")),"CSS4stage"] <- 2
dataset[which(str_detect(dataset$Daily.New.KSW.CSS,"중증")),"CSS4stage"] <- 3
dataset[which(str_detect(dataset$Daily.New.KSW.CSS,"최중증")),"CSS4stage"] <- 4

check <- dataset[which(dataset$모든.기관..입원일 != dataset$모든.기관..입원일),]
sum(is.na(dataset$모든.기관..입원일)) #0
sum(is.na(dataset$모든.기관..퇴원일)) #0 ## 이걸로 구해도 되는지 확인 


dataset$모든.기관..입원일<- as.Date(dataset$모든.기관..입원일, "%Y.%m.%d")
dataset$모든.기관..퇴원일 <- as.Date(dataset$모든.기관..퇴원일, "%Y.%m.%d")
dataset$평가일.Date.of.assessment. <- as.Date(dataset$평가일.Date.of.assessment., "%Y.%m.%d")
length(unique(dataset$질본.번호)) #2263

##사망환자 제외 및 censoring 환자
# 사망 환자 기준
# 퇴원시 상태: 4 = 사망
# 사망일 = 데이터 존재할 경우 사망
# 상태(8.3기준)  = ‘사망;
# 이중 하나라도 해당
# 
# Censoring 기준
# 격리해제 날짜가 없음


treatPat <- dataset[-which(dataset$퇴원.시.상태 ==4 | dataset$상태.8.3기준.=="사망"),]
length(unique(treatPat$질본.번호)) #2094
table(treatPat$상태.8.3기준.)
table(treatPat$퇴원.시.상태)

treatPat2 <- NULL
for(id in unique(treatPat$질본.번호)){
  tmp <- treatPat[which(treatPat$질본.번호 %in% id),]
  tmp <- tmp[order(tmp$평가일.Date.of.assessment.),]
  tmp$time <- tmp$평가일.Date.of.assessment. - tmp[1,"모든.기관..입원일"]
  tmp$event <- ifelse(is.na(tmp[1,"모든.기관..퇴원일"]),0,1)
  treatPat2 <- rbind(treatPat2, tmp)
}
sum(is.na(treatPat2$time)) #0
table(treatPat2$time)
a <- treatPat2[which(treatPat2$time<0),c("모든.기관..입원일","평가일.Date.of.assessment.")]
table(treatPat2$event)
treatPat3 <- treatPat2[-which(treatPat2$time<0),]
length(unique(treatPat3$질본.번호)) #2094


med_list <- c("질본.번호","모든.기관..입원일","평가일.Date.of.assessment.", "CSS4stage", "성별", "나이", "치료.약제.Kaletra"                                                                            
              , "치료.약제.Interferon.beta"                                                                    
              , "치료.약제.Interferon.alpha"                                                                   
              , "치료.약제.Ribavirin"                                                                          
              , "치료.약제.Hydroxychloroquine"                                                                 
              , "치료.약제.Steroid"                                                                            
              , "치료.약제.Remdesivir"                                                                         
              , "치료.약제.Favipiravir"                                                                        
              , "치료.약제.Antibiotic"                                                                         
              , "치료.약제.Other","time","event" )
treatPat4 <- treatPat3[,med_list]

for (i in med_list){
  print(i)
  print(table(treatPat4[,i]))
  print("==========================================")
}

##Kaletra, Hydroxychlroquine 이외의 약물을 처방 받은 환자 제외 (Remdesivir: 전체 샘플 중 2명 처방)
removeSamp <- unique(treatPat4[which(treatPat4$치료.약제.Interferon.beta==2 | treatPat4$치료.약제.Ribavirin==4 | treatPat4$치료.약제.Remdesivir==8 | treatPat4$치료.약제.Other==7),"질본.번호"])
length(removeSamp) #236

treatPat5 <- treatPat4[-which(treatPat4$질본.번호 %in% removeSamp),]
dim(treatPat5) #41082 18
length(unique(treatPat5$질본.번호)) #1858

treatPat5$kal <- ifelse(!is.na(treatPat5$치료.약제.Kaletra==1),1,0)
treatPat5$hyd <- ifelse(!is.na(treatPat5$치료.약제.Hydroxychloroquine==5),1,0)
treatPat5$str <- ifelse(!is.na(treatPat5$치료.약제.Steroid==6),1,0)

table(treatPat5$kal)
table(treatPat5$hyd)
table(treatPat5$str)


tCovData <- NULL
for(id in unique(treatPat5$질본.번호)){
  tmp <- treatPat5[which(treatPat5$질본.번호 %in% id),]
  tmp$time1 <- tmp$time
  tmp$time2 <- tmp$time + 1
  tCovData <- rbind(tCovData, tmp)
}
tCovData <- tCovData[,c("질본.번호","time1","time2","event","kal","hyd","str","성별","나이","CSS4stage")]
length(unique(tCovData$질본.번호)) #1858
tCovData$나이 <- scale(tCovData$나이)
table(tCovData$나이)

write.csv(tCovData, "dataset_COX_timedependent_V3.1.csv", row.names = F)


tCovData <- read.csv("dataset_COX_timedependent_V3.1.csv", header = T, stringsAsFactors = F)

fit <- coxph(Surv(time1, time2, event)~factor(kal)+factor(hyd)+factor(str)+factor(성별)+나이+CSS4stage, data = tCovData)
summary(fit)
zph <- cox.zph(fit)


###########################################################
##########################contrast matrix############################################
##################################################################################
library(multcomp)

K <- matrix(c(1, -1, 0, 0, 0, 0), 1)
t <- glht(fit, linfct = K)
summary(t)

K <- matrix(c(1, 0, -1, 0, 0, 0), 1)
t <- glht(fit, linfct = K)
summary(t)


K <- matrix(c(0, 1, -1, 0, 0, 0), 1)
t <- glht(fit, linfct = K)
summary(t)


###########################################################
################Proportional Hazzard Assumption 확인######
###########################################################
tCovData <- read.csv("dataset_COX_timedependent_V3.1.csv", header = T, stringsAsFactors = F)
dim(tCovData) #41082
length(unique(tCovData$질본.번호)) #1858
str_samp <- tCovData[which(tCovData$str==1),"질본.번호"]
tCovData <- tCovData[-which(tCovData$질본.번호 %in% str_samp),]
dim(tCovData) #38217
length(unique(tCovData$질본.번호)) #1775
tCovData <- tCovData[,-7]
names(tCovData)


kal_samp <- tCovData[which(tCovData$kal==1),"질본.번호"]
kal_samp <- unique(kal_samp)
length(kal_samp) #873


hyd_samp <- tCovData[which(tCovData$hyd==1),"질본.번호"]
hyd_samp <- unique(hyd_samp)
length(hyd_samp) #762

# str_samp <- tCovData[which(tCovData$str==1),"질본.번호"]
# str_samp <- unique(str_samp)
# length(str_samp)


tCovData2 <- NULL
for(id in unique(tCovData$질본.번호)){
  tmp <- tCovData[which(tCovData$질본.번호 %in% id),]
  tmp[1,"time"] <- tmp[nrow(tmp),"time2"]
  tCovData2 <- rbind(tCovData2, tmp[1,])
}
dim(tCovData2) #1775 10
tCovData2$kal<- ifelse(tCovData2$질본.번호 %in% kal_samp, 1,0)
tCovData2$hyd<- ifelse(tCovData2$질본.번호 %in% hyd_samp, 1,0)
# tCovData2$str<- ifelse(tCovData2$질본.번호 %in% str_samp, 1,0)
table(tCovData2$kal)
table(tCovData2$hyd)
table(tCovData2$str)

fit <- coxph(Surv(time,event)~kal+hyd+성별+나이+CSS4stage, data = tCovData2)
summary(fit)
zph <- cox.zph(fit)
zph
png("plot_time_kal.png")
plot(zph[1], lwd = 2)
abline(0,0,col=1, lty=3, lwd =2)
abline(h = fit$coefficients[1], col =3, lwd = 2, lty=2)
legend("bottomright", legend = c('Reference line for null effect', "Average hazard over time","Time-varying hazzard"), lty = c(3,2,1), col = c(1,3,1), lwd = 2)
dev.off()
## cutpoint : 27


png("plot_time_hyd.png")
plot(zph[2], lwd = 2)
abline(0,0,col=1, lty=3, lwd =2)
abline(h = fit$coefficients[1], col =3, lwd = 2, lty=2)
legend("bottomright", legend = c('Reference line for null effect', "Average hazard over time","Time-varying hazzard"), lty = c(3,2,1), col = c(1,3,1), lwd = 2)
dev.off()
## cutpoint : 19

## time-varying coefficients
kal_split <- survSplit(Surv(time, event)~.,data = tCovData2, cut = c(27), episode = "tgroup", id="id")
head(kal_split)

fit_split <- coxph(Surv(tstart, time, event)~kal:strata(tgroup)+hyd+성별+나이+CSS4stage, data = kal_split)
fit_split

hyd_split <- survSplit(Surv(time, event)~.,data = tCovData2, cut = c(19), episode = "tgroup", id="id")
head(hyd_split)

fit_split <- coxph(Surv(tstart, time, event)~kal+hyd:strata(tgroup)+성별+나이+CSS4stage, data = hyd_split)
fit_split

##########################################################
########################## GEE ###########################
###########################################################

setwd("/Users/taeyeong/Desktop/BIBS/코로나임상/대구데이터/약물반응")
library(multgee)

dataset <- read.csv('dataset_COX_timedependent_V3.1.csv', header = T, stringsAsFactors = F)
dim(dataset) #41082 10

dataset <- dataset[,-match(c("time2","event"),colnames(dataset))]



##time이 중복된 sample 제거
dataset2 <- NULL
Samp <- NULL
for(id in unique(dataset$질본.번호)){
  idx<-NULL
  tmp <- dataset[which(dataset$질본.번호 %in% id),]
  if (nrow(tmp)>2){
    for(i in nrow(tmp)-1){
      if (tmp[i,"time1"]==tmp[(i+1),"time1"]){
        idx <- i+1
        Samp <- c(Samp, id)
      }
    }
    ifelse(tmp$질본.번호 %in% Samp, tmp2 <- tmp[-idx,],tmp2 <- tmp)
    dataset2 <- rbind(dataset2, tmp2)
  }
  else{
    dataset2 <- rbind(dataset2,tmp)
  }
}
dim(dataset2) #40417 8

write.csv(dataset2,"datasetGEE.csv", row.names = F)

####################################################################################
##reduced form dataset##############################################################
####################################################################################
dataset <- read.csv("datasetGEE.csv", header = T, stringsAsFactors = F)
dim(dataset) #40417 8



dataset2 <- NULL
for (id in unique(dataset$질본.번호)){
  tmp <- dataset[which(dataset$질본.번호 %in% id),]
  remove_idx<- c()
  if (nrow(tmp)>2){
    for (i in 1:(nrow(tmp)-1)){
      if(tmp[i,"CSS4stage"] == tmp[(i+1),"CSS4stage"] & tmp[i,"str"] == tmp[(i+1),"str"] & tmp[i,"kal"] == tmp[(i+1),"kal"] & tmp[i,"hyd"] == tmp[(i+1),"hyd"]){
        remove_idx <- c(remove_idx, i+1)
      }
    }
  }
  if(length(remove_idx)>0){
    tmp <- tmp[-c(remove_idx),]
  }
  dataset2 <- rbind(dataset2, tmp)
}
dim(dataset2) #6188 8
length(unique(dataset2$질본.번호)) #1858
dataset3 <- dataset2[-which(dataset2$질본.번호=="hhh" & dataset2$CSS4stage==4),]
dim(dataset3) #6187 8

write.csv(dataset3, "dataset_GEE_reduced_V3.1.csv", row.names = F)

dataset <- read.csv("dataset_GEE_reduced_V3.1.csv", header = T, stringsAsFactors = F)

for (id in unique(dataset$질본.번호)){
  tmp <- dataset[which(dataset$질본.번호 %in% id),]
  tmp1 <- which(duplicated(tmp$time1))
  print(tmp1)
  print(id)
}
fit <- ordLORgee(CSS4stage~factor(kal)+factor(hyd)+factor(str)+factor(성별)+나이+time1+factor(kal)*factor(hyd)+factor(kal)*factor(str) +factor(hyd)*factor(str)+factor(kal)*factor(hyd)*factor(str), id = 질본.번호, repeated = time1, link = "logit", data = dataset, LORstr = "independence")
fit <- ordLORgee(CSS4stage~나이+factor(성별)+time1, id = 질본.번호, repeated = time1, link = "logit", data = dataset, LORstr = "independence")
fit <- ordLORgee(CSS4stage~나이+factor(성별), id = 질본.번호, repeated = time1, link = "logit", data = dataset, LORstr = "independence")
summary(fit)
write.csv(summary(fit)$coefficients, "Result_GEE_indp_v3.1.csv")


##GEE for server
setwd("/home2/tyjung/Corona_Daegu")
library(multgee)

check <- NULL
for(id in unique(dataset$질본.번호)){
  tmp <- dataset[which(dataset$질본.번호 %in% id),]
  tmp$maxCSS <- max(tmp$CSS4stage)
  check <- rbind(check, tmp[1,])
}
fit <- lm(maxCSS~나이, data = check)
summary(fit)
fit <- lm(CSS4stage~나이, data = dataset)
summary(fit)



