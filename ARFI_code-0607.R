setwd('~')
rm(list = ls())

library(tidyverse)
library(tableone)
library(reshape2)
library(survival)
library(survminer)
library(cmprsk)
library(ggalluvial)
library(ggsci)
library(patchwork)
library(ComplexUpset)
library(igraph)
library(ggraph)
library(tidygraph)
library(showtext)
library(forestplot)
library(pheatmap)

# Initial processing####
multimerge <- function(df = list(), ...){
  mergedf <- df[[1]]
  df[[1]] <- NULL
  for (i in df){
    mergedf <- merge(mergedf, i, all = T)
  }
  return(mergedf)
}
permanant <- function(df, ...){
  for (i in 1:nrow(df)){
    for (j in 1:(ncol(df)-1)){
      if (!is.na(df[i, j]) & df[i, j] == 1){
        df[i, j+1] <-  1
      }
    }
  }
  return(df)
} #once for permanant
id_match <- function(df, ...){
  df <- mutate(df, HHIDPN=paste0(HHID, PN))
  df[substring(df$HHIDPN, 1, 5)=='00000',]$HHIDPN <- substring(df[substring(df$HHIDPN, 1, 5)=='00000',]$HHIDPN, 6)
  df[substring(df$HHIDPN, 1, 1)==0,]$HHIDPN <- substring(df[substring(df$HHIDPN, 1, 1)==0,]$HHIDPN, 2)
  return(df)
} #link the 'HHID' and 'PN'

hfat00 <- haven::read_sas("./h00f1d.sas7bdat", NULL)[c('HHIDPN', 'G1361', 'G1369', 'G1339', 'G1442', 'G1326', 'G2698', 'G2716')]
hfat02 <- haven::read_sas("./h02f2c.sas7bdat", NULL)[c('HHIDPN', 'HC095', 'HC103', 'HC079', 'HC148', 'HG005', 'HG011')]
hfat04 <- haven::read_sas("./h04f1c.sas7bdat", NULL)[c('HHIDPN', 'JC095', 'JC103', 'JC079', 'JC148', 'JG005', 'JG011')]
hfat06 <- haven::read_sas("./h06f3a.sas7bdat", NULL)[c('HHIDPN', 'KC095', 'KC103', 'KC079', 'KC148', 'KG005', 'KG011')]
hfat08 <- haven::read_sas("./h08f3a.sas7bdat", NULL)[c('HHIDPN', 'LC095', 'LC103', 'LC079', 'LC148', 'LG005', 'LG011')]
hfat10 <- haven::read_sas("./hd10f5f.sas7bdat", NULL)[c('HHIDPN', 'MC095', 'MC103', 'MC079', 'MC148', 'MG005', 'MG011')]
hfat12 <- haven::read_sas("./h12f3a.sas7bdat", NULL)[c('HHIDPN', 'NC095', 'NC103', 'NC079', 'NC148', 'NG005', 'NG011')]
hfat14 <- haven::read_sas("./h14f2b.sas7bdat", NULL)[c('HHIDPN', 'OC095', 'OC103', 'OC079', 'OC148', 'OG005', 'OG011')]
hfat16 <- haven::read_sas("./h16f2b.sas7bdat", NULL)[c('HHIDPN', 'PC095', 'PC103', 'PC079', 'PC148', 'PG005', 'PG011')]
hfat18 <- haven::read_sas("./h18e1a.sas7bdat", NULL)[c('HHIDPN', 'QC095', 'QC103', 'QC079', 'QC148', 'QG005', 'QG011')]
hfat20_c <- id_match(haven::read_sas("./h20c_r.sas7bdat", NULL))[c('HHIDPN', 'RC095', 'RC103', 'RC079', 'RC148', 'RC141', 'RC139')]
hfat20_g <- id_match(haven::read_sas("./h20g_r.sas7bdat", NULL))[c('HHIDPN', 'RG005', 'RG011')]

age_var <- paste0('R', 5:15, 'AGEY_E')
EI_var <- paste0('R', 5:15, 'EI')
HI_var <- paste0('R', 5:15, 'HI')
CI_var <- paste0('R', 5:15, 'CI')
DEP_var <- paste0('R', 5:15, 'DEP')
FRA_var <- paste0('R', 5:15, 'FRAIL')
SLE_var <- paste0('R', 5:15, 'SLEEPR')

# BMI caculation
hfat20 <- merge(hfat20_c, hfat20_g, by = 'HHIDPN', all.x = T)
r15bmi <- hfat20[c('RC141', 'RC139')]
r15bmi[r15bmi == 998|r15bmi==999|r15bmi==98|r15bmi==-8] <- NA
hfat20[c('RC141', 'RC139')] <- r15bmi
hfat20$R15BMI <- hfat20$RC139*0.4539 / (0.3048*hfat20$RC141)^2
hfat20$R15BMI <- ifelse(hfat20$R15BMI > 60, NA, hfat20$R15BMI)

df_locfatraw <- multimerge(list(hfat00, hfat02, hfat04, hfat06, hfat08, hfat10, hfat12, hfat14, hfat16, hfat18, hfat20))


## Visual impairment & hearing impairment
eye <- df_locfatraw[c('G1361', 'HC095', 'JC095', 'KC095', 'LC095', 'MC095', 'NC095', 'OC095', 'PC095', 'QC095', 'RC095')]
hear <- df_locfatraw[c('G1369', 'HC103', 'JC103', 'KC103', 'LC103', 'MC103', 'NC103', 'OC103', 'PC103', 'QC103', 'RC103')]
fall <- df_locfatraw[c('G1339', 'HC079', 'JC079', 'KC079', 'LC079', 'MC079', 'NC079', 'OC079', 'PC079', 'QC079', 'RC079')]
fati <- df_locfatraw[c('G1442', 'HC148', 'JC148', 'KC148', 'LC148', 'MC148', 'NC148', 'OC148', 'PC148', 'QC148', 'RC148')]
chair <- df_locfatraw[c('G2698', 'HG005', 'JG005', 'KG005', 'LG005', 'MG005', 'NG005', 'OG005', 'PG005', 'QG005', 'RG005')]
lift <- df_locfatraw[c('G2716', 'HG011', 'JG011', 'KG011', 'LG011', 'MG011', 'NG011', 'OG011', 'PG011', 'QG011', 'RG011')]

names(eye) <- paste0('R', 5:15, 'EI')
names(hear) <- paste0('R', 5:15, 'HI')
names(fall) <- paste0('R', 5:15, 'FAL')
names(fati) <- paste0('R', 5:15, 'FAT')
names(lift) <- paste0('R', 5:15, 'LIFT')
names(chair) <- paste0('R', 5:15, 'CHAIR')

eye[eye > 5 | eye < 1] <- NA
eye[eye < 4 & eye > 0] <- 0
eye[eye == 4 | eye == 5] <- 1
hear[hear > 5 | hear < 1] <- NA
hear[hear < 4 & hear > 0] <- 0
hear[hear == 4 | hear == 5] <- 1

fall[fall != 1 & fall != 5] <- NA
fall[fall == 5] <- 0
fati[fati != 1 & fati != 5] <- NA
fati[fati == 5] <- 0
chair[chair != 1 & chair != 5] <- NA
chair[chair == 5] <- 0
lift[lift != 1 & lift != 5] <- NA
lift[lift == 5] <- 0

eye <- permanant(eye)
hear <- permanant(hear)

df_locfatraw <- cbind(df_locfatraw, eye, hear, fall, fati)

df_origin <- haven::read_sas("./randhrs1992_2018v1.sas7bdat", NULL)
age <- paste0('R', 5:14, 'AGEY_E')
wave_status <- paste0('R', 5:14, 'IWSTAT')
sleeping_disorders <- paste0('R', 5:14, 'SLEEPR')
depression_cesd <- paste0('R', 5:14, 'CESD')

#Sensitivity analysis
sensi_var <- c('R5STROKE', 'R5CANCRE', 'R5PSYCHE', 'R5ARTHRE', 'R5HIBPE', 'R5DIABE', 'R5LUNGE', 'R5HEARTE')

#Covariates
bmi <- paste0('R', 5:14, 'BMI')
drink <- paste0('R', 5:14, 'DRINK')
drinkd <- paste0('R', 5:14, 'DRINKD')
income <- paste0('R', 5:14, 'ICAP')
smoken <- paste0('R', 5:14, 'SMOKEN')
smokev <- paste0('R', 5:14, 'SMOKEV')
exe <- paste0('R', 5:14, 'VIGACT')
INW <- paste0('INW', 5:14)

df_locrandraw <- df_origin[df_origin$HACOHORT <= 4 & df_origin$R4IWSTAT == 1, ]
df_locrandraw <- data.frame(df_locrandraw[c('HHIDPN', INW, 'RAHRSAMP', 'RAAHDSMP', 'RAGENDER', 'RARACEM', 'RAEDUC', 'HACOHORT', 'RAWTSAMP',
                                            age, wave_status, sleeping_disorders, depression_cesd, 
                                            bmi, drink, drinkd, income, smoken, smokev, 'RADYEAR', exe, sensi_var)])
# define the smoking status and drinking status
df_locrandraw <- within(df_locrandraw, {
  smoke <- NA
  drink <- NA
  smoke[R5SMOKEV==0] <- 0
  smoke[R5SMOKEV==1] <- 1
  smoke[R5SMOKEN==1] <- 2
  drink[R5DRINK==0] <- 0
  drink[R5DRINK==1] <- 1
  drink[R5DRINKD>0] <- 2
})

df_localldata <- merge(df_locrandraw, df_locfatraw, by = 'HHIDPN', all.x = T)

#Define Frailty Index
bmi <- paste0('R', 4:15, 'BMI')
WAS <- df_localldata[bmi]
for (i in 1:nrow(WAS)){
  for (j in 2:ncol(WAS)) {
    if (!is.na(WAS[i, j-1]) & !is.na(WAS[i, j]) & (WAS[i, j] - WAS[i, j-1])/WAS[i, j-1] >= 0.1){
      WAS[i, j-1] <- 1
    }
    else if (!is.na(WAS[i, j-1]) & !is.na(WAS[i, j]) & (WAS[i, j] - WAS[i, j-1])/WAS[i, j-1] < 0.1){
      WAS[i, j-1] <- 0
    }
    else {WAS[i, j-1] <- NA}
  }
}
WAS <- WAS[, -ncol(WAS)]
WAS_var <- paste0('R', 5:15, 'WAS')
names(WAS) <- WAS_var
df_localldata <- cbind(df_localldata, WAS)

was <- df_localldata[c('HHIDPN', WAS_var)]
other <- cbind(df_locfatraw$HHIDPN,fall,fati,chair,lift)
names(other) <- c('HHIDPN', names(other[,-1]))
FRA <- merge(was, other, by = 'HHIDPN', all.x = T)

for (n in 1:11){
  FRA[paste0('R', (n + 4), 'FRAIL')] <- rowSums(FRA[,seq(n+1, n+45, 11)], na.rm = T)
  FRA[paste0('test', n)] <- ifelse((is.na(FRA[, n+1])+is.na(FRA[, n+12])+is.na(FRA[, n+23])+is.na(FRA[, n+34])+is.na(FRA[, n+45]))>2, NA, 0)
}
for (i in seq(57,77,2)){
  FRA[, i] <- FRA[, i]+FRA[, i+1]
}

FRA_var <- paste0('R', 5:15, 'FRAIL')
FRA <- FRA[FRA_var]
FRA[FRA < 3] <- 0
FRA[FRA >= 3] <- 1

df_localldata <- cbind(df_localldata, FRA)


#Define depression by CESD
hfat20_d <- id_match(haven::read_sas("./h20d_r.sas7bdat", NULL))[c('HHIDPN', paste0('RD', 142:146), 'RD174', 'RD184', 'RD124', 'RD129', paste0('RD', 110:117))]

cesd_origin <- hfat20_d[paste0('RD',110:117)]
cesd_origin[cesd_origin != 1 & cesd_origin != 5] <- NA
cesd_origin[cesd_origin == 5] <- 0
cesd_origin$RD113 <- 1-cesd_origin$RD113
cesd_origin$RD115 <- 1-cesd_origin$RD115
cesd_origin$test <- ifelse((is.na(cesd_origin$RD110) + is.na(cesd_origin$RD111) + is.na(cesd_origin$RD113) + is.na(cesd_origin$RD114) + is.na(cesd_origin$RD115) + is.na(cesd_origin$RD116) + is.na(cesd_origin$RD117))>3, NA, 0)

cesd15 <- data.frame(HHIDPN = hfat20_d$HHIDPN, R15CESD = rowSums(cesd_origin[,-which(names(cesd_origin)=='RD112')], na.rm = T), R15SLEEPR = cesd_origin$RD112)
cesd15$R15CESD <- cesd15$R15CESD+cesd_origin$test
cesd15$R15DEP <- ifelse(cesd15$R15CESD < 4, 0, 1)

CESD <- paste0('R', 5:14, 'CESD')
sleep_cesd <-  paste0('R', 5:14, 'SLEEPR')
DEP <- df_localldata[CESD] - df_localldata[sleep_cesd]
DEP[DEP < 4] <- 0
DEP[DEP >= 4] <- 1
names(DEP) <- paste0('R', 5:14, 'DEP')
df_localldata <- cbind(df_localldata, DEP)

df_localldata <- merge(df_localldata, cesd15[, -2],  by = 'HHIDPN', all.x = T)


#Define CI
df_cog <- haven::read_sas("./cogfinalimp_9518wide.sas7bdat", NULL)
df_cog$HHIDPN <- paste0(df_cog$HHID, df_cog$PN)
cog1 <- df_cog[c('HHIDPN', paste0('cogtot27_imp', seq(1998, 2018, 2)))]
cog1[substring(cog1$HHIDPN, 1, 5)=='00000',]$HHIDPN <- substring(cog1[substring(cog1$HHIDPN, 1, 5)=='00000',]$HHIDPN, 6)
cog1[substring(cog1$HHIDPN, 1, 1)==0,]$HHIDPN <- substring(cog1[substring(cog1$HHIDPN, 1, 1)==0,]$HHIDPN, 2)
names(cog1) <- c('HHIDPN', paste0('R', 4:14, 'CI'))

r15delay <- ifelse(hfat20_d$RD124 == 1, 2, 
                   ifelse(hfat20_d$RD124 == 5, 0, 
                          ifelse(hfat20_d$RD124 == 6 & hfat20_d$RD129 == 1, 1, 
                                 ifelse(hfat20_d$RD124 == 6 & hfat20_d$RD129 == 5, 0, NA))))
r15minus <- data.frame(
  r15minus1 = ifelse(hfat20_d$RD142 == 93, 1, 
                     ifelse(hfat20_d$RD142 == -8 | hfat20_d$RD142 == 998 | hfat20_d$RD142 == 999, NA, 0)),
  r15minus2 = ifelse(hfat20_d$RD143 == 86, 1, 
                     ifelse(hfat20_d$RD143 == -8 | hfat20_d$RD143 == 998 | hfat20_d$RD143 == 999, NA, 0)),
  r15minus3 = ifelse(hfat20_d$RD144 == 79, 1, 
                     ifelse(hfat20_d$RD144 == -8 | hfat20_d$RD144 == 998 | hfat20_d$RD144 == 999, NA, 0)),
  r15minus4 = ifelse(hfat20_d$RD145 == 72, 1, 
                     ifelse(hfat20_d$RD145 == -8 | hfat20_d$RD145 == 998 | hfat20_d$RD145 == 999, NA, 0)),
  r15minus5 = ifelse(hfat20_d$RD146 == 65, 1, 
                     ifelse(hfat20_d$RD146 == -8 | hfat20_d$RD146 == 998 | hfat20_d$RD146 == 999, NA, 0)))
r15minus$r15minus_total <- rowSums(r15minus, na.rm = T)
r15minus$test <- ifelse((is.na(r15minus$r15minus1)+is.na(r15minus$r15minus2)+is.na(r15minus$r15minus3)+is.na(r15minus$r15minus4)+is.na(r15minus$r15minus5))>2, NA, 0)
r15minus$r15minus_total <- r15minus$r15minus_total + r15minus$test
cog20 <- data.frame(HHIDPN = hfat20_d$HHIDPN, R15CI = r15delay + r15minus$r15minus_total + hfat20_d$RD174 + hfat20_d$RD184)

cog_all <- merge(cog1, cog20, by = 'HHIDPN', all.x = T)

cog_all[, 2:13][cog_all[, 2:13] <= 11] <- 1
cog_all[, 2:13][cog_all[, 2:13] > 11] <- 0
cog_all[, 3:13] <- permanant(cog_all[, 3:13])

df_localldata <- merge(df_localldata, cog_all[, c(1,3:13)], by = 'HHIDPN')

## mark the final status of ARFI and first time ocurrred the ARFI
timede <- function(var_list){
  tmp <- df_localldata[var_list]
  tmp$timevar <- NA
  tmp$outvar <- NA
  for (i in 1:nrow(tmp)){
    for (j in 1:(ncol(tmp))){
      if (!is.na(tmp[i,j])){
        if (tmp[i,j] == 1){
          tmp$timevar[i] <- j+4
          tmp$outvar[i] <- 1
          break
        }
        else {
          tmp$timevar[i] <- -1
          tmp$outvar[i] <- 0
        }
      }
    }
  }
  df_localldata <<- cbind(df_localldata, tmp[, (ncol(tmp)-1):ncol(tmp)])
}

timede(c(paste0('R', 5:15, 'DEP')))
timede(c(paste0('R', 5:15, 'FRAIL')))
timede(c(paste0('R', 5:15, 'CI')))
timede(c(paste0('R', 5:15, 'EI')))
timede(c(paste0('R', 5:15, 'HI')))
timede(c(paste0('R', 5:15, 'SLEEPR')))
names(df_localldata) <- c(names(df_localldata[, 1:(ncol(df_localldata)-12)]), 'DEPTIME', 'DEPOC', 'FRAILTIME', 'FRAILOC', 'CITIME', 'CIOC', 'EITIME', 'EIOC', 'HITIME', 'HIOC', 'SLEEPTIME', 'SLEEPOC')

write.csv(df_localldata, file = 'loc_alldata.csv', row.names = F)

# Table one####
df_localldata <- filter(read.csv('./loc_alldata.csv', header = T), R5AGEY_E>50, R5AGEY_E<91)
df_localldata <- within(df_localldata, {
  cohort <- NA
  cohort[HACOHORT == 3] <- 'HRS'
  cohort[HACOHORT == 0|HACOHORT == 1] <- 'AHEAD'
  cohort[HACOHORT == 2] <- 'CODA'
  cohort[HACOHORT == 4] <- 'WB'
})

var <- c('R5AGEY_E', 'RAGENDER', 'RARACEM', 'RAEDUC', 'R5BMI', 'drink', 'H5ICAP', 'smoke', 'R5VIGACT', 
         'R5EI', 'R5HI', 'R5CI', 'R5FRAIL', 'R5SLEEPR', 'R5DEP')
catvars <- c('RAGENDER', 'RARACEM', 'RAEDUC', 'drink' ,'smoke', 'R5VIGACT', 
             'R5EI', 'R5HI', 'R5CI', 'R5FRAIL', 'R5DEP', 'R5SLEEPR')

table1 <- CreateTableOne(vars = var,  strata = 'cohort', addOverall = T, factorVars = catvars, data = df_localldata)
tableout <- print(table1, showAllLevels = T, formatOptions = list(big.mark = ','), quote = F, noSpaces = T, printToggle = F)
write.csv(tableout, file = 'table1.csv')


# Step 1: Calculate the incidence rate of different ages####
df_localldata <- filter(read.csv('./loc_alldata.csv', header = T), R5AGEY_E>50, R5AGEY_E<91)

ir <- function(wave1, wave2, wave3, pf = F){
  temp <- df_localldata[!is.na(df_localldata[wave1]) & !is.na(df_localldata$R5AGEY_E) & df_localldata[wave1] == 0, ]
  temp$age <- temp$R5AGEY_E
  temp$outcome <- NA
  temp$year <- NA
  for (i in 1:nrow(temp)){
    if (!is.na(temp[i, wave2]) & temp[i, wave2] == 1 & pf == F){
      temp[i, 'outcome'] <- 1
      temp[i, 'year'] <- 2
    } else if(!is.na(temp[i, wave3])){
      temp[i, 'outcome'] <- temp[i, wave3]
      temp[i, 'year'] <- 4
    }
  }
  temp <- filter(temp, !is.na(temp$outcome) & !is.na(temp$age)) %>%
    select(age, outcome, RAGENDER, year)
  if (pf ==T){
    temp <- filter(temp, age > 65)
  }
  names(temp) <- c('X1', 'X2', 'X3', 'X4')
  return(temp)
}
eye_ie <- ir('R5EI', 'R6EI', 'R7EI')
hear_ie <- ir('R5HI', 'R6HI', 'R7HI')
CI_ie <- ir('R5CI', 'R6CI', 'R7CI')
sleep_ie <- ir('R5SLEEPR','R6SLEEPR','R7SLEEPR')
fra_ie <- ir('R5FRAIL', 'R6FRAIL', 'R7FRAIL', pf = T)
dep_ie <- ir('R5DEP', 'R6DEP', 'R7DEP')

# calculate the age- and sex- incidence rate of ARFI
age_group <- c('51-55', '56-60', '61-65', '66-70', '71-75', '76-80', '81-85', '86-90')
groupir <- function(ie){
  ir <- c(sum(1000 * ie[ie$X1 <= 55 & ie$X1 > 50,]$X2) / sum(ie[ie$X1 <= 55 & ie$X1 > 50,]$X4),
          sum(1000 * ie[ie$X1 <= 60 & ie$X1 > 55,]$X2) / sum(ie[ie$X1 <= 60 & ie$X1 > 55,]$X4),
          sum(1000 * ie[ie$X1 <= 65 & ie$X1 > 60,]$X2) / sum(ie[ie$X1 <= 65 & ie$X1 > 60,]$X4),
          sum(1000 * ie[ie$X1 <= 70 & ie$X1 > 65,]$X2) / sum(ie[ie$X1 <= 70 & ie$X1 > 65,]$X4),
          sum(1000 * ie[ie$X1 <= 75 & ie$X1 > 70,]$X2) / sum(ie[ie$X1 <= 75 & ie$X1 > 70,]$X4),
          sum(1000 * ie[ie$X1 <= 80 & ie$X1 > 75,]$X2) / sum(ie[ie$X1 <= 80 & ie$X1 > 75,]$X4),
          sum(1000 * ie[ie$X1 <= 85 & ie$X1 > 80,]$X2) / sum(ie[ie$X1 <= 85 & ie$X1 > 80,]$X4),
          sum(1000 * ie[ie$X1 <= 90 & ie$X1 > 85,]$X2) / sum(ie[ie$X1 <= 90 & ie$X1 > 85,]$X4))
  ir_total <- sum(1000 * ie[ie$X1 <= 90 & ie$X1 > 50,]$X2) / sum(ie[ie$X1 <= 90 & ie$X1 > 50,]$X4)
  print(ir_total)
  print(sum(ie[ie$X1 <= 90 & ie$X1 > 50,]$X2))
  return(ir)
}

df_groupir <- data.frame(age_group, groupir(eye_ie), groupir(hear_ie), 
                         groupir(CI_ie), groupir(fra_ie), groupir(sleep_ie), groupir(dep_ie))
df_groupir$age_group <- factor(df_groupir$age_group, levels = c('51-55', '56-60', '61-65', '66-70', '71-75', '76-80', '81-85', '86-90'), ordered = T)
names(df_groupir) <- c('age_group', 'Visual impairment', 'Hearing impairment', 'Cognitive impairment', 'Physical frailty', 'Restless sleep', 'Depression')
write.csv(df_groupir, file = 'irboth.csv', row.names = F)
df_longir <- reshape2::melt(df_groupir, id.vars = 'age_group', value.name = 'inpariment')

#male
df_groupir_male <- data.frame(age_group, groupir(eye_ie[eye_ie$X3==1,]), groupir(hear_ie[hear_ie$X3==1,]), 
                              groupir(CI_ie[CI_ie$X3==1,]), groupir(fra_ie[fra_ie$X3==1,]), groupir(sleep_ie[sleep_ie$X3==1,]), groupir(dep_ie[dep_ie$X3==1,]))
df_groupir_male$age_group <- factor(df_groupir_male$age_group, levels = c('51-55', '56-60', '61-65', '66-70', '71-75', '76-80', '81-85', '86-90'), ordered = T)
names(df_groupir_male) <- c('age_group', 'Visual impairment', 'Hearing impairment', 'Cognitive impairment', 'Physical frailty', 'Restless sleep', 'Depression')
write.csv(df_groupir_male, file = 'irmale.csv', row.names = F)
df_longir_male <- reshape2::melt(df_groupir_male, id.vars = 'age_group', value.name = 'inpariment')

#female
df_groupir_female <- data.frame(age_group, groupir(eye_ie[eye_ie$X3==2,]), groupir(hear_ie[hear_ie$X3==2,]), 
                                groupir(CI_ie[CI_ie$X3==2,]), groupir(fra_ie[fra_ie$X3==2,]), groupir(sleep_ie[sleep_ie$X3==2,]), groupir(dep_ie[dep_ie$X3==2,]))
df_groupir_female$age_group <- factor(df_groupir_female$age_group, levels = c('51-55', '56-60', '61-65', '66-70', '71-75', '76-80', '81-85', '86-90'), ordered = T)
names(df_groupir_female) <- c('age_group', 'Visual impairment', 'Hearing impairment', 'Cognitive impairment', 'Physical frailty', 'Restless sleep', 'Depression')
write.csv(df_groupir_female, file = 'irfemale.csv', row.names = F)
df_longir_female <- reshape2::melt(df_groupir_female, id.vars = 'age_group', value.name = 'inpariment')

figure1 <- ggplot(data = df_longir, aes(x = age_group, y = inpariment, group = variable, color = variable, shape = variable)) + 
  geom_point() + geom_line(size = 0.8) + theme_classic() + 
  theme(legend.position = 'none') + 
  labs(x = 'Age groups', y = 'Incidence rate (case/1000 person years)', title = 'A. Both sexes') + 
  coord_cartesian(ylim = c(0,200)) + 
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 50, 100, 150, 200)) +
  ggsci::scale_color_lancet()

figure1_male <- ggplot(data = df_longir_male, aes(x = age_group, y = inpariment, group = variable, color = variable, shape = variable)) + 
  geom_point() + geom_line(size = 0.8) + theme_classic() + 
  theme(legend.title=element_blank()) + 
  labs(x = 'Age groups', y = 'Incidence rate (case/1000 person years)', title = 'C. Males') + 
  coord_cartesian(ylim = c(0,210)) + 
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 50, 100, 150, 200)) +
  ggsci::scale_color_lancet()
figure1_female <- ggplot(data = df_longir_female, aes(x = age_group, y = inpariment, group = variable, color = variable, shape = variable)) + 
  geom_point() + geom_line(size = 0.8) + theme_classic() + 
  theme(legend.title=element_blank()) + 
  labs(x = 'Age groups', y = 'Incidence rate (case/1000 person years)', title = 'B. Females') + 
  coord_cartesian(ylim = c(0,210)) + 
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 50, 100, 150, 200)) +
  ggsci::scale_color_lancet()

p1 <- figure1 + (figure1_female / figure1_male) + plot_layout(guides = 'collect', widths = c(2, 1))
ggsave(file = 'figure1.pdf', plot = p1, width = 12, height = 6.4)

# Step 2: ARFI coexistence and progression####
## Upset plotting
df_localldata <- filter(read.csv('loc_alldata.csv', header = T), R5AGEY_E >50 & R5AGEY_E < 91)

for (j in list(EI_var, HI_var, CI_var, DEP_var, FRA_var, SLE_var)){
  for (i in 1:(length(j)-1)){
    df_localldata[is.na(df_localldata[j[i+1]]), j[i+1]] <- df_localldata[is.na(df_localldata[j[i+1]]), j[i]]
  }
}

df_matrix <- df_localldata %>%
  select(HACOHORT, R5AGEY_E, R5EI, R5HI, R5CI, R5DEP, R5FRAIL, R5SLEEPR) %>%
  mutate(R5FRAIL = ifelse(R5AGEY_E < 65, 0, R5FRAIL)) %>%
  filter(!is.na(R5EI), !is.na(R5HI), !is.na(R5CI), !is.na(R5DEP), !is.na(R5FRAIL), !is.na(R5SLEEPR))
names(df_matrix) <- c('cohort', 'age', 'Visual impairment', 'Hearing impairment', 'Cognitive impairment', 'Depression', 'Physical frailty', 'Restless sleep')

dfall <- df_matrix %>%
  mutate(`Physical frailty` = ifelse(cohort==4, 0, `Physical frailty`)) %>%
  select(`Visual impairment`:`Restless sleep`)
dfwb <- filter(df_matrix, cohort == 4)[,-c(1:2,7)]
dfhrs <- filter(df_matrix, cohort == 3)[,-(1:2)]
dfcoda <- filter(df_matrix, cohort == 2)[,-(1:2)]
dfahead <- filter(df_matrix, cohort == 0|cohort == 1)[,-(1:2)]

percal <- function(df){
  n <- ncol(df)
  for (i in 1:n){
    for (j in i:n){
      tmp <- df[df[,i]==1&df[,j]==1,]
      print(paste(colnames(df)[i],colnames(df)[j]))
      print(paste0(nrow(tmp), ' (', round(nrow(tmp)/nrow(df)*100,1), ')'))
    }
  }
}
percal(dfall)
percal(dfwb)
percal(dfhrs)
percal(dfcoda)
percal(dfahead)

pall <- ComplexUpset::upset(dfall, colnames(dfall), name='Co-occurrence of ARFIs', width_ratio=0.1, min_size=11, 
                            themes=upset_modify_themes(
                              list('intersections_matrix'=theme(text=element_text(size=25)),
                                   "Intersection size"=theme(text=element_text(size=25)),
                                   'overall_sizes'=theme(text=element_text(size=25), axis.text.x=element_text(angle=90)))))


p1 <- ComplexUpset::upset(dfwb, colnames(dfwb), name='Co-occurrence of ARFIs', width_ratio=0.25, min_size=18, 
                          
                          base_annotations=list(
                            'Intersection size'=intersection_size(
                              text_mapping=aes(label=paste0(
                                !!get_size_mode('exclusive_intersection'), 
                                '\n',
                                '(',
                                round(!!get_size_mode('exclusive_intersection')/nrow(dfwb) * 100, digits = 1),
                                '%)'
                              )))
                          ),
                          
                          themes=upset_modify_themes(
                            list('intersections_matrix'=theme(text=element_text(size=30)),
                                 "Intersection size"=theme(text=element_text(size=30)),
                                 'overall_sizes'=theme(text=element_text(size=30), axis.text.x=element_text(angle=90)))))

p2 <- ComplexUpset::upset(dfhrs, colnames(dfhrs), name='Co-occurrence of ARFIs', width_ratio=0.25, min_size=83, 
                          
                          base_annotations=list(
                            'Intersection size'=intersection_size(
                              text_mapping=aes(label=paste0(
                                !!get_size_mode('exclusive_intersection'), 
                                '\n',
                                '(',
                                round(!!get_size_mode('exclusive_intersection')/nrow(dfhrs) * 100, digits = 1),
                                '%)'
                              )))
                          ),
                          
                          themes=upset_modify_themes(
                            list('intersections_matrix'=theme(text=element_text(size=30)),
                                 "Intersection size"=theme(text=element_text(size=30)),
                                 'overall_sizes'=theme(text=element_text(size=30), axis.text.x=element_text(angle=90)))))
p3 <- ComplexUpset::upset(dfcoda, colnames(dfcoda), name='Co-occurrence of ARFIs', width_ratio=0.25, min_size=18,
                          
                          base_annotations=list(
                            'Intersection size'=intersection_size(
                              text_mapping=aes(label=paste0(
                                !!get_size_mode('exclusive_intersection'), 
                                '\n',
                                '(',
                                round(!!get_size_mode('exclusive_intersection')/nrow(dfcoda) * 100, digits = 1),
                                '%)'
                              )))
                          ),
                          
                          themes=upset_modify_themes(
                            list('intersections_matrix'=theme(text=element_text(size=30)),
                                 "Intersection size"=theme(text=element_text(size=30)),
                                 'overall_sizes'=theme(text=element_text(size=30), axis.text.x=element_text(angle=90)))))
p4 <- ComplexUpset::upset(dfahead, colnames(dfahead), name='Co-occurrence of ARFIs', width_ratio=0.25, min_size=37,
                          
                          base_annotations=list(
                            'Intersection size'=intersection_size(
                              text_mapping=aes(label=paste0(
                                !!get_size_mode('exclusive_intersection'), 
                                '\n',
                                '(',
                                round(!!get_size_mode('exclusive_intersection')/nrow(dfahead) * 100, digits = 1),
                                '%)'
                              )))
                          ),
                          
                          themes=upset_modify_themes(
                            list('intersections_matrix'=theme(text=element_text(size=30)),
                                 "Intersection size"=theme(text=element_text(size=30)),
                                 'overall_sizes'=theme(text=element_text(size=30), axis.text.x=element_text(angle=90)))))

ggsave(file = 'upset_p1.pdf', plot = p1, width = 16, height = 9.6)
ggsave(file = 'upset_p2.pdf', plot = p2, width = 16, height = 9.6)
ggsave(file = 'upset_p3.pdf', plot = p3, width = 16, height = 9.6)
ggsave(file = 'upset_p4.pdf', plot = p4, width = 16, height = 9.6)

ggsave(file = 'upset_all.pdf', plot = pall, width = 20, height = 12)

## Sankey diagram plotting 
df_localldata <- filter(read.csv('D:/Dataset/HRS/loc_alldata.csv', header = T), R5AGEY_E>50 & R5AGEY_E<91)

hfat20_pr <- id_match(haven::read_sas("./h20pr_r.sas7bdat", NULL))[c('HHIDPN', 'RCSR', 'HHID')]
hfat20_a <- haven::read_sas("./h20a_h.sas7bdat", NULL)[c('HHID', 'RA044')]
hfat_base <-  merge(hfat20_pr, hfat20_a[!duplicated(hfat20_a$HHID), ], all.x = T, by = 'HHID')[, -1]
names(hfat_base) <- c('HHIDPN', 'INW15', 'R15AGEY_E')
hfat_base$INW15 <- ifelse(hfat_base$INW15==5 | hfat_base$INW15==3, 0, hfat_base$INW15)
hfat_base$HHIDPN <- as.integer(hfat_base$HHIDPN)

df_localldata <- merge(df_localldata, hfat_base, all.x = T, by = 'HHIDPN')
df_localldata$R15AGEY_E <- ifelse(is.na(df_localldata$R15AGEY_E), df_localldata$R13AGEY_E+4, df_localldata$R15AGEY_E)
df_localldata$INW15 <- ifelse(is.na(df_localldata$INW15), 0, df_localldata$INW15)

sk <- filter(df_localldata, R5AGEY_E>50 & R5AGEY_E<91)
sk <- within(sk, {
  R5FRAIL[R5AGEY_E < 65] <- 0
  R7FRAIL[R7AGEY_E < 65] <- 0
  R9FRAIL[R9AGEY_E < 65] <- 0
  R11FRAIL[R11AGEY_E < 65] <- 0
  R13FRAIL[R13AGEY_E < 65] <- 0
})
sk1 <- filter(sk, HACOHORT == 3)
sk2 <- filter(sk, HACOHORT == 0|HACOHORT == 1)
sk3 <- filter(sk, HACOHORT == 2)
sk4 <- filter(sk, HACOHORT == 4)

sk_plot <- function(sk_, skname){
  sk <- sk_
  for (j in list(EI_var, HI_var, CI_var, DEP_var, FRA_var, SLE_var)){
    for (i in 1:(length(j)-1)){
      sk[is.na(sk[j[i+1]]), j[i+1]] <- sk[is.na(sk[j[i+1]]), j[i]]
    }
  }
  
  sk <- mutate(sk, wave5 = R5EI + R5HI + R5CI + R5DEP + R5FRAIL + R5SLEEPR, 
               wave7 = R7EI + R7HI + R7CI + R7DEP + R7FRAIL + R7SLEEPR, 
               wave9 = R9EI + R9HI + R9CI + R9DEP + R9FRAIL + R9SLEEPR, 
               wave11 = R11EI + R11HI + R11CI + R11DEP + R11FRAIL + R11SLEEPR, 
               wave13 = R13EI + R13HI + R13CI + R13DEP + R13FRAIL + R13SLEEPR, 
               wave15 = R15EI + R15HI + R15CI + R15DEP + R15FRAIL + R15SLEEPR) %>%
    filter(!is.na(wave5))
  sk <- within(sk, {
    wave7[RADYEAR < 2004] <- 'Death'
    wave9[RADYEAR < 2008 & !is.na(wave7)] <- 'Death'
    wave11[RADYEAR < 2012 & !is.na(wave9)] <- 'Death'
    wave13[RADYEAR < 2016 & !is.na(wave11)] <- 'Death'
    wave15[RADYEAR < 2018 & !is.na(wave13)] <- 'Death'
  })
  print(nrow(sk))
  sk_long <- melt(sk, id.vars = 'HHIDPN', measure.vars = c('wave5', 'wave7', 'wave9', 'wave11', 'wave13', 'wave15'), value.name = 'num', variable.name = 'wave')
  sta <- melt(sk, id.vars = 'HHIDPN', measure.vars = c('INW5','INW7','INW9','INW11','INW13', 'INW15'), value.name = 'status', variable.name = 'wave')
  sk_age <- melt(sk, id.vars = 'HHIDPN', measure.vars = c('R5AGEY_E', 'R7AGEY_E', 'R9AGEY_E', 'R11AGEY_E', 'R13AGEY_E', 'R15AGEY_E'), value.name = 'age', variable.name = 'wave')
  sk_long <- cbind(sk_long, sta[, 3], sk_age[,3])
  
  sk_long[sk_long[,4] == 0 & sk_long[,3] != 'Death', 3] <- 'Not visited'
  sk_long <- within(sk_long, {
    Status <- num
    Time <- NA
    Status[num == 0] <- 'No functional impairments existing'
    Status[num == 1] <- 'One functional impairment existing'
    Status[num == 2] <- 'Two functional impairments coexisting'
    Status[num == 3] <- 'Three functional impairments coexisting'
    Status[num == 4] <- 'Four functional impairments coexisting'
    Status[num == 5] <- 'Five functional impairments coexisting'
    Status[num == 6] <- 'Six functional impairments coexisting'
    Time[wave == 'wave5'] <- paste0('2000 (', sprintf('%0.1f', mean(sk$R5AGEY_E, na.rm = T)), ' y)')
    Time[wave == 'wave7'] <- paste0('2004 (', sprintf('%0.1f', mean(sk$R7AGEY_E, na.rm = T)), ' y)')
    Time[wave == 'wave9'] <- paste0('2008 (', sprintf('%0.1f', mean(sk$R9AGEY_E, na.rm = T)), ' y)')
    Time[wave == 'wave11'] <- paste0('2012 (', sprintf('%0.1f', mean(sk$R11AGEY_E, na.rm = T)), ' y)')
    Time[wave == 'wave13'] <- paste0('2016 (', sprintf('%0.1f', mean(sk$R13AGEY_E, na.rm = T)), ' y)')
    Time[wave == 'wave15'] <- paste0('2020 (', sprintf('%0.1f', mean(sk$R15AGEY_E, na.rm = T)), ' y)')
  })
  sk_long$Status <- factor(sk_long$Status, levels = c('No functional impairments existing', 'One functional impairment existing', 'Two functional impairments coexisting', 'Three functional impairments coexisting', 'Four functional impairments coexisting', 'Five functional impairments coexisting', 'Six functional impairments coexisting', 'Not visited', 'Death'))
  sk_long$Time <- factor(sk_long$Time)
  sk_long$id <- rep(1:nrow(sk), 6)
  
  p_sk <- ggplot(sk_long, aes(x = Time, stratum = Status, alluvium = id, fill = Status, label = Status)) +
    geom_flow(color ="darkgray") + 
    geom_stratum() + labs(x = '', y = 'Number of participants', title = skname) + 
    theme(legend.title = element_text(size=10), legend.text=element_text(size=8), legend.key.size = unit(35, "pt")) + 
    guides(shape = guide_legend(override.aes = list(size = 0.2))) + theme_minimal() +
    ggsci::scale_fill_lancet()
  
  sk_long <- mutate(sk_long, Status = ifelse(Status=='Death'|Status=='Not visited', NA, Status))
  
  table <- CreateTableOne(vars = 'Status', strata = 'Time', data = sk_long, factorVars = 'Status')
  tableout <- print(table, showAllLevels = T, formatOptions = list(big.mark = ','), quote = F, noSpaces = T, printToggle = F)
  write.csv(tableout, file = paste0(skname, '.csv'))
  
  return(p_sk)
}

psk1 <- sk_plot(sk1, 'B. HRS original cohort (born in 1931-1941)')
psk2 <- sk_plot(sk2, 'D. AHEAD cohort (born before 1924)')
psk3 <- sk_plot(sk3, 'C. CODA cohort (born in 1924-1930)')
psk4 <- sk_plot(sk4, 'A. WB cohort (born in 1942-1947)')
all <- sk_plot(sk, 'All')

psk <- wrap_plots(psk4, psk1, psk3, psk2) + plot_layout(guides = 'collect')
ggsave(file = 'sk_new.pdf', plot = psk, width = 16, height = 12)

# Step 3: Trajectory network construction####
###Identification of significant associations
df_localldata <- filter(read.csv('./loc_alldata.csv', header = T), R5AGEY_E >50 & R5AGEY_E < 91)
#Sensitivity analysis
#Step 1
#df_localldata <- filter(df_localldata, R5STROKE != 1 & R5CANCRE != 1 & R5PSYCHE != 1 & G1326 != 1)
#Step 2
#df_localldata <- filter(df_localldata, (R5EI == 0 & R5HI == 0 & R5CI == 0 & R5DEP == 0 & R5FRAIL == 0 & R5SLEEPR == 0))
age <- paste0('R', 5:14, 'AGEY_E')
bmi <- paste0('R', 5:14, 'BMI')
drink <- paste0('R', 5:14, 'DRINKD')

income <- paste0('R', 5:14, 'ICAP')

smoke <- paste0('R', 5:14, 'SMOKEN')
exe <- paste0('R', 5:14, 'VIGACT')
  
INW <- paste0('INW', 5:14)

exe_var <- paste0('R', 7:14, 'VGACTX')
df_tra <- df_localldata[c('HHIDPN', INW, age, 'RAGENDER', 'RARACEM', 'RAEDUC', bmi, drink, income, smoke, exe, 'drink', 'smoke',
                          'DEPTIME', 'DEPOC', 'FRAILTIME', 'FRAILOC', 'CITIME', 'CIOC', 'EITIME', 'EIOC', 'HITIME', 'HIOC', 'SLEEPTIME', 'SLEEPOC', 'RADYEAR')]
df_tra[exe_var][df_tra[exe_var] <= 3] <- 0
df_tra[exe_var][df_tra[exe_var] > 3] <- 1

## Calculate t1 for ARFI when treating as exposure
T1 <- function(var_list, t, oc, var_name){
  df_tmp <- df_localldata[c(var_list, t, oc)]
  df_tmp$t1 <- NA
  for (i in 1: nrow(df_tmp)){
    if (!is.na(df_tmp[i,oc]) & df_tmp[i,oc] == 1){
      df_tmp[i, 't1'] <- 2 * df_tmp[i, t] + 1990
    } else if (!is.na(df_tmp[i,oc]) & df_tmp[i,oc] == 0){
      for (j in 1:length(var_list)){
        if (!is.na(df_tmp[i,j])){
          df_tmp[i, 't1'] <- 2 * j + 1998
          break
        }
      }
    }
  }
  df_tra <<- cbind(df_tra, df_tmp$t1)
  names(df_tra) <<- c(names(df_tra[,1:(ncol(df_tra)-1)]), var_name)
}
T1(EI_var, 'EITIME', 'EIOC', 't1_vi')
T1(HI_var, 'HITIME', 'HIOC', 't1_hi')
T1(CI_var, 'CITIME', 'CIOC', 't1_ci')
T1(DEP_var, 'DEPTIME', 'DEPOC', 't1_dep')
T1(FRA_var, 'FRAILTIME', 'FRAILOC', 't1_fra')
T1(SLE_var, 'SLEEPTIME', 'SLEEPOC', 't1_sle')

## Calculate t2 for ARFI when treating as outcome
T2 <- function(var_list, t, oc, var_name){
  df_tmp <- df_localldata[c(var_list, t, oc)]
  df_tmp$t2 <- NA
  for (i in 1: nrow(df_tmp)){
    if (!is.na(df_tmp[i,oc]) & df_tmp[i,oc] == 1){
      df_tmp[i, 't2'] <- 2 * df_tmp[i, t] + 1990
    } else if (!is.na(df_tmp[i,oc]) & df_tmp[i,oc] == 0){
      for (j in length(var_list):1){
        if (!is.na(df_tmp[i,j])){
          df_tmp[i, 't2'] <- 2 * j + 1998
          break
        }
      }
    }
  }
  df_tra <<- cbind(df_tra, df_tmp$t2)
  names(df_tra) <<- c(names(df_tra[,1:(ncol(df_tra)-1)]), var_name)
}
T2(EI_var, 'EITIME', 'EIOC', 't2_vi')
T2(HI_var, 'HITIME', 'HIOC', 't2_hi')
T2(CI_var, 'CITIME', 'CIOC', 't2_ci')
T2(DEP_var, 'DEPTIME', 'DEPOC', 't2_dep')
T2(FRA_var, 'FRAILTIME', 'FRAILOC', 't2_fra')
T2(SLE_var, 'SLEEPTIME', 'SLEEPOC', 't2_sle')

# Define the mortality outcome and time
df_tra <- within(df_tra, {
  DEATHOC <- NA
  DEATHOC[!is.na(RADYEAR)] <- 1
  DEATHOC[is.na(RADYEAR)] <- 0
  t2_dead <- RADYEAR
  t2_dead[is.na(RADYEAR)] <- 2020
})
df_tra$DEATHTIME <- (df_tra$RADYEAR-1990)/2

write.csv(df_tra, file = 'tra.csv', row.names = F)
df_tra <- read.csv('./tra.csv', header = T)

###Testing mediation roles
median_role <- function(expo_, t1_){
  t1_list = c('t1_vi', 't1_hi', 't1_ci', 't1_dep', 't1_fra', 't1_sle')
  t2_list = c('t2_vi', 't2_hi', 't2_ci', 't2_dep', 't2_fra', 't2_sle')
  oc_list = c('EIOC', 'HIOC', 'CIOC', 'DEPOC', 'FRAILOC', 'SLEEPOC')
  
  index <- which(t1_list == t1_)
  
  t1_list <- t1_list[-index]
  t2_list <- t2_list[-index]
  oc_list <- oc_list[-index]
  
  for (i in 1:5){
    tmp <- df_tra[c(expo_, oc_list[i], t1_, t2_list[i], t1_list[-i], oc_list[-i], 'R5AGEY_E', 'RAGENDER', 'R5BMI', 'drink', 'H5ICAP', 'smoke', 'R5VIGACT', 'RARACEM', 'RAEDUC')]
    names(tmp) <- c('expo', 'oc', 't1', 't2', paste0('m', 1:4), paste0('moc', 1:4), names(tmp[, 13:ncol(tmp)]))
    tmp <- filter(tmp, !is.na(expo) & !is.na(oc)) %>%
      filter(!(oc == 1 & t2 <= t1)) %>%
      filter(t2 > t1)
    tmp <- within(tmp, {
      d1 <- NA
      d1[m1<=t1] <- 1
      d1[m1>t1 & m1<=t2] <- 2
      d1[m1>t2] <- 0
      d1[moc1==0] <- 0
      d2 <- NA
      d2[m2<=t1] <- 1
      d2[m2>t1 & m2<=t2] <- 2
      d2[m2>t2] <- 0
      d2[moc2==0] <- 0
      d3 <- NA
      d3[m3<=t1] <- 1
      d3[m3>t1 & m3<=t2] <- 2
      d3[m3>t2] <- 0
      d3[moc3==0] <- 0
      d4 <- NA
      d4[m4<=t1] <- 1
      d4[m4>t1 & m4<=t2] <- 2
      d4[m4>t2] <- 0
      d4[moc4==0] <- 0
    })
    tmp_cox <- coxph(Surv(t1,t2,oc)~expo+R5AGEY_E+RAGENDER+t1+
                       R5BMI+drink+H5ICAP+smoke+R5VIGACT+RARACEM+RAEDUC+d1+d2+d3+d4, data=tmp)
    print(ShowRegTable(tmp_cox, printToggle = F)[1,1:2])
    print(paste0('ori: ', summary(tmp_cox)$coefficients[1,5]))
    print(paste0(nrow(filter(tmp, expo==1 & oc==1)), '/', nrow(filter(tmp, expo==1))))
  }
  tmp <- df_tra[c(expo_, 'DEATHOC', t1_, 't2_dead', t1_list, oc_list, 'R5AGEY_E', 'RAGENDER', 'R5BMI', 'drink', 'H5ICAP', 'smoke', 'R5VIGACT', 'RARACEM', 'RAEDUC')]
  names(tmp) <- c('expo', 'oc', 't1', 't2', paste0('m', 1:5), paste0('moc', 1:5), names(tmp[, 15:ncol(tmp)]))
  tmp <- filter(tmp, !is.na(expo) & !is.na(oc)) %>%
    filter(!(oc == 1 & t2 <= t1)) %>%
    filter(t2 > t1)
  tmp <- within(tmp, {
    d1 <- NA
    d1[m1<=t1] <- 1
    d1[m1>t1 & m1<=t2] <- 2
    d1[m1>t2] <- 0
    d1[moc1==0] <- 0
    d2 <- NA
    d2[m2<=t1] <- 1
    d2[m2>t1 & m2<=t2] <- 2
    d2[m2>t2] <- 0
    d2[moc2==0] <- 0
    d3 <- NA
    d3[m3<=t1] <- 1
    d3[m3>t1 & m3<=t2] <- 2
    d3[m3>t2] <- 0
    d3[moc3==0] <- 0
    d4 <- NA
    d4[m4<=t1] <- 1
    d4[m4>t1 & m4<=t2] <- 2
    d4[m4>t2] <- 0
    d4[moc4==0] <- 0
    d5 <- NA
    d5[m5<=t1] <- 1
    d5[m5>t1 & m5<=t2] <- 2
    d5[m5>t2] <- 0
    d5[moc5==0] <- 0
  })
  motality_cox <- coxph(Surv(t1,t2,oc)~expo+R5AGEY_E+RAGENDER+t1+
                          R5BMI+drink+H5ICAP+smoke+R5VIGACT+RARACEM+RAEDUC+
                          d1+d2+d3+d4+d5, data=tmp)
  print(ShowRegTable(motality_cox, printToggle = F)[1,1:2])
  print(paste0('ori: ', summary(motality_cox)$coefficients[1,5]))
  print(paste0(nrow(filter(tmp, expo==1 & oc==1)), '/', nrow(filter(tmp, expo==1))))
}

median_role(expo_ = 'EIOC', t1_ = 't1_vi')
median_role(expo_ = 'HIOC', t1_ = 't1_hi')
median_role(expo_ = 'CIOC', t1_ = 't1_ci')
median_role(expo_ = 'DEPOC', t1_ = 't1_dep')
median_role(expo_ = 'FRAILOC', t1_ = 't1_fra')
median_role(expo_ = 'SLEEPOC', t1_ = 't1_sle')

###Plotting the trajectory network
network_plot <- function(pathway, filename){
  label <- c('Visual impairment', 'Hearing impairment', 'Cognitive impairment', 'Physical frailty', 'Depression', 'Restless sleep', 'Mortality')
  ARFI <- c('Visual impairment (VI)', 'Hearing impairment (HI)', 'Cognitive impairment (CI)', 'Physical frailty (PF)', 'Depression', 'Restless sleep (RS)', 'Mortality')
  n <- data.frame(label,ARFI)
  hrplot <- read.csv(pathway, header = T)
  color_var1 <- c(rep('#ad002a', 100*nrow(hrplot[hrplot$from == 'Visual impairment',])), 
                  rep('#42b540', 100*nrow(hrplot[hrplot$from == 'Hearing impairment',])), 
                  rep('#00468b', 100*nrow(hrplot[hrplot$from == 'Cognitive impairment',])), 
                  rep('#ed0000', 100*nrow(hrplot[hrplot$from == 'Depression',])), 
                  rep('#925e9f', 100*nrow(hrplot[hrplot$from == 'Physical frailty',])), 
                  rep('#fdaf91', 100*nrow(hrplot[hrplot$from == 'Restless sleep',])))
  
  color_var2 <- c(rep('#ad002a', 100*nrow(hrplot[hrplot$from == 'Visual impairment',])), 
                  rep('#42b540', 100*nrow(hrplot[hrplot$from == 'Hearing impairment',])), 
                  rep('#00468b', 100*nrow(hrplot[hrplot$from == 'Cognitive impairment',])), 
                  rep('#925e9f', 100*nrow(hrplot[hrplot$from == 'Physical frailty',])), 
                  rep('#ed0000', 100*nrow(hrplot[hrplot$from == 'Depression',])), 
                  rep('#fdaf91', 100*nrow(hrplot[hrplot$from == 'Restless sleep',]))) 
  
  hrplot <- graph_from_data_frame(filter(hrplot, HR>1), vertices=n, directed = TRUE)
  d1 <- as_tbl_graph(hrplot)
  
  p1 <- ggraph(d1, layout = "circle")+
    geom_edge_parallel(aes(alpha=logP, width=HR), color = color_var1, arrow = arrow(length = unit(2.5, 'mm')), 
                       start_cap = circle(3, 'mm'), end_cap = circle(4, 'mm')) +
    labs(edge_alpha="-log P", edge_width="Hazard Ratio") +
    geom_node_point(aes(fill = ARFI, color = ARFI), 
                    shape=21, size = 5) + 
    scale_edge_width(range=c(0.5,1.5)) +
    theme_graph() + theme(legend.position = 'none') + ggsci::scale_fill_lancet() + ggsci::scale_color_lancet()
  
  p2 <- ggraph(d1, layout = "circle")+
    geom_edge_link(aes(alpha=logP, width=HR), color = color_var2, arrow = arrow(length = unit(2.5, 'mm')), end_cap = circle(2, 'mm')) +
    labs(edge_alpha="-log P", edge_width="Hazard Ratio") +
    geom_node_point(aes(fill = ARFI, color = ARFI), shape=21, size = 2) + scale_edge_width(range=c(0.5,1.5)) +
    #geom_node_text(aes(label = name)) +
    theme_graph() + facet_edges(~from) + ggsci::scale_fill_lancet() + ggsci::scale_color_lancet()
  
  p_net <- p1 + p2 + plot_layout(guides = 'collect', widths = c(1.25, 2))
  
  font_add('Arial','/Library/Fonts/Arial.ttf')
  showtext_auto() 
  ggsave(file = filename, plot = p_net, width = 16, height = 6.4)
}

network_plot('./hrplot.csv', 'figure_tra.pdf')
network_plot('./hrplot1.csv', 'figure_tra1.pdf')
network_plot('./hrplot2.csv', 'figure_tra2.pdf')
network_plot('./hrplot3.csv', 'figure_tra3.pdf')

#All
label <- c('Visual impairment', 'Hearing impairment', 'Cognitive impairment', 'Physical frailty', 'Depression', 'Restless sleep', 'Mortality')
ARFI <- c('Visual impairment (VI)', 'Hearing impairment (HI)', 'Cognitive impairment (CI)', 'Physical frailty (PF)', 'Depression', 'Restless sleep (RS)', 'Mortality')
n <- data.frame(label, ARFI)
hrplot <- read.csv('./hrplot.csv', header = T)
hrplot <- graph_from_data_frame(filter(hrplot, HR>1), vertices=n, directed = TRUE)
d1 <- as_tbl_graph(hrplot)

p_final <- ggraph(d1, layout = "circle")+
  geom_edge_parallel(aes(alpha=logP, width=HR), 
                     color = c(rep('#ad002a', 600), rep('#42b540', 500), rep('#00468b', 600), rep('#ed0000', 600), rep('#925e9f', 600), rep('#fdaf91', 500)), 
                     arrow = arrow(length = unit(2.5, 'mm')), 
                     start_cap = circle(3, 'mm'), end_cap = circle(4, 'mm')) +
  labs(edge_alpha="-log P", edge_width="Hazard Ratio") +
  geom_node_point(aes(fill = ARFI, color = ARFI), 
                  shape=21, size = 5) + 
  scale_edge_width(range=c(0.5,1.5)) +
  theme_graph() + ggsci::scale_fill_lancet() + ggsci::scale_color_lancet()
font_add('Arial','/Library/Fonts/Arial.ttf')
showtext_auto()
ggsave(file = 'tra_final.pdf', plot = p_final, width = 8, height = 6.4)

###The association of number of ARFI with mortality
df_localldata <- filter(read.csv('./loc_alldata.csv', header = T), R5AGEY_E >50 & R5AGEY_E < 91) %>%
  mutate(R5FRAIL = ifelse(R5AGEY_E<65, 0, R5FRAIL),
         noa = R5EI+R5HI+R5CI+R5DEP+R5SLEEPR+R5FRAIL)
df_tra <- read.csv('./tra.csv', header = T) %>%
  left_join(df_localldata[c('HHIDPN','noa')], by = 'HHIDPN') %>%
  mutate(t1 = 2000, `Each additional ARFI`=noa,
         noa = ifelse(noa==6,5,noa)) %>%
  filter(!is.na(noa) & !is.na(DEATHOC)) %>%
  filter(!(DEATHOC==1 & t2_dead<=t1)) %>% filter(t2_dead>t1) %>%
  mutate(`Number of ARFIs`=factor(noa, labels = c('Without any ARFI', 'With one ARFI', 'With two ARFIs', 'With three ARFIs', 'With four ARFIs', 'With more than four ARFIs')),
         Age=R5AGEY_E, BMI=R5BMI, Income=H5ICAP, 
         `Drinking status`=factor(drink, labels = c('Never', 'Former', 'Current')),
         `Smoking status`=factor(smoke,  labels = c('Never', 'Former', 'Current')),
         `Vigorous physical activity`=factor(R5VIGACT, labels = c('No','Yes')),
         Race=factor(RARACEM, labels = c('White race', 'Black race', 'Others')),
         Education=factor(RAEDUC, labels = c('Lower than High-School','GED','High school','Some college','College and above')),
         Gender=factor(RAGENDER, labels = c('Male', 'Female')))


num_cox <- coxph(Surv(t1,t2_dead,DEATHOC)~`Number of ARFIs`+Age+Gender+Race+BMI+Income+Education+
                   `Drinking status`+`Smoking status`+`Vigorous physical activity`, data=df_tra)
each_cox <- coxph(Surv(t1,t2_dead,DEATHOC)~`Each additional ARFI`+Age+Gender+Race+BMI+Income+Education+
                    `Drinking status`+`Smoking status`+`Vigorous physical activity`, data=df_tra)

forest_results <- data.frame(rbind(c(NA,NA,NA),
                                   c(1,1,1),
                                   cbind(summary(num_cox)[["conf.int"]][1:5, c(1,3,4)]),
                                   cbind(summary(each_cox)[["conf.int"]][, c(1,3,4)]))[1:8,])
colnames(forest_results) <- c('mean', 'lower', 'upper')
rownames(forest_results) <- NULL
fordata <- data.frame(
  df_tra %>%
    mutate(t = t2_dead - t1) %>%
    group_by(`Number of ARFIs`) %>%
    summarize(sum(DEATHOC), sum(t))
)
df_tra %>%
  mutate(t = t2_dead - t1) %>%
  summarize(sum(DEATHOC), sum(t))

tabletext <- cbind(
  c('','Without any ARFI', 'With one ARFI', 'With two ARFIs', 'With three ARFIs', 'With four ARFIs', 'With more than four ARFIs', 'Each additional ARFI'),
  paste0(c('Cases', fordata[,2], 7981),'/',c('Person-years', fordata[,3], 227657)),
  c('HR [95% CI]',"Ref",
    as.vector(ShowRegTable(num_cox, printToggle = F)[1,1]),
    as.vector(ShowRegTable(num_cox, printToggle = F)[2,1]),
    as.vector(ShowRegTable(num_cox, printToggle = F)[3,1]),
    as.vector(ShowRegTable(num_cox, printToggle = F)[4,1]),
    as.vector(ShowRegTable(num_cox, printToggle = F)[5,1]),
    as.vector(ShowRegTable(each_cox, printToggle = F)[1,1])),
  c(paste('P-value'), '', rep('<0.0001',6))
)

pdf('forest.pdf', width = 10, height = 3.5)
forestplot(labeltext = tabletext, forest_results, graph.pos=4, graphwidth = unit(7, "cm"), lineheight = unit(0.8, 'cm'),
           hrzl_lines = list('1' = gpar(lwd=2),
                             '2' = gpar(lwd=1)), new_page = T, align = 'l', colgap = unit(0.01, "npc"), boxsize = 0.2, 
           is.summary=c(T, rep(F,7)), xlog=T, vertices = T, title = 'The association of the number of ARFIs with mortality')
dev.off()

# Exploration analysis####
df_localldata <- filter(read.csv('D:/Dataset/HRS/loc_alldata.csv', header = T), R5AGEY_E >50 & R5AGEY_E < 91)



dfpair <- df_localldata %>% 
  mutate(T1 = 2000, 
         T2 = ifelse(is.na(RADYEAR), 2019, RADYEAR), 
         OC = ifelse(is.na(RADYEAR), 0, 1)) %>% 
  filter(T2>2000) %>% 
  select(R5AGEY_E, RAGENDER, R5BMI, drink, H5ICAP, smoke, R5VIGACT, RARACEM, RAEDUC, 
         R5EI, R5HI, R5CI, R5DEP, R5FRAIL, R5SLEEPR, 
         T1, T2, OC)

ARFI_var <- c('R5EI', 'R5HI', 'R5CI', 'R5DEP', 'R5FRAIL', 'R5SLEEPR')


pair_HR <- matrix(data = NA, nrow = 5, ncol = 10)

for (i in 1:5) {
  
  row_HR <- rep('', 10)
  
  for (j in (i+1):6){
    
    
    tmp <- dfpair %>% 
      mutate(e1 = dfpair[[ARFI_var[i]]], 
             e2 = dfpair[[ARFI_var[j]]], 
             EXPO = e1+e2, 
             OtherARFIs = rowSums(dfpair[ARFI_var[-c(i,j)]]), 
             OtherARFIs = ifelse(OtherARFIs>0, 1, OtherARFIs))
    
    
    pair_cox <- coxph(Surv(T1, T2, OC)~factor(EXPO)+
                        R5AGEY_E+RAGENDER+R5BMI+drink+H5ICAP+smoke+R5VIGACT+RARACEM+RAEDUC+OtherARFIs, data=tmp)
    
    
    row_HR[2*(j-1)-1] <- ShowRegTable(pair_cox, printToggle = F)[2,1]
    row_HR[2*(j-1)] <- summary(pair_cox)$coefficients[2,5]
    
    
  }
  
  pair_HR[i, ] <- row_HR
  
}

rownames(pair_HR) <- c('VI', 'HI', 'CI', 'Depression', 'PF')
colnames(pair_HR) <- c('HI', 'P-value', 'CI', 'P-value', 'Depression', 'P-value', 'PF', 'P-value', 'RS', 'P-value')
pair_HR

write.csv(x = pair_HR, file = 'pair_HR.csv')

hr_m <- read.csv('HR.csv', row.names = 1)


heatplot <- pheatmap(hr_m, cluster_row = F, cluster_cols = F, 
                     border='black', color = colorRampPalette(c("white","red"))(100))

ggsave(filename = 'Hheatplot.pdf', plot = heatplot, height = 3.5, width = 5)



# Sensitive analysis####
df_tra <- read.csv('./tra.csv', header = T)
FG_test <- function(expo_, t1_, t1_list, t2_list, oc_list){
  for (i in 1:5){
    tmp <- df_tra[c(expo_, oc_list[i], t1_, t2_list[i], t1_list[-i], oc_list[-i], 'R5AGEY_E', 'RAGENDER', 'R5BMI', 'drink', 'H5ICAP', 'smoke', 'R5VIGACT', 'RARACEM', 'RAEDUC', 'RADYEAR')]
    names(tmp) <- c('expo', 'oc', 't1', 't2', 'm1', 'm2', 'm3', 'm4', 'moc1', 'moc2', 'moc3', 'moc4', names(tmp[, 13:ncol(tmp)]))
    tmp <- filter(tmp, !is.na(expo) & !is.na(oc)) %>%
      filter(!(oc == 1 & t2 <= t1)) %>%
      filter(t2 > t1)
    tmp <- within(tmp, {
      d1 <- NA
      d1[m1<=t1] <- 1
      d1[m1>t1 & m1<=t2] <- 2
      d1[m1>t2] <- 0
      d1[moc1==0] <- 0
      d2 <- NA
      d2[m2<=t1] <- 1
      d2[m2>t1 & m2<=t2] <- 2
      d2[m2>t2] <- 0
      d2[moc2==0] <- 0
      d3 <- NA
      d3[m3<=t1] <- 1
      d3[m3>t1 & m3<=t2] <- 2
      d3[m3>t2] <- 0
      d3[moc3==0] <- 0
      d4 <- NA
      d4[m4<=t1] <- 1
      d4[m4>t1 & m4<=t2] <- 2
      d4[m4>t2] <- 0
      d4[moc4==0] <- 0
      oc[oc==0 & RADYEAR-t2 <= 2 & t2 != 2020] <- 2
      t2[oc==2 & !is.na(RADYEAR) & !is.na(oc)] <- RADYEAR[oc==2 & !is.na(RADYEAR) & !is.na(oc)]
    })
    cova <- tmp[c('expo','R5AGEY_E','RAGENDER','t1','R5BMI','drink','H5ICAP','smoke','R5VIGACT','RARACEM','RAEDUC','d1','d2','d3','d4')]
    fit_crr <- crr(ftime = (tmp$t2-tmp$t1), fstatus = tmp$oc, cova, failcode = 1, cencode = 0)
    fit_crr[['converged']] <- T
    print(summary(fit_crr))
  }
}
FG_test(expo_ = 'EIOC', t1_ = 't1_vi',
        t1_list = c('t1_hi', 't1_ci', 't1_dep', 't1_fra', 't1_sle'),
        t2_list = c('t2_hi', 't2_ci', 't2_dep', 't2_fra', 't2_sle'),
        oc_list = c('HIOC', 'CIOC', 'DEPOC', 'FRAILOC', 'SLEEPOC'))
FG_test(expo_ = 'HIOC', t1_ = 't1_hi',
        t1_list = c('t1_vi', 't1_ci', 't1_dep', 't1_fra', 't1_sle'),
        t2_list = c('t2_vi', 't2_ci', 't2_dep', 't2_fra', 't2_sle'),
        oc_list = c('EIOC', 'CIOC', 'DEPOC', 'FRAILOC', 'SLEEPOC'))
FG_test(expo_ = 'CIOC', t1_ = 't1_ci',
        t1_list = c('t1_vi', 't1_hi', 't1_dep', 't1_fra', 't1_sle'),
        t2_list = c('t2_vi', 't2_hi', 't2_dep', 't2_fra', 't2_sle'),
        oc_list = c('EIOC', 'HIOC', 'DEPOC', 'FRAILOC', 'SLEEPOC'))
FG_test(expo_ = 'DEPOC', t1_ = 't1_dep',
        t1_list = c('t1_vi', 't1_hi', 't1_ci', 't1_fra', 't1_sle'),
        t2_list = c('t2_vi', 't2_hi', 't2_ci', 't2_fra', 't2_sle'),
        oc_list = c('EIOC', 'HIOC', 'CIOC', 'FRAILOC', 'SLEEPOC'))
FG_test(expo_ = 'FRAILOC', t1_ = 't1_fra',
        t1_list = c('t1_vi', 't1_hi', 't1_ci', 't1_dep', 't1_sle'),
        t2_list = c('t2_vi', 't2_hi', 't2_ci', 't2_dep', 't2_sle'),
        oc_list = c('EIOC', 'HIOC', 'CIOC', 'DEPOC', 'SLEEPOC'))
FG_test(expo_ = 'SLEEPOC', t1_ = 't1_sle',
        t1_list = c('t1_vi', 't1_hi', 't1_ci', 't1_dep', 't1_fra'),
        t2_list = c('t2_vi', 't2_hi', 't2_ci', 't2_dep', 't2_fra'),
        oc_list = c('EIOC', 'HIOC', 'CIOC', 'DEPOC', 'FRAILOC'))
