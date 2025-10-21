library(dplyr)
library(readxl)
library(nleqslv)

# Загрузка и трансформация данных

col_types = rep("guess", 36)
col_types[2] <-"text"   
col_types[3] <- "text"
col_types[4] <- "text"
col_types[5] <- "numeric"
col_types[10] <- "date"

data_ED <- read_excel("data_ED.xls", col_types = col_types)
data_ED$DEBQsumm <- data_ED$DEBQРезультатЭк + data_ED$Эм +data_ED$Ог

data_men <- data_ED %>% filter(Пол == "Мужской")
data_women <- data_ED %>% filter(Пол == "Женский")


data_DEBQ1 <- data_women %>% filter(DEBQsumm<7)
data_DEBQ2 <- data_women %>% filter((DEBQsumm>=7) &(DEBQsumm<8))
data_DEBQ3 <- data_women %>% filter((DEBQsumm>=8) &(DEBQsumm<9))
data_DEBQ4 <- data_women %>% filter((DEBQsumm>=9) &(DEBQsumm<10))
data_DEBQ5 <- data_women %>% filter(DEBQsumm>=10)

data_age0 <- data_women %>% filter(age_group == 0)
data_age1 <- data_women %>% filter(age_group == 1)
data_age2 <- data_women %>% filter(age_group == 2)
data_age3 <- data_women %>% filter(age_group == 3)



