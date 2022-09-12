---
title: "01_Time Since Vaccination"
date: "9/12/2023, code originally written 9/21/21"
output: "png files. Rates of VB (cases/100k), Estimated protection against hospitalization
        among vaccine breakthrough"
---
  
library(dplyr)
library(ciisr)
library(stringr)
library(lubridate)
library(tidyr)
library(dbplyr)
library(covidr)
library(ggplot2)

# person-time analysis for vaccine-breakthrough cases, only 
# categories are 2 weeks - 3 mo out from final dose, 
# 3 mo - 6 mo out 
# 6+ mo out 

# stratify by age group 
# 12-17, 18-39, 40-49, 50-59, 60-69, 70-79, 80+ 


connect("covid19", server = "DPHE144")
cv <- connect("covid_vaccine", global_env = F)

# population of fully vaccinated individuals 
iz <- tbl(cv, in_schema("tab", "Patient_UTD_Status")) %>%
  filter(UTD_Flag == 1) %>% 
  select(patient_id, UTD_on) %>% 
  left_join(tbl(cv, in_schema("tab", "LPHA_Patients")) %>% 
              select(patient_id, age_at_1stvaccination) %>% 
              distinct()) %>% 
  rename_all(tolower) %>%
  collect() %>%
  mutate(utd_on = ymd(str_sub(utd_on, 1, 10)))

# .. ----------------------------------------------------------------------

# CIIS patient_ids for breakthrough cases 
breakthrough_cases <- tbl(conn, in_schema("cases", "covid19_cedrs_dashboard_constrained")) %>% 
  filter(breakthrough == 1) %>% 
  select(eventid, earliest_collectiondate, breakthrough) %>% 
  left_join(tbl(conn, in_schema("ciis", "case_patientids")) %>% 
              select(eventid, patient_id)) %>% 
  filter(!is.na(patient_id)) %>% 
  collect()

iz2 <- iz %>% 
  left_join(breakthrough_cases) %>% 
  select(-eventid) %>% 
  mutate(earliest_collectiondate = as.character(earliest_collectiondate), 
         final_date = ifelse(is.na(earliest_collectiondate), 
                             as.character(Sys.Date() - days(8)), 
                             earliest_collectiondate)) %>% 
  mutate_at(vars(final_date, earliest_collectiondate), ymd) %>% 
  mutate(fully_vax = utd_on + days(14), 
         threemo = utd_on + days(90), 
         sixmo = utd_on + days(180)) %>% 
  select(patient_id, fully_vax, threemo, sixmo, final_date, 
         age_at_1stvaccination, breakthrough) %>% 
  filter(fully_vax <= final_date)

iz3 <- iz2 %>% 
  mutate(replace_threemo = threemo > final_date, 
         replace_sixmo = sixmo > final_date) %>% 
  mutate_at(vars(threemo, sixmo, final_date), as.character) %>% 
  mutate(threemo = ifelse(replace_threemo, final_date, threemo), 
         sixmo = ifelse(replace_sixmo, final_date, sixmo)) %>% 
  mutate_at(vars(threemo, sixmo, final_date), ymd) %>% 
  select(-replace_threemo, -replace_sixmo)

iz4 <- iz3 %>% 
  mutate(days_to_2wks_3mo = interval(fully_vax, threemo) / days(1), 
         days_to_3mo_6mo = interval(threemo, sixmo) / days(1), 
         days_to_6moplus = interval(sixmo, final_date) / days(1), 
         age_group = case_when(
           age_at_1stvaccination >= 12 & age_at_1stvaccination <= 17 ~ "12-17", 
           age_at_1stvaccination >= 18 & age_at_1stvaccination <= 39 ~ "18-39", 
           age_at_1stvaccination >= 40 & age_at_1stvaccination <= 49 ~ "40-49", 
           age_at_1stvaccination >= 50 & age_at_1stvaccination <= 59 ~ "50-59", 
           age_at_1stvaccination >= 60 & age_at_1stvaccination <= 69 ~ "60-69", 
           age_at_1stvaccination >= 70 & age_at_1stvaccination <= 79 ~ "70-79", 
           age_at_1stvaccination >= 80 ~ "80+", 
         )) %>% 
  filter(!is.na(age_group))

numerator <- iz4 %>% 
  filter(breakthrough) %>% 
  mutate(group = case_when(
    final_date <= threemo ~ "persontime2wks_3mo", 
    final_date > threemo & final_date <= sixmo ~ "persontime3mo_6mo", 
    final_date > sixmo ~ "persontime6moplus"
  ))

# 2weeks-3mo VB cases / 2weeks-3mo fully vaccinated person-days 
# 3mo-6mo VB cases / 3mo-6mo fully vaccinated person-days
# 6mo+ VB cases / 6mo+ fully vaccinated person-days

final <- iz4 %>% 
  group_by(age_group) %>% 
  summarise(persontime2wks_3mo = sum(days_to_2wks_3mo), 
            persontime3mo_6mo = sum(days_to_3mo_6mo), 
            persontime6moplus = sum(days_to_6moplus)) %>% 
  gather("group", "denominator", -age_group) %>% 
  arrange(age_group) %>% 
  left_join(numerator %>% 
              group_by(age_group) %>% 
              count(group) %>% 
              rename(numerator = n) %>% 
              ungroup()) %>% 
  mutate(rate = (numerator / denominator)*100000) %>% 
  mutate(group = case_when(
    group == "persontime2wks_3mo" ~ "0-3 months", 
    group == "persontime3mo_6mo" ~ "3-6 months", 
    group == "persontime6moplus" ~ "6+ months"
  ), 
  rate_rounded = round(rate, 1))

protection <- final %>% 
  select(age_group, group, rate) %>% 
  spread(group, rate) %>% 
  mutate(`Protection of 0-3mo vs. 6+mo` = (`6+ months` - `0-3 months`) / `6+ months`, 
         `Protection of 3-6mo vs. 6+mo` = (`3-6 months` - `0-3 months`) / `3-6 months`) %>% 
  select(age_group, `Protection of 0-3mo vs. 6+mo`, `Protection of 3-6mo vs. 6+mo`) %>% 
  gather("group", "protection", -age_group) %>% 
  mutate(p2 = paste0(round(protection*100, 1), "%"))

cols <- rev(covid_sequential(color = "blue", n = 3))
cols[1] <- "#A4AECB"
names(cols) <- c("0-3 months", "3-6 months", "6+ months")

final %>% 
  rename(`Time since final vaccination` = group) %>% 
  ggplot(aes(x = age_group, y = rate, fill = `Time since final vaccination`)) + 
  geom_col(position = "dodge") + 
  geom_text(aes(label = rate_rounded, color = `Time since final vaccination`), 
            position = position_dodge(width = 1), 
            vjust = -0.75, size = 3, 
            fontface = "bold") + 
  theme_covid() + 
  scale_color_manual(values = cols) + 
  scale_fill_manual(values = cols) + 
  scale_y_continuous(limits = c(0, max(final$rate) + 5)) + 
  theme(legend.position = "top", legend.justification = "left", 
        legend.text = element_text(face = "bold")) + 
  labs(x = "Age Group", y = "Rate of\nBreakthrough\nCases\nper\n100,000")

ggsave("~/../Downloads/rates.png", width = 10.5, height = 4.5)

cols2 <- rev(covid_sequential(color = "blue", n = 3))[2:3]
names(cols2) <- c("Protection of 0-3mo vs. 6+mo", 
                  "Protection of 3-6mo vs. 6+mo")

protection %>% 
  ggplot(aes(x = age_group, y = protection, fill = group)) + 
  geom_col(position = "dodge", width = 0.75) + 
  geom_text(aes(label = p2, color = group), 
            position = position_dodge(width = 1), 
            vjust = -0.75, size = 3, 
            fontface = "bold") + 
  theme_covid() + 
  scale_color_manual(values = cols2) + 
  scale_fill_manual(values = cols2) + 
  scale_y_continuous(limits = c(0, 1), 
                     labels = scales::percent) + 
  theme(legend.position = "top", legend.justification = "left", 
        legend.text = element_text(face = "bold")) + 
  labs(x = "Age Group", y = "Estimated\nProtection\nAgainst\nVaccine\nBreakthrough", 
       color = "", fill = "")

ggsave("~/../Downloads/protection.png", width = 10.5, height = 4.5)

# move breakthroughs up to July 1 ------------------------------------------

# CIIS patient_ids for breakthrough cases 
breakthrough_cases <- tbl(conn, in_schema("cases", "covid19_cedrs_dashboard_constrained")) %>% 
  filter(breakthrough == 1) %>% 
  select(eventid, earliest_collectiondate, breakthrough) %>% 
  left_join(tbl(conn, in_schema("ciis", "case_patientids")) %>% 
              select(eventid, patient_id)) %>% 
  filter(earliest_collectiondate >= "2021-07-01") %>% 
  collect()

iz2 <- iz %>% 
  left_join(breakthrough_cases) %>% 
  select(-eventid) %>% 
  mutate(earliest_collectiondate = as.character(earliest_collectiondate), 
         final_date = ifelse(is.na(earliest_collectiondate), 
                             as.character(Sys.Date() - days(8)), 
                             earliest_collectiondate)) %>% 
  mutate_at(vars(final_date, earliest_collectiondate), ymd) %>% 
  mutate(fully_vax = utd_on + days(14), 
         threemo = utd_on + days(90), 
         sixmo = utd_on + days(180)) %>% 
  # filter(fully_vax >= ymd("2021-07-01")) %>% 
  select(patient_id, fully_vax, threemo, sixmo, final_date, 
         age_at_1stvaccination, breakthrough) %>% 
  filter(fully_vax <= final_date)

iz3 <- iz2 %>% 
  mutate(replace_threemo = threemo > final_date, 
         replace_sixmo = sixmo > final_date) %>% 
  mutate_at(vars(threemo, sixmo, final_date), as.character) %>% 
  mutate(threemo = ifelse(replace_threemo, final_date, threemo), 
         sixmo = ifelse(replace_sixmo, final_date, sixmo)) %>% 
  mutate_at(vars(threemo, sixmo, final_date), ymd) %>% 
  select(-replace_threemo, -replace_sixmo)

iz4 <- iz3 %>% 
  mutate(days_to_2wks_3mo = interval(fully_vax, threemo) / days(1), 
         days_to_3mo_6mo = interval(threemo, sixmo) / days(1), 
         days_to_6moplus = interval(sixmo, final_date) / days(1), 
         age_group = case_when(
           age_at_1stvaccination >= 12 & age_at_1stvaccination <= 17 ~ "12-17", 
           age_at_1stvaccination >= 18 & age_at_1stvaccination <= 39 ~ "18-39", 
           age_at_1stvaccination >= 40 & age_at_1stvaccination <= 49 ~ "40-49", 
           age_at_1stvaccination >= 50 & age_at_1stvaccination <= 59 ~ "50-59", 
           age_at_1stvaccination >= 60 & age_at_1stvaccination <= 69 ~ "60-69", 
           age_at_1stvaccination >= 70 & age_at_1stvaccination <= 79 ~ "70-79", 
           age_at_1stvaccination >= 80 ~ "80+", 
         )) %>% 
  filter(!is.na(age_group))

numerator <- iz4 %>% 
  filter(breakthrough) %>% 
  mutate(group = case_when(
    final_date <= threemo ~ "persontime2wks_3mo", 
    final_date > threemo & final_date <= sixmo ~ "persontime3mo_6mo", 
    final_date > sixmo ~ "persontime6moplus"
  ))

# 2weeks-3mo VB cases / 2weeks-3mo fully vaccinated person-days 
# 3mo-6mo VB cases / 3mo-6mo fully vaccinated person-days
# 6mo+ VB cases / 6mo+ fully vaccinated person-days

final <- iz4 %>% 
  group_by(age_group) %>% 
  summarise(persontime2wks_3mo = sum(days_to_2wks_3mo), 
            persontime3mo_6mo = sum(days_to_3mo_6mo), 
            persontime6moplus = sum(days_to_6moplus)) %>% 
  gather("group", "denominator", -age_group) %>% 
  left_join(numerator %>% 
              group_by(age_group) %>% 
              count(group) %>% 
              rename(numerator = n) %>% 
              ungroup()) %>% 
  mutate(rate = (numerator / denominator)*100000) %>% 
  mutate(group = case_when(
    group == "persontime2wks_3mo" ~ "0-3 months", 
    group == "persontime3mo_6mo" ~ "3-6 months", 
    group == "persontime6moplus" ~ "6+ months"
  ), 
  rate_rounded = round(rate, 1))

protection <- final %>% 
  select(age_group, group, rate) %>% 
  spread(group, rate) %>% 
  mutate(`Protection of 0-3mo vs. 6+mo` = (`6+ months` - `0-3 months`) / `6+ months`, 
         `Protection of 3-6mo vs. 6+mo` = (`3-6 months` - `0-3 months`) / `3-6 months`) %>% 
  select(age_group, `Protection of 0-3mo vs. 6+mo`, `Protection of 3-6mo vs. 6+mo`) %>% 
  gather("group", "protection", -age_group) %>% 
  mutate(p2 = paste0(round(protection*100, 1), "%"))

cols <- rev(covid_sequential(color = "purple", n = 3))
cols[1] <- "#CAB9C3"
names(cols) <- c("0-3 months", "3-6 months", "6+ months")

final %>% 
  rename(`Time since final vaccination` = group) %>% 
  ggplot(aes(x = age_group, y = rate, fill = `Time since final vaccination`)) + 
  geom_col(position = "dodge") + 
  geom_text(aes(label = rate_rounded, color = `Time since final vaccination`), 
            position = position_dodge(width = 1), 
            vjust = -0.75, size = 3, 
            fontface = "bold") + 
  theme_covid() + 
  scale_color_manual(values = cols) + 
  scale_fill_manual(values = cols) + 
  scale_y_continuous(limits = c(0, max(final$rate) + 5)) + 
  theme(legend.position = "top", legend.justification = "left", 
        legend.text = element_text(face = "bold")) + 
  labs(x = "Age Group", y = "Rate of\nBreakthrough\nCases\nper\n100,000")

ggsave("~/../Downloads/rates2.png", width = 10.5, height = 4.5)

cols2 <- rev(covid_sequential(color = "purple", n = 3))[2:3]
names(cols2) <- c("Protection of 0-3mo vs. 6+mo", 
                  "Protection of 3-6mo vs. 6+mo")

protection %>% 
  ggplot(aes(x = age_group, y = protection, fill = group)) + 
  geom_col(position = "dodge", width = 0.75) + 
  geom_text(aes(label = p2, color = group), 
            position = position_dodge(width = 1), 
            vjust = -0.75, size = 3, 
            fontface = "bold") + 
  theme_covid() + 
  scale_color_manual(values = cols2) + 
  scale_fill_manual(values = cols2) + 
  scale_y_continuous(limits = c(0, 1), 
                     labels = scales::percent) + 
  theme(legend.position = "top", legend.justification = "left", 
        legend.text = element_text(face = "bold")) + 
  labs(x = "Age Group", y = "Estimated\nProtection\nAgainst\nVaccine\nBreakthrough", 
       color = "", fill = "")

ggsave("~/../Downloads/protection2.png", width = 10.5, height = 4.5)


# hospitalizations ------------------------------------------

# CIIS patient_ids for breakthrough cases 
breakthrough_cases <- tbl(conn, in_schema("cases", "covid19_cedrs_dashboard_constrained")) %>% 
  filter(breakthrough == 1, 
         hospitalized_cophs == 1) %>% 
  select(eventid, earliest_collectiondate, breakthrough) %>% 
  left_join(tbl(conn, in_schema("ciis", "case_patientids")) %>% 
              select(eventid, patient_id)) %>% 
  collect()

iz2 <- iz %>% 
  left_join(breakthrough_cases) %>% 
  select(-eventid) %>% 
  mutate(earliest_collectiondate = as.character(earliest_collectiondate), 
         final_date = ifelse(is.na(earliest_collectiondate), 
                             as.character(Sys.Date() - days(8)), 
                             earliest_collectiondate)) %>% 
  mutate_at(vars(final_date, earliest_collectiondate), ymd) %>% 
  mutate(fully_vax = utd_on + days(14), 
         threemo = utd_on + days(90), 
         sixmo = utd_on + days(180)) %>% 
  # filter(fully_vax >= ymd("2021-07-01")) %>% 
  select(patient_id, fully_vax, threemo, sixmo, final_date, 
         age_at_1stvaccination, breakthrough) %>% 
  filter(fully_vax <= final_date)

iz3 <- iz2 %>% 
  mutate(replace_threemo = threemo > final_date, 
         replace_sixmo = sixmo > final_date) %>% 
  mutate_at(vars(threemo, sixmo, final_date), as.character) %>% 
  mutate(threemo = ifelse(replace_threemo, final_date, threemo), 
         sixmo = ifelse(replace_sixmo, final_date, sixmo)) %>% 
  mutate_at(vars(threemo, sixmo, final_date), ymd) %>% 
  select(-replace_threemo, -replace_sixmo)

iz4 <- iz3 %>% 
  mutate(days_to_2wks_3mo = interval(fully_vax, threemo) / days(1), 
         days_to_3mo_6mo = interval(threemo, sixmo) / days(1), 
         days_to_6moplus = interval(sixmo, final_date) / days(1), 
         age_group = case_when(
           age_at_1stvaccination >= 12 & age_at_1stvaccination <= 17 ~ "12-17", 
           age_at_1stvaccination >= 18 & age_at_1stvaccination <= 39 ~ "18-39", 
           age_at_1stvaccination >= 40 & age_at_1stvaccination <= 49 ~ "40-49", 
           age_at_1stvaccination >= 50 & age_at_1stvaccination <= 59 ~ "50-59", 
           age_at_1stvaccination >= 60 & age_at_1stvaccination <= 69 ~ "60-69", 
           age_at_1stvaccination >= 70 & age_at_1stvaccination <= 79 ~ "70-79", 
           age_at_1stvaccination >= 80 ~ "80+", 
         )) %>% 
  filter(!is.na(age_group))

numerator <- iz4 %>% 
  filter(breakthrough) %>% 
  mutate(group = case_when(
    final_date <= threemo ~ "persontime2wks_3mo", 
    final_date > threemo & final_date <= sixmo ~ "persontime3mo_6mo", 
    final_date > sixmo ~ "persontime6moplus"
  ))

# 2weeks-3mo VB cases / 2weeks-3mo fully vaccinated person-days 
# 3mo-6mo VB cases / 3mo-6mo fully vaccinated person-days
# 6mo+ VB cases / 6mo+ fully vaccinated person-days

final <- iz4 %>% 
  group_by(age_group) %>% 
  summarise(persontime2wks_3mo = sum(days_to_2wks_3mo), 
            persontime3mo_6mo = sum(days_to_3mo_6mo), 
            persontime6moplus = sum(days_to_6moplus)) %>% 
  gather("group", "denominator", -age_group) %>% 
  left_join(numerator %>% 
              group_by(age_group) %>% 
              count(group) %>% 
              rename(numerator = n) %>% 
              ungroup()) %>% 
  mutate(rate = (numerator / denominator)*1000000) %>% 
  mutate(group = case_when(
    group == "persontime2wks_3mo" ~ "0-3 months", 
    group == "persontime3mo_6mo" ~ "3-6 months", 
    group == "persontime6moplus" ~ "6+ months"
  ), 
  rate_rounded = round(rate, 2))

protection <- final %>% 
  select(age_group, group, rate) %>% 
  spread(group, rate) %>% 
  mutate(`Protection of 0-3mo vs. 6+mo` = (`6+ months` - `0-3 months`) / `6+ months`, 
         `Protection of 3-6mo vs. 6+mo` = (`3-6 months` - `0-3 months`) / `3-6 months`) %>% 
  select(age_group, `Protection of 0-3mo vs. 6+mo`, `Protection of 3-6mo vs. 6+mo`) %>% 
  gather("group", "protection", -age_group) %>% 
  mutate(p2 = paste0(round(protection*100, 1), "%"))

cols <- rev(covid_sequential(color = "green", n = 3))
cols[1] <- "#C4D4CB"
names(cols) <- c("0-3 months", "3-6 months", "6+ months")

final %>% 
  rename(`Time since final vaccination` = group) %>% 
  ggplot(aes(x = age_group, y = rate, fill = `Time since final vaccination`)) + 
  geom_col(position = "dodge") + 
  geom_text(aes(label = rate_rounded, color = `Time since final vaccination`), 
            position = position_dodge(width = 1), 
            vjust = -0.75, size = 3, 
            fontface = "bold") + 
  theme_covid() + 
  scale_color_manual(values = cols) + 
  scale_fill_manual(values = cols) + 
  scale_y_continuous(limits = c(0, max(final$rate) + 5)) + 
  theme(legend.position = "top", legend.justification = "left", 
        legend.text = element_text(face = "bold")) + 
  labs(x = "Age Group", y = "Rate of\nBreakthrough\nCases\nper\n1 Million")

ggsave("~/../Downloads/rates3.png", width = 10.5, height = 4.5)

cols2 <- rev(covid_sequential(color = "green", n = 3))[2:3]
names(cols2) <- c("Protection of 0-3mo vs. 6+mo", 
                  "Protection of 3-6mo vs. 6+mo")

protection %>% 
  ggplot(aes(x = age_group, y = protection, fill = group)) + 
  geom_col(position = "dodge", width = 0.75) + 
  geom_text(aes(label = p2, color = group), 
            position = position_dodge(width = 1), 
            vjust = -0.75, size = 3, 
            fontface = "bold") + 
  theme_covid() + 
  scale_color_manual(values = cols2) + 
  scale_fill_manual(values = cols2) + 
  scale_y_continuous(limits = c(0, 1), 
                     labels = scales::percent) + 
  theme(legend.position = "top", legend.justification = "left", 
        legend.text = element_text(face = "bold")) + 
  labs(x = "Age Group", y = "Estimated\nProtection\nAgainst\nVaccine\nBreakthrough", 
       color = "", fill = "")

ggsave("~/../Downloads/protection3.png", width = 10.5, height = 4.5)

# hospitalizations since july 1 ------------------------------------------

# CIIS patient_ids for breakthrough cases 
breakthrough_cases <- tbl(conn, in_schema("cases", "covid19_cedrs_dashboard_constrained")) %>% 
  filter(breakthrough == 1, 
         hospitalized_cophs == 1) %>% 
  select(eventid, earliest_collectiondate, breakthrough) %>% 
  left_join(tbl(conn, in_schema("ciis", "case_patientids")) %>% 
              select(eventid, patient_id)) %>% 
  filter(earliest_collectiondate >= "2021-07-01") %>% 
  collect()

iz2 <- iz %>% 
  left_join(breakthrough_cases) %>% 
  select(-eventid) %>% 
  mutate(earliest_collectiondate = as.character(earliest_collectiondate), 
         final_date = ifelse(is.na(earliest_collectiondate), 
                             as.character(Sys.Date() - days(8)), 
                             earliest_collectiondate)) %>% 
  mutate_at(vars(final_date, earliest_collectiondate), ymd) %>% 
  mutate(fully_vax = utd_on + days(14), 
         threemo = utd_on + days(90), 
         sixmo = utd_on + days(180)) %>% 
  # filter(fully_vax >= ymd("2021-07-01")) %>% 
  select(patient_id, fully_vax, threemo, sixmo, final_date, 
         age_at_1stvaccination, breakthrough) %>% 
  filter(fully_vax <= final_date)

iz3 <- iz2 %>% 
  mutate(replace_threemo = threemo > final_date, 
         replace_sixmo = sixmo > final_date) %>% 
  mutate_at(vars(threemo, sixmo, final_date), as.character) %>% 
  mutate(threemo = ifelse(replace_threemo, final_date, threemo), 
         sixmo = ifelse(replace_sixmo, final_date, sixmo)) %>% 
  mutate_at(vars(threemo, sixmo, final_date), ymd) %>% 
  select(-replace_threemo, -replace_sixmo)

iz4 <- iz3 %>% 
  mutate(days_to_2wks_3mo = interval(fully_vax, threemo) / days(1), 
         days_to_3mo_6mo = interval(threemo, sixmo) / days(1), 
         days_to_6moplus = interval(sixmo, final_date) / days(1), 
         age_group = case_when(
           age_at_1stvaccination >= 12 & age_at_1stvaccination <= 17 ~ "12-17", 
           age_at_1stvaccination >= 18 & age_at_1stvaccination <= 39 ~ "18-39", 
           age_at_1stvaccination >= 40 & age_at_1stvaccination <= 49 ~ "40-49", 
           age_at_1stvaccination >= 50 & age_at_1stvaccination <= 59 ~ "50-59", 
           age_at_1stvaccination >= 60 & age_at_1stvaccination <= 69 ~ "60-69", 
           age_at_1stvaccination >= 70 & age_at_1stvaccination <= 79 ~ "70-79", 
           age_at_1stvaccination >= 80 ~ "80+", 
         )) %>% 
  filter(!is.na(age_group))

numerator <- iz4 %>% 
  filter(breakthrough) %>% 
  mutate(group = case_when(
    final_date <= threemo ~ "persontime2wks_3mo", 
    final_date > threemo & final_date <= sixmo ~ "persontime3mo_6mo", 
    final_date > sixmo ~ "persontime6moplus"
  ))

# 2weeks-3mo VB cases / 2weeks-3mo fully vaccinated person-days 
# 3mo-6mo VB cases / 3mo-6mo fully vaccinated person-days
# 6mo+ VB cases / 6mo+ fully vaccinated person-days

final <- iz4 %>% 
  group_by(age_group) %>% 
  summarise(persontime2wks_3mo = sum(days_to_2wks_3mo), 
            persontime3mo_6mo = sum(days_to_3mo_6mo), 
            persontime6moplus = sum(days_to_6moplus)) %>% 
  gather("group", "denominator", -age_group) %>% 
  left_join(numerator %>% 
              group_by(age_group) %>% 
              count(group) %>% 
              rename(numerator = n) %>% 
              ungroup()) %>% 
  mutate(rate = (numerator / denominator)*1000000) %>% 
  mutate(group = case_when(
    group == "persontime2wks_3mo" ~ "0-3 months", 
    group == "persontime3mo_6mo" ~ "3-6 months", 
    group == "persontime6moplus" ~ "6+ months"
  ), 
  rate_rounded = round(rate, 2))

protection <- final %>% 
  select(age_group, group, rate) %>% 
  spread(group, rate) %>% 
  mutate(`Protection of 0-3mo vs. 6+mo` = (`6+ months` - `0-3 months`) / `6+ months`, 
         `Protection of 3-6mo vs. 6+mo` = (`3-6 months` - `0-3 months`) / `3-6 months`) %>% 
  select(age_group, `Protection of 0-3mo vs. 6+mo`, `Protection of 3-6mo vs. 6+mo`) %>% 
  gather("group", "protection", -age_group) %>% 
  mutate(p2 = paste0(round(protection*100, 1), "%"))

cols <- rev(covid_sequential(color = "green", n = 3))
cols[1] <- "#C4D4CB"
names(cols) <- c("0-3 months", "3-6 months", "6+ months")

final %>% 
  rename(`Time since final vaccination` = group) %>% 
  ggplot(aes(x = age_group, y = rate, fill = `Time since final vaccination`)) + 
  geom_col(position = "dodge") + 
  geom_text(aes(label = rate_rounded, color = `Time since final vaccination`), 
            position = position_dodge(width = 1), 
            vjust = -0.75, size = 3, 
            fontface = "bold") + 
  theme_covid() + 
  scale_color_manual(values = cols) + 
  scale_fill_manual(values = cols) + 
  scale_y_continuous(limits = c(0, max(final$rate) + 5)) + 
  theme(legend.position = "top", legend.justification = "left", 
        legend.text = element_text(face = "bold")) + 
  labs(x = "Age Group", y = "Rate of\nBreakthrough\nCases\nper\n1 Million")

ggsave("~/../Downloads/rates3_july1.png", width = 10.5, height = 4.5)

# move breakthroughs up to July 1 ------------------------------------------
## by vaccine type 

# population of fully vaccinated individuals 
iz <- tbl(cv, in_schema("tab", "Patient_UTD_Status")) %>%
  filter(UTD_Flag == 1) %>% 
  select(patient_id, UTD_on) %>% 
  left_join(tbl(cv, in_schema("tab", "LPHA_Patients")) %>% 
              select(patient_id, age_at_1stvaccination, 
                     vaccination_date, vaccination_code) %>% 
              distinct() %>% 
              group_by(patient_id) %>% 
              filter(vaccination_date == min(vaccination_date)) %>% 
              ungroup() %>% 
              select(-vaccination_date)) %>% 
  rename_all(tolower) %>%
  collect() %>%
  mutate(utd_on = ymd(str_sub(utd_on, 1, 10)))

iz <- iz %>% 
  mutate(vaccine = case_when(
    vaccination_code == "COVID-19 mRNA (MOD)" ~ "Moderna", 
    vaccination_code == "COVID-19 mRNA (PFR)" ~ "Pfizer", 
    vaccination_code == "COVID-19 Vector-NR (JSN)" ~ "Janssen", 
    T ~ "Other"
  )) %>% 
  select(-vaccination_code)

# CIIS patient_ids for breakthrough cases 
breakthrough_cases <- tbl(conn, in_schema("cases", "covid19_cedrs_dashboard_constrained")) %>% 
  filter(breakthrough == 1) %>% 
  select(eventid, earliest_collectiondate, breakthrough) %>% 
  left_join(tbl(conn, in_schema("ciis", "case_patientids")) %>% 
              select(eventid, patient_id)) %>% 
  filter(earliest_collectiondate >= "2021-07-01") %>% 
  collect() 

iz2 <- iz %>% 
  left_join(breakthrough_cases) %>% 
  select(-eventid) %>% 
  mutate(earliest_collectiondate = as.character(earliest_collectiondate), 
         final_date = ifelse(is.na(earliest_collectiondate), 
                             as.character(Sys.Date() - days(8)), 
                             earliest_collectiondate)) %>% 
  mutate_at(vars(final_date, earliest_collectiondate), ymd) %>% 
  mutate(fully_vax = utd_on + days(14), 
         threemo = utd_on + days(90), 
         sixmo = utd_on + days(180)) %>% 
  # filter(fully_vax >= ymd("2021-07-01")) %>% 
  select(patient_id, fully_vax, threemo, sixmo, final_date, 
         age_at_1stvaccination, breakthrough, vaccine) %>% 
  filter(fully_vax <= final_date)

iz3 <- iz2 %>% 
  mutate(replace_threemo = threemo > final_date, 
         replace_sixmo = sixmo > final_date) %>% 
  mutate_at(vars(threemo, sixmo, final_date), as.character) %>% 
  mutate(threemo = ifelse(replace_threemo, final_date, threemo), 
         sixmo = ifelse(replace_sixmo, final_date, sixmo)) %>% 
  mutate_at(vars(threemo, sixmo, final_date), ymd) %>% 
  select(-replace_threemo, -replace_sixmo)

iz4 <- iz3 %>% 
  mutate(days_to_2wks_3mo = interval(fully_vax, threemo) / days(1), 
         days_to_3mo_6mo = interval(threemo, sixmo) / days(1), 
         days_to_6moplus = interval(sixmo, final_date) / days(1), 
         age_group = case_when(
           age_at_1stvaccination >= 12 & age_at_1stvaccination <= 17 ~ "12-17", 
           age_at_1stvaccination >= 18 & age_at_1stvaccination <= 39 ~ "18-39", 
           age_at_1stvaccination >= 40 & age_at_1stvaccination <= 49 ~ "40-49", 
           age_at_1stvaccination >= 50 & age_at_1stvaccination <= 59 ~ "50-59", 
           age_at_1stvaccination >= 60 & age_at_1stvaccination <= 69 ~ "60-69", 
           age_at_1stvaccination >= 70 & age_at_1stvaccination <= 79 ~ "70-79", 
           age_at_1stvaccination >= 80 ~ "80+", 
         )) %>% 
  filter(!is.na(age_group))

numerator <- iz4 %>% 
  filter(breakthrough) %>% 
  mutate(group = case_when(
    final_date <= threemo ~ "persontime2wks_3mo", 
    final_date > threemo & final_date <= sixmo ~ "persontime3mo_6mo", 
    final_date > sixmo ~ "persontime6moplus"
  ))

# 2weeks-3mo VB cases / 2weeks-3mo fully vaccinated person-days 
# 3mo-6mo VB cases / 3mo-6mo fully vaccinated person-days
# 6mo+ VB cases / 6mo+ fully vaccinated person-days

final <- iz4 %>% 
  group_by(age_group, vaccine) %>% 
  summarise(persontime2wks_3mo = sum(days_to_2wks_3mo), 
            persontime3mo_6mo = sum(days_to_3mo_6mo), 
            persontime6moplus = sum(days_to_6moplus)) %>% 
  ungroup() %>% 
  gather("group", "denominator", -age_group, -vaccine) %>% 
  left_join(numerator %>% 
              group_by(age_group, vaccine) %>% 
              count(group) %>% 
              rename(numerator = n) %>% 
              ungroup()) %>% 
  mutate(numerator = ifelse(is.na(numerator), 0, numerator)) %>% 
  mutate(rate = (numerator / denominator)*100000) %>% 
  mutate(group = case_when(
    group == "persontime2wks_3mo" ~ "0-3 months", 
    group == "persontime3mo_6mo" ~ "3-6 months", 
    group == "persontime6moplus" ~ "6+ months"
  ), 
  rate_rounded = round(rate, 1)) 

protection <- final %>% 
  select(age_group, vaccine, group, rate) %>% 
  spread(group, rate) %>% 
  mutate(`Protection of 0-3mo vs. 6+mo` = (`6+ months` - `0-3 months`) / `6+ months`, 
         `Protection of 0-3mo vs. 3-6mo` = (`3-6 months` - `0-3 months`) / `3-6 months`) %>% 
  mutate_at(vars(`Protection of 0-3mo vs. 6+mo`, `Protection of 0-3mo vs. 3-6mo`), 
            function(x) ifelse(is.nan(x) | is.infinite(x), NA, x)) %>% 
  select(age_group, vaccine, `Protection of 0-3mo vs. 6+mo`, `Protection of 0-3mo vs. 3-6mo`) %>% 
  gather("group", "protection", -age_group, -vaccine) %>% 
  mutate(p2 = paste0(round(protection*100, 1), "%"), 
         p2 = ifelse(p2 == "NA%", NA, p2)) %>% 
  mutate(group = factor(group, 
                        levels = c("Protection of 0-3mo vs. 3-6mo", 
                                   "Protection of 0-3mo vs. 6+mo")))

cols <- rev(covid_sequential(color = "purple", n = 3))
cols[1] <- "#CAB9C3"
names(cols) <- c("0-3 months", "3-6 months", "6+ months")

cols2 <- rev(covid_sequential(color = "purple", n = 3))[2:3]
names(cols2) <- c("Protection of 0-3mo vs. 6+mo", 
                  "Protection of 0-3mo vs. 3-6mo")

f2 <- final %>% 
  filter(vaccine == "Moderna")

f2 %>% 
  rename(`Time since final vaccination` = group) %>% 
  ggplot(aes(x = age_group, y = rate, fill = `Time since final vaccination`)) + 
  geom_col(position = "dodge") + 
  geom_text(aes(label = rate_rounded, color = `Time since final vaccination`), 
            position = position_dodge(width = 1), 
            vjust = -0.75, size = 3, 
            fontface = "bold") + 
  theme_covid() + 
  scale_color_manual(values = cols) + 
  scale_fill_manual(values = cols) + 
  scale_y_continuous(limits = c(0, max(f2$rate) + 5)) + 
  labs(x = "Age Group", y = "Rate of\nBreakthrough\nCases\nper\n100,000", 
       title = "Person-time analysis of vaccine breakthrough since July 1st", 
       subtitle = "by time since vaccination and age group (Moderna)") 

ggsave("~/../Downloads/rates2_moderna.png", width = 10.5, height = 4.5)

protection %>% 
  filter(vaccine == "Moderna") %>% 
  ggplot(aes(x = age_group, y = protection, fill = group)) + 
  geom_col(position = "dodge", width = 0.75) + 
  geom_text(aes(label = p2, color = group), 
            position = position_dodge(width = 1), 
            vjust = -0.75, size = 3, 
            fontface = "bold") + 
  theme_covid() + 
  scale_color_manual(values = cols2) + 
  scale_fill_manual(values = cols2) + 
  scale_y_continuous(limits = c(0, 1), 
                     labels = scales::percent) + 
  theme(legend.position = "top", legend.justification = "left", 
        legend.text = element_text(face = "bold")) + 
  labs(x = "Age Group", y = "Estimated\nProtection\nAgainst\nVaccine\nBreakthrough", 
       color = "", fill = "", 
       title = "Estimated protection against vaccine breakthrough since July 1st", 
       subtitle = "based on time since vaccination and age group (Moderna)")

ggsave("~/../Downloads/protection2_moderna.png", width = 10.5, height = 4.5)

f2 <- final %>% 
  filter(vaccine == "Pfizer")

f2 %>% 
  rename(`Time since final vaccination` = group) %>% 
  ggplot(aes(x = age_group, y = rate, fill = `Time since final vaccination`)) + 
  geom_col(position = "dodge") + 
  geom_text(aes(label = rate_rounded, color = `Time since final vaccination`), 
            position = position_dodge(width = 1), 
            vjust = -0.75, size = 3, 
            fontface = "bold") + 
  theme_covid() + 
  scale_color_manual(values = cols) + 
  scale_fill_manual(values = cols) + 
  scale_y_continuous(limits = c(0, max(f2$rate) + 5)) + 
  labs(x = "Age Group", y = "Rate of\nBreakthrough\nCases\nper\n100,000", 
       title = "Person-time analysis of vaccine breakthrough since July 1st", 
       subtitle = "by time since vaccination and age group (Pfizer)") 

ggsave("~/../Downloads/rates2_pfizer.png", width = 10.5, height = 4.5)

protection %>% 
  filter(vaccine == "Pfizer") %>% 
  ggplot(aes(x = age_group, y = protection, fill = group)) + 
  geom_col(position = "dodge", width = 0.75) + 
  geom_text(aes(label = p2, color = group), 
            position = position_dodge(width = 1), 
            vjust = -0.75, size = 3, 
            fontface = "bold") + 
  theme_covid() + 
  scale_color_manual(values = cols2) + 
  scale_fill_manual(values = cols2) + 
  scale_y_continuous(limits = c(0, 1), 
                     labels = scales::percent) + 
  theme(legend.position = "top", legend.justification = "left", 
        legend.text = element_text(face = "bold")) + 
  labs(x = "Age Group", y = "Estimated\nProtection\nAgainst\nVaccine\nBreakthrough", 
       color = "", fill = "", 
       title = "Estimated protection against vaccine breakthrough since July 1st", 
       subtitle = "based on time since vaccination and age group (Pfizer)")

ggsave("~/../Downloads/protection2_pfizer.png", width = 10.5, height = 4.5)

f2 <- final %>% 
  filter(vaccine == "Janssen")

f2 %>% 
  rename(`Time since final vaccination` = group) %>% 
  ggplot(aes(x = age_group, y = rate, fill = `Time since final vaccination`)) + 
  geom_col(position = "dodge") + 
  geom_text(aes(label = rate_rounded, color = `Time since final vaccination`), 
            position = position_dodge(width = 1), 
            vjust = -0.75, size = 3, 
            fontface = "bold") + 
  theme_covid() + 
  scale_color_manual(values = cols) + 
  scale_fill_manual(values = cols) + 
  scale_y_continuous(limits = c(0, max(f2$rate) + 5)) + 
  labs(x = "Age Group", y = "Rate of\nBreakthrough\nCases\nper\n100,000", 
       title = "Person-time analysis of vaccine breakthrough since July 1st", 
       subtitle = "by time since vaccination and age group (Janssen)") 

ggsave("~/../Downloads/rates2_janssen.png", width = 10.5, height = 4.5)

protection %>% 
  filter(vaccine == "Janssen") %>% 
  ggplot(aes(x = age_group, y = protection, fill = group)) + 
  geom_col(position = "dodge", width = 0.75) + 
  geom_text(aes(label = p2, color = group), 
            position = position_dodge(width = 1), 
            vjust = -0.75, size = 3, 
            fontface = "bold") + 
  theme_covid() + 
  scale_color_manual(values = cols2) + 
  scale_fill_manual(values = cols2) + 
  scale_y_continuous(limits = c(0, 1), 
                     labels = scales::percent) + 
  theme(legend.position = "top", legend.justification = "left", 
        legend.text = element_text(face = "bold")) + 
  labs(x = "Age Group", y = "Estimated\nProtection\nAgainst\nVaccine\nBreakthrough", 
       color = "", fill = "", 
       title = "Estimated protection against vaccine breakthrough since July 1st", 
       subtitle = "based on time since vaccination and age group (Janssen)")

ggsave("~/../Downloads/protection2_janssen.png", width = 10.5, height = 4.5)

# hospitalizations since july 1 ------------------------------------------
## by vaccine type 

# population of fully vaccinated individuals 
iz <- tbl(cv, in_schema("tab", "Patient_UTD_Status")) %>%
  filter(UTD_Flag == 1) %>% 
  select(patient_id, UTD_on) %>% 
  left_join(tbl(cv, in_schema("tab", "LPHA_Patients")) %>% 
              select(patient_id, age_at_1stvaccination, 
                     vaccination_date, vaccination_code) %>% 
              distinct() %>% 
              group_by(patient_id) %>% 
              filter(vaccination_date == min(vaccination_date)) %>% 
              ungroup() %>% 
              select(-vaccination_date)) %>% 
  rename_all(tolower) %>%
  collect() %>%
  mutate(utd_on = ymd(str_sub(utd_on, 1, 10)))

iz <- iz %>% 
  mutate(vaccine = case_when(
    vaccination_code == "COVID-19 mRNA (MOD)" ~ "Moderna", 
    vaccination_code == "COVID-19 mRNA (PFR)" ~ "Pfizer", 
    vaccination_code == "COVID-19 Vector-NR (JSN)" ~ "Janssen", 
    T ~ "Other"
  )) %>% 
  select(-vaccination_code)

# CIIS patient_ids for breakthrough cases 
breakthrough_cases <- tbl(conn, in_schema("cases", "covid19_cedrs_dashboard_constrained")) %>% 
  filter(breakthrough == 1, 
         hospitalized_cophs == 1) %>% 
  select(eventid, earliest_collectiondate, breakthrough) %>% 
  left_join(tbl(conn, in_schema("ciis", "case_patientids")) %>% 
              select(eventid, patient_id)) %>% 
  filter(earliest_collectiondate >= "2021-07-01") %>% 
  collect()

iz2 <- iz %>% 
  left_join(breakthrough_cases) %>% 
  select(-eventid) %>% 
  mutate(earliest_collectiondate = as.character(earliest_collectiondate), 
         final_date = ifelse(is.na(earliest_collectiondate), 
                             as.character(Sys.Date() - days(8)), 
                             earliest_collectiondate)) %>% 
  mutate_at(vars(final_date, earliest_collectiondate), ymd) %>% 
  mutate(fully_vax = utd_on + days(14), 
         threemo = utd_on + days(90), 
         sixmo = utd_on + days(180)) %>% 
  # filter(fully_vax >= ymd("2021-07-01")) %>% 
  select(patient_id, fully_vax, threemo, sixmo, final_date, 
         age_at_1stvaccination, breakthrough, vaccine) %>% 
  filter(fully_vax <= final_date)

iz3 <- iz2 %>% 
  mutate(replace_threemo = threemo > final_date, 
         replace_sixmo = sixmo > final_date) %>% 
  mutate_at(vars(threemo, sixmo, final_date), as.character) %>% 
  mutate(threemo = ifelse(replace_threemo, final_date, threemo), 
         sixmo = ifelse(replace_sixmo, final_date, sixmo)) %>% 
  mutate_at(vars(threemo, sixmo, final_date), ymd) %>% 
  select(-replace_threemo, -replace_sixmo)

iz4 <- iz3 %>% 
  mutate(days_to_2wks_3mo = interval(fully_vax, threemo) / days(1), 
         days_to_3mo_6mo = interval(threemo, sixmo) / days(1), 
         days_to_6moplus = interval(sixmo, final_date) / days(1), 
         age_group = case_when(
           age_at_1stvaccination >= 12 & age_at_1stvaccination <= 17 ~ "12-17", 
           age_at_1stvaccination >= 18 & age_at_1stvaccination <= 39 ~ "18-39", 
           age_at_1stvaccination >= 40 & age_at_1stvaccination <= 49 ~ "40-49", 
           age_at_1stvaccination >= 50 & age_at_1stvaccination <= 59 ~ "50-59", 
           age_at_1stvaccination >= 60 & age_at_1stvaccination <= 69 ~ "60-69", 
           age_at_1stvaccination >= 70 & age_at_1stvaccination <= 79 ~ "70-79", 
           age_at_1stvaccination >= 80 ~ "80+", 
         )) %>% 
  filter(!is.na(age_group))

numerator <- iz4 %>% 
  filter(breakthrough) %>% 
  mutate(group = case_when(
    final_date <= threemo ~ "persontime2wks_3mo", 
    final_date > threemo & final_date <= sixmo ~ "persontime3mo_6mo", 
    final_date > sixmo ~ "persontime6moplus"
  ))

# 2weeks-3mo VB cases / 2weeks-3mo fully vaccinated person-days 
# 3mo-6mo VB cases / 3mo-6mo fully vaccinated person-days
# 6mo+ VB cases / 6mo+ fully vaccinated person-days

final <- iz4 %>% 
  group_by(age_group, vaccine) %>% 
  summarise(persontime2wks_3mo = sum(days_to_2wks_3mo), 
            persontime3mo_6mo = sum(days_to_3mo_6mo), 
            persontime6moplus = sum(days_to_6moplus)) %>% 
  ungroup() %>% 
  gather("group", "denominator", -age_group, -vaccine) %>% 
  left_join(numerator %>% 
              group_by(age_group, vaccine) %>% 
              count(group) %>% 
              rename(numerator = n) %>% 
              ungroup()) %>% 
  mutate(numerator = ifelse(is.na(numerator), 0, numerator)) %>% 
  mutate(rate = (numerator / denominator)*1000000) %>% 
  mutate(group = case_when(
    group == "persontime2wks_3mo" ~ "0-3 months", 
    group == "persontime3mo_6mo" ~ "3-6 months", 
    group == "persontime6moplus" ~ "6+ months"
  ), 
  rate_rounded = round(rate, 2))

protection <- final %>% 
  select(age_group, vaccine, group, rate) %>% 
  spread(group, rate) %>% 
  mutate(`Protection of 0-3mo vs. 6+mo` = (`6+ months` - `0-3 months`) / `6+ months`, 
         `Protection of 0-3mo vs. 3-6mo` = (`3-6 months` - `0-3 months`) / `3-6 months`) %>% 
  mutate_at(vars(`Protection of 0-3mo vs. 6+mo`, `Protection of 0-3mo vs. 3-6mo`), 
            function(x) ifelse(is.nan(x) | is.infinite(x), NA, x)) %>% 
  select(age_group, vaccine, `Protection of 0-3mo vs. 6+mo`, `Protection of 0-3mo vs. 3-6mo`) %>% 
  gather("group", "protection", -age_group, -vaccine) %>% 
  mutate(p2 = paste0(round(protection*100, 1), "%"), 
         p2 = ifelse(p2 == "NA%", NA, p2)) # %>% 
# mutate(group = factor(group, 
#                       levels = c("Protection of 0-3mo vs. 3-6mo", 
#                                  "Protection of 0-3mo vs. 6+mo")))

cols <- rev(covid_sequential(color = "green", n = 3))
cols[1] <- "#C4D4CB"
names(cols) <- c("0-3 months", "3-6 months", "6+ months")

cols2 <- rev(covid_sequential(color = "green", n = 3))[2:3]
names(cols2) <- c("Protection of 0-3mo vs. 6+mo", 
                  "Protection of 0-3mo vs. 3-6mo")

f2 <- final %>% 
  filter(vaccine == "Moderna")

f2 %>% 
  rename(`Time since final vaccination` = group) %>% 
  ggplot(aes(x = age_group, y = rate, fill = `Time since final vaccination`)) + 
  geom_col(position = "dodge") + 
  geom_text(aes(label = rate_rounded, color = `Time since final vaccination`), 
            position = position_dodge(width = 1), 
            vjust = -0.75, size = 3, 
            fontface = "bold") + 
  theme_covid() + 
  scale_color_manual(values = cols) + 
  scale_fill_manual(values = cols) + 
  scale_y_continuous(limits = c(0, max(f2$rate) + 5)) + 
  theme(legend.position = "top", legend.justification = "left", 
        legend.text = element_text(face = "bold")) + 
  labs(x = "Age Group", y = "Rate of\nBreakthrough\nCases\nper\n1 Million", 
       title = "Person-time analysis of hospitalized vaccine breakthrough since July 1", 
       subtitle = "by time since vaccination and age group (Moderna)")

ggsave("~/../Downloads/rates3_july1_moderna.png", width = 10.5, height = 4.5)

protection %>% 
  filter(vaccine == "Moderna") %>% 
  ggplot(aes(x = age_group, y = protection, fill = group)) + 
  geom_col(position = "dodge", width = 0.75) + 
  geom_text(aes(label = p2, color = group), 
            position = position_dodge(width = 1), 
            vjust = -0.75, size = 3, 
            fontface = "bold") + 
  theme_covid() + 
  scale_color_manual(values = cols2) + 
  scale_fill_manual(values = cols2) + 
  scale_y_continuous(limits = c(0, 1), 
                     labels = scales::percent) + 
  theme(legend.position = "top", legend.justification = "left", 
        legend.text = element_text(face = "bold")) + 
  labs(x = "Age Group", y = "Estimated\nProtection\nAgainst\nVaccine\nBreakthrough", 
       color = "", fill = "", 
       title = "Estimated protection against hospitalization among vaccine breakthrough", 
       subtitle = "based on time since vaccination by age group (Moderna)")

ggsave("~/../Downloads/protection3_moderna.png", width = 10.5, height = 4.5)

f2 <- final %>% 
  filter(vaccine == "Pfizer")

f2 %>% 
  rename(`Time since final vaccination` = group) %>% 
  ggplot(aes(x = age_group, y = rate, fill = `Time since final vaccination`)) + 
  geom_col(position = "dodge") + 
  geom_text(aes(label = rate_rounded, color = `Time since final vaccination`), 
            position = position_dodge(width = 1), 
            vjust = -0.75, size = 3, 
            fontface = "bold") + 
  theme_covid() + 
  scale_color_manual(values = cols) + 
  scale_fill_manual(values = cols) + 
  scale_y_continuous(limits = c(0, max(f2$rate) + 5)) + 
  theme(legend.position = "top", legend.justification = "left", 
        legend.text = element_text(face = "bold")) + 
  labs(x = "Age Group", y = "Rate of\nBreakthrough\nCases\nper\n1 Million", 
       title = "Person-time analysis of hospitalized vaccine breakthrough since July 1", 
       subtitle = "by time since vaccination and age group (Pfizer)")

ggsave("~/../Downloads/rates3_july1_pfizer.png", width = 10.5, height = 4.5)

protection %>% 
  filter(vaccine == "Pfizer") %>% 
  ggplot(aes(x = age_group, y = protection, fill = group)) + 
  geom_col(position = "dodge", width = 0.75) + 
  geom_text(aes(label = p2, color = group), 
            position = position_dodge(width = 1), 
            vjust = -0.75, size = 3, 
            fontface = "bold") + 
  theme_covid() + 
  scale_color_manual(values = cols2) + 
  scale_fill_manual(values = cols2) + 
  scale_y_continuous(limits = c(0, 1), 
                     labels = scales::percent) + 
  theme(legend.position = "top", legend.justification = "left", 
        legend.text = element_text(face = "bold")) + 
  labs(x = "Age Group", y = "Estimated\nProtection\nAgainst\nVaccine\nBreakthrough", 
       color = "", fill = "", 
       title = "Estimated protection against hospitalization among vaccine breakthrough", 
       subtitle = "based on time since vaccination by age group (Pfizer)")

ggsave("~/../Downloads/protection3_pfizer.png", width = 10.5, height = 4.5)

f2 <- final %>% 
  filter(vaccine == "Janssen")

f2 %>% 
  rename(`Time since final vaccination` = group) %>% 
  ggplot(aes(x = age_group, y = rate, fill = `Time since final vaccination`)) + 
  geom_col(position = "dodge") + 
  geom_text(aes(label = rate_rounded, color = `Time since final vaccination`), 
            position = position_dodge(width = 1), 
            vjust = -0.75, size = 3, 
            fontface = "bold") + 
  theme_covid() + 
  scale_color_manual(values = cols) + 
  scale_fill_manual(values = cols) + 
  scale_y_continuous(limits = c(0, max(f2$rate) + 5)) + 
  theme(legend.position = "top", legend.justification = "left", 
        legend.text = element_text(face = "bold")) + 
  labs(x = "Age Group", y = "Rate of\nBreakthrough\nCases\nper\n1 Million", 
       title = "Person-time analysis of hospitalized vaccine breakthrough since July 1", 
       subtitle = "by time since vaccination and age group (Janssen)")

ggsave("~/../Downloads/rates3_july1_janssen.png", width = 10.5, height = 4.5)

protection %>% 
  filter(vaccine == "Janssen") %>% 
  ggplot(aes(x = age_group, y = protection, fill = group)) + 
  geom_col(position = "dodge", width = 0.75) + 
  geom_text(aes(label = p2, color = group), 
            position = position_dodge(width = 1), 
            vjust = -0.75, size = 3, 
            fontface = "bold") + 
  theme_covid() + 
  scale_color_manual(values = cols2) + 
  scale_fill_manual(values = cols2) + 
  scale_y_continuous(limits = c(0, 1), 
                     labels = scales::percent) + 
  theme(legend.position = "top", legend.justification = "left", 
        legend.text = element_text(face = "bold")) + 
  labs(x = "Age Group", y = "Estimated\nProtection\nAgainst\nVaccine\nBreakthrough", 
       color = "", fill = "", 
       title = "Estimated protection against hospitalization among vaccine breakthrough", 
       subtitle = "based on time since vaccination by age group (Janssen)")

ggsave("~/../Downloads/protection3_janssen.png", width = 10.5, height = 4.5)
