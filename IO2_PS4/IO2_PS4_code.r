# Industrial Organization 2: Problem Set 4
# Author: Joe Wilske
# Date: March 17, 2024

library(dplyr)
library(knitr)
library(kableExtra)

insdata <- read.csv("C:/Gitdir/Projects/PhD-Year-2-Problem-Sets/IO2_PS4/insdata.csv")

demand_reg <- lm(d ~ p, data = insdata)
alpha <- unname(coef(demand_reg)[1])
se_alpha <- unname(summary(demand_reg)$coefficients[, "Std. Error"][1])
beta <- unname(coef(demand_reg)[2])
se_beta <- unname(summary(demand_reg)$coefficients[, "Std. Error"][2])

H_data <- insdata[insdata$d == 1, ]
L_data <- insdata[insdata$d == 0, ]

AICH_reg <- lm(0.96*m ~ p, data = H_data)
gamma_H <- unname(coef(AICH_reg)[1])
se_gamma_H <- unname(summary(AICH_reg)$coefficients[, "Std. Error"][1])
delta_H <- unname(coef(AICH_reg)[2])
se_delta_H <- unname(summary(AICH_reg)$coefficients[, "Std. Error"][2])

AICL_reg <- lm(0.9*m ~ p, data = L_data)
gamma_L <- unname(coef(AICL_reg)[1])
se_gamma_L <- unname(summary(AICL_reg)$coefficients[, "Std. Error"][1]) 
delta_L <- unname(coef(AICL_reg)[2])
se_delta_L <- unname(summary(AICL_reg)$coefficients[, "Std. Error"][2])

insdata <- insdata %>%
  mutate( c = 0.96*m*d + 0.9*m*(1 - d) ) %>%
  mutate( p2 = p^2 )

TCS_reg <- lm(c ~ p + p2, data = insdata)
phi1 <- unname(coef(TCS_reg)[1])
se_phi1 <- unname(summary(TCS_reg)$coefficients[, "Std. Error"][1])
phi2 <- unname(coef(TCS_reg)[2])
se_phi2 <- unname(summary(TCS_reg)$coefficients[, "Std. Error"][2])
phi3 <- unname(coef(TCS_reg)[3])
se_phi3 <- unname(summary(TCS_reg)$coefficients[, "Std. Error"][3])

table_data <- data.frame(
    Demand = c(alpha, beta, NA),
    SE_Demand = c(se_alpha, se_beta, NA),
    AICH = c(gamma_H, delta_H, NA),
    SE_AICH = c(se_gamma_H, se_delta_H, NA),
    AICL = c(gamma_L, delta_L, NA),
    SE_AICL = c(se_gamma_L, se_delta_L, NA),
    TCS = c(phi1, phi2, phi3),
    SE_TCS = c(se_phi1, se_phi2, se_phi3)
) %>% t()

rownames(table_data) <- c("Demand", "SE_Demand", "AICH", "SE_AICH", "AICL", "SE_AICL", "TCS", "SE_TCS")
colnames(table_data) <- c("Intercept", "Coef_1", "Coef_2")

kable(table_data, align = "c",
    caption = "Estimates and Standard Errors", digits = 3) %>%
    kable_styling(full_width = FALSE) %>%
    row_spec(c(1,3,5,7), bold = TRUE) 
