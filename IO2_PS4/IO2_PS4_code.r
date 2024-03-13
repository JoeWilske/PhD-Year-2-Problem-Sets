# Industrial Organization 2: Problem Set 4
# Author: Joe Wilske
# Date: March 17, 2024


library(dplyr)
library(knitr)
library(kableExtra)
library(ggplot2)


# Read-in dataset
insdata <- read.csv("C:/Gitdir/Projects/PhD-Year-2-Problem-Sets/IO2_PS4/insdata.csv")


#####################################################################################
### --- Task 1: Estimate --- ########################################################
#####################################################################################

# List of all the prices that appear in the dataset in ascending order
P <- insdata$p %>%
  unique() %>%
  sort()

# Creates a vector of the fraction of people who chose H for each price.
find_D <- function(insdata, P) {
  D1 <- c(0,0,0,0,0,0)
  D2 <- c(0,0,0,0,0,0)

  for (i in seq_along(insdata$p)) {

    for (j in seq_along(P)) {

      if (insdata$p[i] == P[j] && insdata$d[i] == 1 ) {
        D1[j] <- D1[j] + 1
      } 
      
      if (insdata$p[i] == P[j]) {
        D2[j] <- D2[j] + 1
      }

    }
  }

  D_fracs <- (D1 / D2)
  insdata <- insdata %>%
    mutate(D = 0)

  for (i in seq_along(insdata$D)) {

    for (j in seq_along(P)) {

      if (insdata$p[i] == P[j]) {
        insdata$D[i] <- D_fracs[j]
      }

    }
  }
  return (insdata)
}

# Augment insdata with newly created D vector
insdata <- find_D(insdata, P)

# Run demand regression and assign variable names to estimates
demand_reg <- lm(D ~ p, data = insdata)
alpha <- unname(coef(demand_reg)[1])
se_alpha <- unname(summary(demand_reg)$coefficients[, "Std. Error"][1])
beta <- unname(coef(demand_reg)[2])
se_beta <- unname(summary(demand_reg)$coefficients[, "Std. Error"][2])

# Create data subsets by plan chosen by consumers
H_data <- insdata[insdata$d == 1, ]
L_data <- insdata[insdata$d == 0, ]

# Run AICH regression and assign variable names
AICH_reg <- lm(0.96*m - 0.9*m ~ p, data = H_data)
gamma_H <- unname(coef(AICH_reg)[1])
se_gamma_H <- unname(summary(AICH_reg)$coefficients[, "Std. Error"][1])
delta_H <- unname(coef(AICH_reg)[2])
se_delta_H <- unname(summary(AICH_reg)$coefficients[, "Std. Error"][2])

# Run AICL regression and assign variable names
AICL_reg <- lm(0.96*m - 0.9*m ~ p, data = L_data)
gamma_L <- unname(coef(AICL_reg)[1])
se_gamma_L <- unname(summary(AICL_reg)$coefficients[, "Std. Error"][1]) 
delta_L <- unname(coef(AICL_reg)[2])
se_delta_L <- unname(summary(AICL_reg)$coefficients[, "Std. Error"][2])

find_AC <- function(insdata, P) {

  ach <- c(0,0,0,0,0,0)
  h_count <- c(0,0,0,0,0,0)
  acl <- c(0,0,0,0,0,0)
  l_count <- c(0,0,0,0,0,0)

  for (i in seq_along(insdata$p)) {

    for (j in seq_along(P)) {

      if (insdata$d[i] == 1 && insdata$p[i] == P[j]) {

        ach[j] <- ach[j] + 0.96*insdata$m[i]
        h_count[j] <- h_count[j] + 1

      } else if (insdata$d[i] == 0 && insdata$p[i] == P[j]) {

        acl[j] <- acl[j] + 0.9*insdata$m[i]
        l_count[j] <- l_count[j] + 1
      }
    }
  }
  ach_fracs <- ach / h_count
  acl_fracs <- acl / l_count
  insdata <- insdata %>%
    mutate(ach = 0) %>%
    mutate(acl = 0)

  for (i in seq_along(insdata$p)) {

    for (j in seq_along(P)) {

      if (insdata$p[i] == P[j]) {

        insdata$ach[i] <- ach_fracs[j]
        insdata$acl[i] <- acl_fracs[j]
      }
    }
  }
  return (insdata)
}

# Add average cost columns,total cost column, and price^2 column to dataset
insdata <- find_AC(insdata, P) %>%
  mutate( tc = ach*D + acl*(1 - D) ) %>%
  mutate( p2 = p^2 )

# Run TCS regression and assign variable names
TCS_reg <- lm(tc ~ p + p2, data = insdata)
phi1 <- unname(coef(TCS_reg)[1])
se_phi1 <- unname(summary(TCS_reg)$coefficients[, "Std. Error"][1])
phi2 <- unname(coef(TCS_reg)[2])
se_phi2 <- unname(summary(TCS_reg)$coefficients[, "Std. Error"][2])
phi3 <- unname(coef(TCS_reg)[3])
se_phi3 <- unname(summary(TCS_reg)$coefficients[, "Std. Error"][3])

# Calculate marginal cost parameters
mu_H <- gamma_H + delta_H*(alpha / beta)
nu_H <- 2*delta_H
mu_L <- gamma_L + delta_L*(alpha / beta)
nu_L <- 2*delta_L
mu_S <- phi2 / beta
nu_S <- 2*phi3 / beta

# Create dataframe for construction of table
table_data <- data.frame(
    Demand = c(alpha, beta, NA),
    SE_Demand = c(se_alpha, se_beta, NA),
    AICH = c(gamma_H, delta_H, NA),
    SE_AICH = c(se_gamma_H, se_delta_H, NA),
    AICL = c(gamma_L, delta_L, NA),
    SE_AICL = c(se_gamma_L, se_delta_L, NA),
    TCS = c(phi1, phi2, phi3),
    SE_TCS = c(se_phi1, se_phi2, se_phi3),
    MCH = c(mu_H, nu_H, NA),
    MCL = c(mu_L, nu_L, NA),
    MCS = c(mu_S, nu_S, NA)
) %>% t()

# Assign row and column names to table_data
rownames(table_data) <- c("Demand", "SE_Demand", "AICH", "SE_AICH", "AICL",
  "SE_AICL", "TCS", "SE_TCS", "MCH", "MCL", "MCS")
colnames(table_data) <- c("Intercept", "Coef_1", "Coef_2")

# Create table
kable(table_data, align = "c", caption = "Estimates and Standard Errors", digits = 3) %>%
    kable_styling(full_width = FALSE) %>%
    row_spec(c(1,3,5,7, 9, 10, 11), bold = TRUE) 


##########################################################################################
### --- Task 2: Initial Plots --- ########################################################
##########################################################################################


ggplot() +
  geom_abline(intercept = - alpha / beta, slope = 1 / beta, color = "#5500ff") + # Demand
  geom_abline(intercept = gamma_H - (delta_H*alpha / beta), slope = delta_H / beta, color = "#ff0000") + # AICH
  geom_abline(intercept = mu_H - (nu_H*alpha / beta), slope = nu_H / beta, color = "#ff02ccbe") + # MCH
  coord_cartesian(xlim = c(0, .5), ylim = c(600, 900)) + # Set window
  labs(title = "Plot of Data and Linear Functions", x = "Quantity", y = "Dollars") + # Add title and axis labels
  theme(aspect.ratio = 1/2, axis.ticks.length  = unit(.5, "cm")) +
  scale_y_continuous(breaks = seq(600, 900, by = 50)) +
  scale_x_continuous(breaks = seq(0, 0.5, by = 0.1))


