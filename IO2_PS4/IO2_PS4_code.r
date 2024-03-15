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
  return (D1 / D2)
}

# Execute above function
D <- find_D(insdata, P)

# Add D to insdata
insdata <- insdata %>%
    mutate(D = 0)
for (i in seq_along(insdata$D)) {

  for (j in seq_along(P)) {

    if (insdata$p[i] == P[j]) {
      insdata$D[i] <- D[j]
    }
  }
}

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

# Create vector of average incremental costs for each price
find_AICH <- function(insdata, P) {

  AICH <- c(0,0,0,0,0,0)
  counter <- c(0,0,0,0,0,0)

  for (i in seq_along(insdata$p)) {

    for (j in seq_along(P)) {

      if (insdata$d[i] == 1 && P[j] == insdata$p[i]) {

        AICH[j] <- AICH[j] + 0.06*insdata$m[i]
        counter[j] <- counter[j] + 1
      }
    }
  }
  return (AICH / counter)
}
# Execute above function
AICH <- find_AICH(insdata, P)

# Create vectors of average total costs for each price. 
# Append to insdata
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

# Execute above function and add total cost and price^2 columns to insdata
insdata <- find_AC(insdata, P) %>%
  mutate( tc = 0.96*m*d + 0.9*m*(1 - d) ) %>%
#  mutate( tc = ach*D + acl*(1 - D) ) %>%
#  mutate( tc = 0.96*m*D + 0.9*(1 - D) ) %>%
#  mutate( tc = ach*d + acl*(1 - d) ) %>%
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
mu_L <- gamma_L + delta_L*((1 - alpha) / beta)
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
kable(table_data, align = "c", caption = "Estimates and Standard Errors", 
      digits = 3) %>%
    kable_styling(full_width = FALSE) %>%
    row_spec(c(1,3,5,7, 9, 10, 11), bold = TRUE)
    

################################################################################################
### --- Task 2: Initial Plots --- ##############################################################
################################################################################################


# Create point-data to graph
point_data <- data.frame(
  prices = P,
  quants = D,
  costs = AICH,
  num_consumers = c(6819, 1282, 257, 4168, 1926, 548)
)

# Set coordinates for DWL area to graph
dwl_coords <- data.frame(
  x = c(0.0, 0.0, .375), 
  y = c(845, 905, 772)
)


# Set graph window parameters
x_min <- 0
x_max <- 1
y_min <- 400
y_max <- 1000

# Create a beautiful graph
ggplot() +
  geom_abline(intercept = - alpha / beta, slope = 1 / beta, color = "#5500ff") + # Demand
  geom_abline(intercept = gamma_H - (delta_H*alpha / beta), slope = delta_H / beta, color = "#ff0000") + # AICH
  geom_abline(intercept = mu_H - (nu_H*alpha / beta), slope = nu_H / beta, color = "#ff02ccbe") + # MCH
  geom_point(data = point_data, aes(x = quants, y = prices, size = num_consumers), color = "#0600ab80", stroke = 2) + # Demand data points
  geom_point(data = point_data, aes(x = quants, y = costs, size = num_consumers), color = "#b6000080", stroke = 2) + # Cost data points
  scale_size_continuous(range = c(2, 10)) + # Control size range of circles
  coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max)) + # Set window
  labs(title = "Plot of Data and Linear Functions", x = "Quantity", y = "Dollars") + # Add title and axis labels
  theme(aspect.ratio = 1/2, axis.ticks.length  = unit(.1, "cm")) + # Control aspect ratio and unit ticks
  scale_y_continuous(breaks = seq(y_min, y_max, by = 100)) + # Set y-axis number labels
  scale_x_continuous(breaks = seq(x_min, x_max, by = 0.1)) + # Set x-axis number labels
  geom_text(aes(x = 0.05, y = 925, label = "Demand"), color = "#5500ff") +
  geom_text(aes(x = 0.8, y = 725, label = "AIC_H"), color = "#ff0000") +
  geom_text(aes(x = 0.92, y = 460, label = "MC_H"), color = "#ff02ccbe") +
  geom_polygon(data = dwl_coords, aes(x = x, y = y), fill = "#d4ff0080") +
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 1100), color = "black", linetype = "solid") +
  geom_segment(aes(x = 0.378, y = 300, xend = 0.378, yend = 772.09), color = "#000000", linetype = "dashed") +
  geom_segment(aes(x = -.1, y = 772.09, xend = 0.38, yend = 772.09), color = "#000000", linetype = "dashed")


###################################################################################################
### --- Task 3: Plot with MCL --- #################################################################
###################################################################################################


# Set graph window parameters
x_min <- 0
x_max <- 1
y_min <- 400
y_max <- 1400

# Create another beautiful graph, this time with the MCL curve
ggplot() +
  geom_abline(intercept = - alpha / beta, slope = 1 / beta, color = "#5500ff") + # Demand
  geom_abline(intercept = gamma_H - (delta_H*alpha / beta), slope = delta_H / beta, color = "#ff0000") + # AICH
  geom_abline(intercept = mu_H - (nu_H*alpha / beta), slope = nu_H / beta, color = "#ff02ccbe") + # MCH
  geom_abline(intercept = mu_L - (nu_L*alpha / beta), slope = nu_L / beta, color = "darkgreen") + # MCL
  geom_point(data = point_data, aes(x = quants, y = prices, size = num_consumers), color = "#0600ab80", stroke = 2) + # Demand data points
  geom_point(data = point_data, aes(x = quants, y = costs, size = num_consumers), color = "#b6000080", stroke = 2) + # Cost data points
  scale_size_continuous(range = c(2, 10)) + # Control size range of circles
  coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max)) + # Set window
  labs(title = "Plot of Data and Linear Functions", x = "Quantity", y = "Dollars") + # Add title and axis labels
  theme(aspect.ratio = 1/2, axis.ticks.length  = unit(.1, "cm")) + # Control aspect ratio and unit ticks
  scale_y_continuous(breaks = seq(y_min, y_max, by = 100)) + # Set y-axis number labels
  scale_x_continuous(breaks = seq(x_min, x_max, by = 0.1)) + # Set x-axis number labels
  geom_text(aes(x = 0.05, y = 925, label = "Demand"), color = "#5500ff") +
  geom_text(aes(x = 0.8, y = 725, label = "AIC_H"), color = "#ff0000") +
  geom_text(aes(x = 0.92, y = 460, label = "MC_H"), color = "#ff02ccbe") +
  geom_text(aes(x = 0.45, y = 1275, label = "MC_L"), color = "darkgreen") +
  geom_polygon(data = dwl_coords, aes(x = x, y = y), fill = "#d4ff0080") +
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 1600), color = "black", linetype = "solid") +
  geom_segment(aes(x = 0.378, y = 300, xend = 0.378, yend = 772.09), color = "#000000", linetype = "dashed") +
  geom_segment(aes(x = -.1, y = 772.09, xend = 0.38, yend = 772.09), color = "#000000", linetype = "dashed")


###########################################################################################################
### --- Task 4: Plot with MCS --- #########################################################################
###########################################################################################################


# Set graph window parameters
x_min <- 0
x_max <- .9
y_min <- 550
y_max <- 925

# Set coordinates for DWL area to graph
dwl_coords <- data.frame(
  x = c(2.176, .378, .378), 
  y = c(131.68, 758, 772)
)

# Create the final beautiful graph, this time with the MCS curve
ggplot() +
  geom_abline(intercept = - alpha / beta, slope = 1 / beta, color = "#5500ff") + # Demand
  geom_abline(intercept = gamma_H - (delta_H*alpha / beta), slope = delta_H / beta, color = "#ff0000") + # AICH
  geom_abline(intercept = mu_H - (nu_H*alpha / beta), slope = nu_H / beta, color = "#ff02ccbe") + # MCH
  geom_abline(intercept = mu_S - (nu_S*alpha / beta), slope = nu_S / beta, color = "darkorange") + # MCS
  geom_point(data = point_data, aes(x = quants, y = prices, size = num_consumers), color = "#0600ab80", stroke = 2) + # Demand data points
  geom_point(data = point_data, aes(x = quants, y = costs, size = num_consumers), color = "#b6000080", stroke = 2) + # Cost data points
  scale_size_continuous(range = c(2, 10)) + # Control size range of circles
  coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max)) + # Set window
  labs(title = "Plot of Data and Linear Functions", x = "Quantity", y = "Dollars") + # Add title and axis labels
  theme(aspect.ratio = 1/2, axis.ticks.length  = unit(.1, "cm")) + # Control aspect ratio and unit ticks
  scale_y_continuous(breaks = seq(y_min, y_max, by = 50)) + # Set y-axis number labels
  scale_x_continuous(breaks = seq(x_min, x_max, by = 0.1)) + # Set x-axis number labels
  geom_text(aes(x = 0.05, y = 913, label = "Demand"), color = "#5500ff") +
  geom_text(aes(x = 0.8, y = 710, label = "AIC_H"), color = "#ff0000") +
  geom_text(aes(x = 0.6, y = 590, label = "MC_H"), color = "#ff02ccbe") +
  geom_text(aes(x = 0.48, y = 700, label = "MC_S"), color = "darkorange") +
  geom_polygon(data = dwl_coords, aes(x = x, y = y), fill = "#d4ff0080") +
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 1600), color = "black", linetype = "solid") +
  geom_segment(aes(x = 0.378, y = 300, xend = 0.378, yend = 772.09), color = "#000000", linetype = "dashed") +
  geom_segment(aes(x = -.1, y = 772.09, xend = 0.38, yend = 772.09), color = "#000000", linetype = "dashed")

