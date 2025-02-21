# Load required library
library(caracas)

## Example model
# -----------------------------------------------------------------------------
# Define Model Parameters
# -----------------------------------------------------------------------------

p_y2_m2 <- symbol("beta")
p_m2_u <- symbol("delta3")
p_y2_u <- symbol("delta4")

# Fixed path coefficients (set to zero)
p_m2_m1 <- 0
p_y1_m1 <- 0
p_y2_y1 <- 0
p_y2_m1 <- 0
p_m2_y1 <- 0

p_m1_u <- symbol("delta1")
p_y1_u <- symbol("delta4")
p_m2_x <- symbol("alpha")
p_y2_x <- symbol("gamma")

# -----------------------------------------------------------------------------
# Compute Intermediate Terms
# -----------------------------------------------------------------------------

rm2y2 <- p_y2_m2 + p_m2_u * p_y2_u + p_m2_x * p_y2_x
rm1m2 <- p_m1_u * p_m2_u
rm1y2 <- p_m1_u * p_y2_u + p_m1_u * p_m2_u * p_y2_m2
ry1m2 <- p_y1_u * p_m2_u
ry1y2 <- p_y1_u * p_y2_u + p_y1_u * p_m2_u * p_y2_m2
rm1y1 <- p_m1_u * p_y1_u

# Fixed correlations (set to zero)
rm1x <- 0
ry1x <- 0
rm2x <- p_m2_x
ry2x <- p_m2_x * p_y2_m2 + p_y2_x

# -----------------------------------------------------------------------------
# Compute Conditional Correlations and Standard Deviations
# -----------------------------------------------------------------------------

ry1m2_m1 <- (ry1m2 - rm1y1 * rm1m2) / sqrt((1 - rm1y1^2) * (1 - rm1m2^2))
ry1y2_m1 <- (ry1y2 - rm1y1 * rm1y2) / sqrt((1 - rm1y1^2) * (1 - rm1y2^2))
rm2y2_m1 <- (rm2y2 - rm1m2 * rm1y2) / sqrt((1 - rm1m2^2) * (1 - rm1y2^2))
ry1x_m1  <- (ry1x - rm1y1 * rm1x) / sqrt((1 - rm1y1^2) * (1 - rm1x^2))
rxy2_m1  <- (ry2x - rm1x * rm1y2) / sqrt((1 - rm1x^2) * (1 - rm1y2^2))
rxm2_m1  <- (rm2x - rm1x * rm1m2) / sqrt((1 - rm1x^2) * (1 - rm1m2^2))

sd_m2_m1 <- sqrt(1 - rm1m2^2)
sd_y2_m1 <- sqrt(1 - rm1y2^2)
sd_y1_m1 <- sqrt(1 - rm1y1^2)

sd_m2_m1y1 <- sqrt(1 - rm1y1^2 - (ry1m2^2 + rm1m2^2 - 2 * ry1m2 * rm1m2 * rm1y1))
sd_y2_m1y1 <- sqrt(1 - rm1y1^2 - (ry1y2^2 + rm1y2^2 - 2 * ry1y2 * rm1y2 * rm1y1))

rm2y2_m1y1 <- (rm2y2_m1 - ry1m2_m1 * ry1y2_m1) / sqrt((1 - ry1m2_m1^2) * (1 - ry1y2_m1^2))
rxy2_m1y1  <- (rxy2_m1 - ry1x_m1 * ry1y2_m1) / sqrt((1 - ry1x_m1^2) * (1 - ry1y2_m1^2))
rxm2_m1y1  <- (rxm2_m1 - ry1x_m1 * ry1m2_m1) / sqrt((1 - ry1x_m1^2) * (1 - ry1m2_m1^2))

# -----------------------------------------------------------------------------
# Compute Final Expression
# -----------------------------------------------------------------------------

expr <- (((rm2y2_m1y1 - rxm2_m1y1 * rxy2_m1y1) /
          (1 - rxm2_m1y1^2)) * sd_y2_m1y1 / sd_m2_m1y1)

# Output the expression
expr

# Convert to LaTeX format
latex_expr <- as.character(expr, format = "latex")
cat(latex_expr)


## Set delta_3 = 0
# -----------------------------------------------------------------------------
# Define Model Parameters
# -----------------------------------------------------------------------------

p_y2_m2 <- symbol("beta")
p_m2_u <- 0
p_y2_u <- symbol("delta4")

# Fixed path coefficients (set to zero)
p_m2_m1 <- 0
p_y1_m1 <- 0
p_y2_y1 <- 0
p_y2_m1 <- 0
p_m2_y1 <- 0

p_m1_u <- symbol("delta1")
p_y1_u <- symbol("delta4")
p_m2_x <- symbol("alpha")
p_y2_x <- symbol("gamma")

# -----------------------------------------------------------------------------
# Compute Intermediate Terms
# -----------------------------------------------------------------------------

rm2y2 <- p_y2_m2 + p_m2_u * p_y2_u + p_m2_x * p_y2_x
rm1m2 <- p_m1_u * p_m2_u
rm1y2 <- p_m1_u * p_y2_u + p_m1_u * p_m2_u * p_y2_m2
ry1m2 <- p_y1_u * p_m2_u
ry1y2 <- p_y1_u * p_y2_u + p_y1_u * p_m2_u * p_y2_m2
rm1y1 <- p_m1_u * p_y1_u

# Fixed correlations (set to zero)
rm1x <- 0
ry1x <- 0
rm2x <- p_m2_x
ry2x <- p_m2_x * p_y2_m2 + p_y2_x

# -----------------------------------------------------------------------------
# Compute Conditional Correlations and Standard Deviations
# -----------------------------------------------------------------------------

ry1m2_m1 <- (ry1m2 - rm1y1 * rm1m2) / sqrt((1 - rm1y1^2) * (1 - rm1m2^2))
ry1y2_m1 <- (ry1y2 - rm1y1 * rm1y2) / sqrt((1 - rm1y1^2) * (1 - rm1y2^2))
rm2y2_m1 <- (rm2y2 - rm1m2 * rm1y2) / sqrt((1 - rm1m2^2) * (1 - rm1y2^2))
ry1x_m1  <- (ry1x - rm1y1 * rm1x) / sqrt((1 - rm1y1^2) * (1 - rm1x^2))
rxy2_m1  <- (ry2x - rm1x * rm1y2) / sqrt((1 - rm1x^2) * (1 - rm1y2^2))
rxm2_m1  <- (rm2x - rm1x * rm1m2) / sqrt((1 - rm1x^2) * (1 - rm1m2^2))

sd_m2_m1 <- sqrt(1 - rm1m2^2)
sd_y2_m1 <- sqrt(1 - rm1y2^2)
sd_y1_m1 <- sqrt(1 - rm1y1^2)

sd_m2_m1y1 <- sqrt(1 - rm1y1^2 - (ry1m2^2 + rm1m2^2 - 2 * ry1m2 * rm1m2 * rm1y1))
sd_y2_m1y1 <- sqrt(1 - rm1y1^2 - (ry1y2^2 + rm1y2^2 - 2 * ry1y2 * rm1y2 * rm1y1))

rm2y2_m1y1 <- (rm2y2_m1 - ry1m2_m1 * ry1y2_m1) / sqrt((1 - ry1m2_m1^2) * (1 - ry1y2_m1^2))
rxy2_m1y1  <- (rxy2_m1 - ry1x_m1 * ry1y2_m1) / sqrt((1 - ry1x_m1^2) * (1 - ry1y2_m1^2))
rxm2_m1y1  <- (rxm2_m1 - ry1x_m1 * ry1m2_m1) / sqrt((1 - ry1x_m1^2) * (1 - ry1m2_m1^2))

# -----------------------------------------------------------------------------
# Compute Final Expression
# -----------------------------------------------------------------------------

expr2 <- (((rm2y2_m1y1 - rxm2_m1y1 * rxy2_m1y1) /
          (1 - rxm2_m1y1^2)) * sd_y2_m1y1 / sd_m2_m1y1)

# Output the expression
expr2

# Convert to LaTeX format
latex_expr2 <- as.character(expr2, format = "latex")
cat(latex_expr2)

## Set delta_4 = 0
# -----------------------------------------------------------------------------
# Define Model Parameters
# -----------------------------------------------------------------------------

p_y2_m2 <- symbol("beta")
p_m2_u <- symbol("delta3")
p_y2_u <- 0

# Fixed path coefficients (set to zero)
p_m2_m1 <- 0
p_y1_m1 <- 0
p_y2_y1 <- 0
p_y2_m1 <- 0
p_m2_y1 <- 0

p_m1_u <- symbol("delta1")
p_y1_u <- 0
p_m2_x <- symbol("alpha")
p_y2_x <- symbol("gamma")

# -----------------------------------------------------------------------------
# Compute Intermediate Terms
# -----------------------------------------------------------------------------

rm2y2 <- p_y2_m2 + p_m2_u * p_y2_u + p_m2_x * p_y2_x
rm1m2 <- p_m1_u * p_m2_u
rm1y2 <- p_m1_u * p_y2_u + p_m1_u * p_m2_u * p_y2_m2
ry1m2 <- p_y1_u * p_m2_u
ry1y2 <- p_y1_u * p_y2_u + p_y1_u * p_m2_u * p_y2_m2
rm1y1 <- p_m1_u * p_y1_u

# Fixed correlations (set to zero)
rm1x <- 0
ry1x <- 0
rm2x <- p_m2_x
ry2x <- p_m2_x * p_y2_m2 + p_y2_x

# -----------------------------------------------------------------------------
# Compute Conditional Correlations and Standard Deviations
# -----------------------------------------------------------------------------

ry1m2_m1 <- (ry1m2 - rm1y1 * rm1m2) / sqrt((1 - rm1y1^2) * (1 - rm1m2^2))
ry1y2_m1 <- (ry1y2 - rm1y1 * rm1y2) / sqrt((1 - rm1y1^2) * (1 - rm1y2^2))
rm2y2_m1 <- (rm2y2 - rm1m2 * rm1y2) / sqrt((1 - rm1m2^2) * (1 - rm1y2^2))
ry1x_m1  <- (ry1x - rm1y1 * rm1x) / sqrt((1 - rm1y1^2) * (1 - rm1x^2))
rxy2_m1  <- (ry2x - rm1x * rm1y2) / sqrt((1 - rm1x^2) * (1 - rm1y2^2))
rxm2_m1  <- (rm2x - rm1x * rm1m2) / sqrt((1 - rm1x^2) * (1 - rm1m2^2))

sd_m2_m1 <- sqrt(1 - rm1m2^2)
sd_y2_m1 <- sqrt(1 - rm1y2^2)
sd_y1_m1 <- sqrt(1 - rm1y1^2)

sd_m2_m1y1 <- sqrt(1 - rm1y1^2 - (ry1m2^2 + rm1m2^2 - 2 * ry1m2 * rm1m2 * rm1y1))
sd_y2_m1y1 <- sqrt(1 - rm1y1^2 - (ry1y2^2 + rm1y2^2 - 2 * ry1y2 * rm1y2 * rm1y1))

rm2y2_m1y1 <- (rm2y2_m1 - ry1m2_m1 * ry1y2_m1) / sqrt((1 - ry1m2_m1^2) * (1 - ry1y2_m1^2))
rxy2_m1y1  <- (rxy2_m1 - ry1x_m1 * ry1y2_m1) / sqrt((1 - ry1x_m1^2) * (1 - ry1y2_m1^2))
rxm2_m1y1  <- (rxm2_m1 - ry1x_m1 * ry1m2_m1) / sqrt((1 - ry1x_m1^2) * (1 - ry1m2_m1^2))

# -----------------------------------------------------------------------------
# Compute Final Expression
# -----------------------------------------------------------------------------

expr3 <- (((rm2y2_m1y1 - rxm2_m1y1 * rxy2_m1y1) /
          (1 - rxm2_m1y1^2)) * sd_y2_m1y1 / sd_m2_m1y1)

# Output the expression
expr3

# Convert to LaTeX format
latex_expr3 <- as.character(expr3, format = "latex")
cat(latex_expr3)
