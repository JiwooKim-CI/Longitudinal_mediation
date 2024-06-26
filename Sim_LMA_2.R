# the longitudinal mediation model
library(tictoc)
library(plyr)
library(doParallel)
cl<-1
registerDoParallel(cl)
tic()
library(dagitty) # to draw DAGs
library(ggplot2)
library(ggthemes)
library(patchwork)
library(tikzDevice)
library(matrixcalc)
library(data.table)
library(dplyr)
library(tidyr)
library(jtools)
setDTthreads(percent = 90)
#### analytic bias from u1, u2, and u3 (delta_1, delta_2 and delta_3)
#### 5 grid points will only run on supercomputer
#### lower number of gridpoints yield similar results
ng <- seq(-.3, .3, length = 3)

p_m2_x <- ng # a path-fixed
p_y2_x <- ng 
p_y2_m2 <- ng 
p_m2_m1 <- ng 
p_y2_y1 <- ng 
p_y2_u2 <- ng 
p_m2_u2 <- ng 

p_u2_u1 <- ng
p_m1_u1 <- p_u2_u1*p_m2_u2 
p_y1_u1 <- p_u2_u1*p_y2_u2 
p_y1_m1 <- ng 
p_gm_m2 <- 1
p_gm_m1 <- -1
p_gy_y2 <- 1
p_gy_y1 <- -1


#using data.table here, in the hope that it is faster
#it is called grid1, because it is the 1st grid that we created
grid1 <- CJ(p_m2_x, p_y2_x, p_y2_m2, p_m2_m1, p_y2_y1, p_y2_u2, p_m2_u2, p_y1_u1,
            p_m1_u1, p_u2_u1, p_y1_m1, p_gm_m2, p_gm_m1, p_gy_y2, p_gy_y1)

grid1[,`:=`(eps_m1 = sqrt(1-p_m1_u1^2),
            eps_y1 = sqrt(1-p_y1_u1^2-
                            p_y1_m1^2-
                            2*p_y1_u1*p_m1_u1*p_y1_m1),
            eps_m2 = sqrt(1-p_m2_m1^2-
                            p_m2_u2^2-
                            p_m2_x^2-
                            2*p_m2_m1*p_m2_u2*p_m1_u1*p_u2_u1),
            eps_y2 = sqrt(1-p_y2_y1^2-
                            p_y2_m2^2-
                            p_y2_u2^2-
                            p_y2_x^2-
                            2*p_y2_x*p_m2_x*p_y2_x-
                            2*p_y2_m2*p_m2_m1*(p_y1_m1*p_y2_y1+p_m1_u1*p_y2_u2*p_u2_u1+p_m1_u1*p_y1_u1*p_y2_y1)-
                            2*p_y2_u2*p_y1_u1*p_u2_u1*p_y2_y1))]

grid1 <- grid1[eps_m1 >= 0 & eps_m2 >= 0 & eps_y1 >= 0 & eps_y2 >= 0]

# computing estimators 
grid1[,`:=`(rm2y2 = p_y2_m2 + p_m2_u2*p_y2_u2 + p_m2_m1*p_y1_m1*p_y2_y1 + 
              p_m2_m1*p_m1_u1*p_y1_u1*p_y2_y1 + p_m2_m1*p_m1_u1*p_m2_m1*p_y2_u2 + 
              p_m2_u2*p_u2_u1*p_y1_u1*p_y2_y1,
            rm1m2 = p_m2_m1 + p_m1_u1*p_u2_u1*p_m2_u2,
            rm1y2 = p_m2_m1*p_y2_m2 + p_m1_u1*p_u2_u1*p_y2_u2 + 
              p_y1_m1*p_y2_y1 + p_m1_u1*p_y1_u1*p_y2_y1,
            ry1m2 = p_m2_m1*p_y1_m1 + p_m2_m1*p_m1_u1*p_y1_u1,
            ry1y2 = p_y2_y1 + p_y1_m1*p_m2_m1*p_y2_m2 + p_y1_u1*p_m1_u1*p_m2_m1*p_y2_m2 + 
              p_y1_u1*p_u2_u1*p_y2_u2 + p_y1_u1*p_u2_u1*p_m2_u2*p_y2_m2,
            rm1y1 = p_y1_m1 + p_m1_u1*p_y1_u1)]

grid1[,`:=`(ry1m2.m1 = (ry1m2 - rm1y1*rm1m2)/sqrt((1-rm1y1^2)*(1-rm1m2^2)),
            ry1y2.m1 = (ry1y2 - rm1y1*rm1y2)/sqrt((1-rm1y1^2)*(1-rm2y2^2)),
            ry2m2.m1 = (rm2y2 - rm1m2*rm1y2)/sqrt((1-rm1m2^2)*(1-rm1y2^2)),
            sd_m2.m1 = sqrt(1 - rm1m2^2),
            sd_y2.m1 = sqrt(1 - rm1y2^2),
            sd_y1.m1 = sqrt(1 - rm1y1^2))]

grid1[,`:=`(change = p_m2_x*(p_gm_m1*p_m2_m1*p_y2_m2*p_gy_y2 +
                               p_gm_m1*p_y1_m1*p_y2_y1*p_gy_y2 +
                               p_gm_m1*p_y1_m1*p_y2_y1*p_gy_y1 +
                               p_gm_m1*p_m1_u1*p_u2_u1*p_m2_u2*p_y2_m2*p_gy_y2 +
                               p_gm_m1*p_m1_u1*p_u2_u1*p_y2_u2*p_gy_y2  +
                               p_gm_m1*p_m1_u1*p_y1_u1*p_gy_y1 +
                               p_gm_m1*p_m1_u1*p_y1_u1*p_y2_y1*p_gy_y2 +
                               p_gm_m2*p_y2_m2*p_gy_y2 +
                               p_gm_m2*p_y2_m2*p_m2_m1*p_y1_m1*p_gy_y1 +
                               p_gm_m2*p_y2_m2*p_m2_m1*p_y1_m1*p_y2_y1*p_gy_y2 +
                               p_gm_m2*p_y2_m2*p_m2_m1*p_m1_u1*p_u2_u1*p_y2_u2*p_gy_y2 +
                               p_gm_m2*p_y2_m2*p_m2_m1*p_m1_u1*p_y1_u1*p_gy_y1 +
                               p_gm_m2*p_y2_m2*p_m2_m1*p_m1_u1*p_y1_u1*p_y2_y1*p_gy_y2 +
                               p_gm_m2*p_m2_u2*p_y2_u2*p_gy_y2 +
                               p_gm_m2*p_m2_u2*p_u2_u1*p_m1_u1*p_y1_m1*p_gy_y1 +
                               p_gm_m2*p_m2_u2*p_u2_u1*p_m1_u1*p_y1_m1*p_y2_y1*p_gy_y2+
                               p_gm_m2*p_m2_u2*p_u2_u1*p_y1_u1*p_gy_y1+
                               p_gm_m2*p_m2_u2*p_u2_u1*p_y1_u1*p_y2_y1*p_gy_y2)/(2-2*(p_m2_m1^2+p_m1_u1*p_u2_u1*p_m2_u2)),
            change_2 = p_m2_x*2*(p_gm_m1*p_m2_m1*p_y2_m2*p_gy_y2 +
                                 p_gm_m1*p_y1_m1*p_y2_y1*p_gy_y2 +
                                 p_gm_m1*p_y1_m1*p_y2_y1*p_gy_y1 +
                                 p_gm_m1*p_m1_u1*p_u2_u1*p_m2_u2*p_y2_m2*p_gy_y2 +
                                 p_gm_m1*p_m1_u1*p_u2_u1*p_y2_u2*p_gy_y2  +
                                 p_gm_m1*p_m1_u1*p_y1_u1*p_gy_y1 +
                                 p_gm_m1*p_m1_u1*p_y1_u1*p_y2_y1*p_gy_y2 +
                                 p_gm_m2*p_y2_m2*p_gy_y2 +
                                 p_gm_m2*p_y2_m2*p_m2_m1*p_y1_m1*p_gy_y1 +
                                 p_gm_m2*p_y2_m2*p_m2_m1*p_y1_m1*p_y2_y1*p_gy_y2 +
                                 p_gm_m2*p_y2_m2*p_m2_m1*p_m1_u1*p_u2_u1*p_y2_u2*p_gy_y2 +
                                 p_gm_m2*p_y2_m2*p_m2_m1*p_m1_u1*p_y1_u1*p_gy_y1 +
                                 p_gm_m2*p_y2_m2*p_m2_m1*p_m1_u1*p_y1_u1*p_y2_y1*p_gy_y2 +
                                 p_gm_m2*p_m2_u2*p_y2_u2*p_gy_y2 +
                                 p_gm_m2*p_m2_u2*p_u2_u1*p_m1_u1*p_y1_m1*p_gy_y1 +
                                 p_gm_m2*p_m2_u2*p_u2_u1*p_m1_u1*p_y1_m1*p_y2_y1*p_gy_y2+
                                 p_gm_m2*p_m2_u2*p_u2_u1*p_y1_u1*p_gy_y1+
                                 p_gm_m2*p_m2_u2*p_u2_u1*p_y1_u1*p_y2_y1*p_gy_y2)/(2-2*(p_m2_m1^2+p_m1_u1*p_u2_u1*p_m2_u2)),
            change_1 = p_m2_x*(p_y2_m2*p_gy_y2 +
                                 p_m2_m1*p_y1_m1*p_gy_y1 +
                                 p_m2_m1*p_y1_m1*p_y2_y1*p_gy_y2 +
                                 p_m2_m1*p_m1_u1*p_u2_u1*p_y2_u2*p_gy_y2 +
                                 p_m2_m1*p_m1_u1*p_y1_u1*p_gy_y1 +
                                 p_m2_m1*p_m1_u1*p_y1_u1*p_y2_y1*p_gy_y2 +
                                 p_m2_u2*p_y2_u2*p_gy_y2 +
                                 p_m2_u2*p_u2_u1*p_m1_u1*p_y1_m1*p_gy_y1 +
                                 p_m2_u2*p_u2_u1*p_m1_u1*p_y1_m1*p_y2_y1*p_gy_y2 +
                                 p_m2_u2*p_u2_u1*p_y1_u1*p_gy_y1 +
                                 p_m2_u2*p_u2_u1*p_y1_u1*p_y2_y1*p_gy_y2),
            ancova = p_m2_x*(((ry2m2.m1*sd_y2.m1/sd_m2.m1)-(ry1m2.m1*sd_m2.m1/sd_y1.m1)*(ry1y2.m1*sd_y2.m1/sd_y1.m1))/
                               (1-(ry1m2.m1*sd_m2.m1/sd_y1.m1)^2)),
            truth = p_m2_x*p_y2_m2,
            naive = p_m2_x*(p_y2_m2+
                             p_m2_m1*p_y1_m1*p_y2_y1+
                             p_m2_m1*p_m1_u1*p_y1_u1*p_y2_y1+
                             p_m2_m1*p_m1_u1*p_u2_u1*p_y2_u2+
                             p_m2_u2*p_u2_u1*p_m1_u1*p_y1_m1*p_y2_y1+
                             p_m2_u2*p_u2_u1*p_y1_u1*p_y2_y1+
                             p_m2_u2*p_y2_u2))]

# computing bias
grid1[,`:=`(bias_change = change - truth,
            bias_change_2 = change_2 - truth,
            bias_ancova = ancova - truth,
            bias_naive = naive - truth,
            bias_change_1 = change_1 - truth)]

# compare the abs biases, both larger and equal, and strictly larger
grid1$bias_ind <- abs(grid1$bias_change) >= abs(grid1$bias_ancova)
table(grid1$bias_ind)

grid1$bias_ind_eq <- abs(grid1$bias_change) > abs(grid1$bias_ancova)
table(grid1$bias_ind_eq)

# in order to fix floating point errors we use the dplyr near function
# as a safer way to compare numbers
# So, we try to correct the issue using "near" function
eq <- near(grid1$bias_change, grid1$bias_ancova)
table(eq)

# for the comparison, |Bias^N| >= |Bias^F|
# return 1 if the two biases are near equal
grid1$bias_ind <- ifelse(eq == 1, 1, grid1$bias_ind)
table(grid1$bias_ind)

# for the comparison, |Bias^N| > |Bias^F|
# return 0 if the two biases are near equal
grid1$bias_ind_eq <- ifelse(eq == 1, 0, grid1$bias_ind_eq)
table(grid1$bias_ind_eq)

# compare the abs biases, both larger and equal, and strictly larger
grid1$bias_ind_2 <- abs(grid1$bias_change_2) >= abs(grid1$bias_ancova)
table(grid1$bias_ind_2)

grid1$bias_ind_eq_2 <- abs(grid1$bias_change_2) > abs(grid1$bias_ancova)
table(grid1$bias_ind_eq_2)

# in order to fix floating point errors we use the dplyr near function
# as a safer way to compare numbers
# So, we try to correct the issue using "near" function
eq_2 <- near(grid1$bias_change_2, grid1$bias_ancova)
table(eq_2)

# for the comparison, |Bias^N| >= |Bias^F|
# return 1 if the two biases are near equal
grid1$bias_ind_2 <- ifelse(eq_2 == 1, 1, grid1$bias_ind_2)
table(grid1$bias_ind_2)

# for the comparison, |Bias^N| > |Bias^F|
# return 0 if the two biases are near equal
grid1$bias_ind_eq_2 <- ifelse(eq == 1, 0, grid1$bias_ind_eq_2)
table(grid1$bias_ind_eq_2)

# RESULTS---------------------------------------------------------------------------------------------
# Final number of |Bias_naive| >= |Bias_fd|
# no constraint
(f.tab <- with(grid1, table(bias_ind)))
prop.table(f.tab)

# p_m1_u == p_m2_u
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_m2_u2*p_u2_u1,], 
               table(bias_ind)))
prop.table(f.tab)

# p_y1_u == p_y2_u
(f.tab <- with(grid1[grid1$p_y1_u1 == grid1$p_y2_u2*p_u2_u1,], 
               table(bias_ind)))
prop.table(f.tab)

# p_y1_m1 == 0
(f.tab <- with(grid1[grid1$p_y1_m1 == 0,], 
               table(bias_ind)))
prop.table(f.tab)

# p_m2_m1 == 0
(f.tab <- with(grid1[grid1$p_m2_m1 == 0,], 
               table(bias_ind)))
prop.table(f.tab)

# p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_y2_y1 == 0,], 
               table(bias_ind)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_u == p_y2_u
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_m2_u2*p_u2_u1 &
                       grid1$p_y1_u1 == grid1$p_y2_u2*p_u2_u1 ,], 
               table(bias_ind)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_m1==0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_m2_u2*p_u2_u1 &
                       grid1$p_y1_m1 == 0,], 
               table(bias_ind)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_m2_m1==0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_m2_u2*p_u2_u1 &
                       grid1$p_m2_m1 == 0,], 
               table(bias_ind)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_m2_u2*p_u2_u1 &
                       grid1$p_y2_y1 == 0,], 
               table(bias_ind)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_y1_m1 == 0
(f.tab <- with(grid1[grid1$p_y1_u1 == grid1$p_y2_u2*p_u2_u1 &
                       grid1$p_y1_m1 == 0,], 
               table(bias_ind)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_m2_m1 == 0
(f.tab <- with(grid1[grid1$p_y1_u1 == grid1$p_y2_u2*p_u2_u1 &
                       grid1$p_m2_m1 == 0,], 
               table(bias_ind)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_y1_u1 == grid1$p_y2_u2*p_u2_u1 &
                       grid1$p_y2_y1 == 0,], 
               table(bias_ind)))
prop.table(f.tab)

# p_y1_m1 == 0 & p_m2_m1 == 0
(f.tab <- with(grid1[grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0,], 
               table(bias_ind)))
prop.table(f.tab)

# p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(bias_ind)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_u == p_y2_u & p_y1_m1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_m2_u2*p_u2_u1 &
                       grid1$p_y1_u1 == grid1$p_y2_u2*p_u2_u1 &
                       grid1$p_y1_m1 == 0,], 
               table(bias_ind)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_u == p_y2_u & p_m2_m1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_m2_u2*p_u2_u1 &
                       grid1$p_y1_u1 == grid1$p_y2_u2*p_u2_u1 &
                       grid1$p_m2_m1 == 0,], 
               table(bias_ind)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_u == p_y2_u & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_m2_u2*p_u2_u1 &
                       grid1$p_y1_u1 == grid1$p_y2_u2*p_u2_u1 &
                       grid1$p_y2_y1 == 0,], 
               table(bias_ind)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_m1 == 0 & p_m2_m1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_m2_u2*p_u2_u1 &
                       grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0,], 
               table(bias_ind)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_m2_u2*p_u2_u1 &
                       grid1$p_y1_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(bias_ind)))
prop.table(f.tab)


# p_m1_u == p_m2_u & p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_m2_u2*p_u2_u1 &
                       grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(bias_ind)))
prop.table(f.tab)

# p_y1_m1 == 0 & p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(bias_ind)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_y1_m1 == 0 & p_m2_m1 == 0
(f.tab <- with(grid1[grid1$p_y1_u1 == grid1$p_y2_u2*p_u2_u1 &
                       grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0,], 
               table(bias_ind)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_y1_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_y1_u1 == grid1$p_y2_u2*p_u2_u1 &
                       grid1$p_y1_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(bias_ind)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_y1_u1 == grid1$p_y2_u2*p_u2_u1 &
                       grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(bias_ind)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_m1 == 0 & p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_m2_u2*p_u2_u1 &
                       grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(bias_ind)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_y1_m1 == 0 & p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_y1_u1 == grid1$p_y2_u2*p_u2_u1 &
                       grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(bias_ind)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_m1_u == p_m2_u & p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_m2_u2*p_u2_u1 &
                       grid1$p_y1_u1 == grid1$p_y2_u2*p_u2_u1 &
                       grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(bias_ind)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_m1_u == p_m2_u & p_y1_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_m2_u2*p_u2_u1 &
                       grid1$p_y1_u1 == grid1$p_y2_u2*p_u2_u1 &
                       grid1$p_y1_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(bias_ind)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_m1_u == p_m2_u & p_y1_m1 == 0 & p_m2_m1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_m2_u2*p_u2_u1 &
                       grid1$p_y1_u1 == grid1$p_y2_u2*p_u2_u1 &
                       grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0,], 
               table(bias_ind)))
prop.table(f.tab)

# common trend assumption
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_m2_u2*p_u2_u1 &
                       grid1$p_y1_u1 == grid1$p_y2_u2*p_u2_u1 &
                       grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(bias_ind_2)))
prop.table(f.tab)


####### graphs
#overlaid densities of both biases
graph1 <- pivot_longer(grid1, cols = c(bias_change, bias_naive, bias_ancova, bias_change_1),
                       names_to = "bias_type", values_to = "bias_result")
p11 <- ggplot(graph1,aes(x=(bias_result),y=..count../max(count), color = bias_type)) + 
  geom_density(bw=.13,trim=FALSE) +
  labs(x="Raw bias", y="Percent") + 
  xlim(-1.7,1.7)+
  theme_apa() + 
  theme(legend.position = c(.2,.8))

#tikz("biasdensity.tex",width = 4,height = 3, standAlone = TRUE)
p11
#dev.off()

graph2 <- pivot_longer(grid1, cols = c(bias_change_2, bias_naive, bias_ancova, bias_change_1),
                       names_to = "bias_type", values_to = "bias_result")
p12 <- ggplot(graph2,aes(x=bias_result, y=..count../max(count), color = bias_type)) + 
  geom_density(bw=.13,trim=FALSE) + 
  labs(x="Raw bias", y="Percent") + 
  xlim(-1.7,1.7)+
  theme_apa() + theme(legend.position = c(.2,.8))
#tikz("biasdensity_2.tex",width = 4,height = 3, standAlone = TRUE)
p12
#dev.off()

#numerical values from table
#both strictly larger than 0 (bias_diff_ind)
#and larger or equal than 0 (bias_diff_ind_eq)
ptable1 <- prop.table(table(grid1$bias_ind))
ptable1b <- prop.table(table(grid1$bias_ind_eq))
ptable1_2 <- prop.table(table(grid1$bias_ind_2))
ptable1b_2 <- prop.table(table(grid1$bias_ind_eq_2))

grid1small1 <- grid1 |> filter(p_m1_u1 == p_m2_u2*p_u2_u1 &
                                 p_y1_u1 == p_y2_u2*p_u2_u1 & 
                                 p_y1_m1 == 0 &
                                 p_m2_m1 == 0 &
                                 p_y2_y1 == 0)


#constraining direct effect gamma = 0, no direct effect violation
graph3 <- pivot_longer(grid1small1, cols = c(bias_change_2, bias_naive, bias_ancova, bias_change_1),
                       names_to = "bias_type", values_to = "bias_result")
p2 <- ggplot(graph3,aes(x=bias_result,y=..count../max(count), color = bias_type)) + 
  geom_histogram() + facet_grid(~bias_type)
  labs(x="Raw bias", y="Percent") + 
  xlim(-1.7,1.7) + 
  theme_apa() + theme(legend.position = c(.2,.8))
#tikz("biasdensitynogamma.tex",width = 4,height = 3, standAlone = TRUE)
p2
#dev.off()
table(grid1small1$bias_change_1)
ptable2 <- prop.table(table(grid1small1$bias_ind_2))
ptable2b <- prop.table(table(grid1small1$bias_ind_eq_2))


toc()


stopCluster(cl)