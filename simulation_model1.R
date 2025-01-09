# the longitudinal mediation model
# This code is written assuming no time-varying confounders assumption

# load packages
library(tictoc)
library(plyr)
library(doParallel)
cl<-49
registerDoParallel(cl)
library(dagitty) 
library(ggplot2)
library(ggthemes)
library(patchwork)
library(tikzDevice)
library(matrixcalc)
library(data.table)
library(dplyr)
library(tidyr)
library(jtools)
library(profvis)
library(cgwtools)
library(forestmangr)
library(dendroTools)


#### analytic bias from u1
#### 9 grid points will only run on supercomputer
#### lower number of gridpoints yield similar results
ng <- list()
ng <- append(.1, ng)
ng <- append(sqrt(ng[[1]]),ng)
ng <- append(sqrt(ng[[1]]),ng)
ng <- append(sqrt(ng[[1]]),ng)
ng <- append(0,ng)
ng <- append(-1*ng[[2]],ng)
ng <- append(-1*ng[[4]],ng)
ng <- append(-1*ng[[6]],ng)
ng <- append(-1*ng[[8]],ng)
ng <- as.numeric(ng)

p_m2_x <- .2 # a path-fixed
p_y2_x <- ng 
p_y2_m2 <- ng 
p_m2_m1 <- ng 
p_y2_y1 <- ng 
p_y2_u <- ng 
p_m2_u <- ng 
p_y1_u <- ng 
p_m1_u <- ng 
p_y1_m1 <- ng 
p_gm_m2 <- 1
p_gm_m1 <- -1
p_gy_y2 <- 1
p_gy_y1 <- -1


#using data.table here, in the hope that it is faster
#it is called grid1, because it is the 1st grid that we created
grid1 <- CJ(p_m2_x, p_y2_x, p_y2_m2, p_m2_m1, p_y2_y1, p_y2_u, p_m2_u, p_y1_u,
            p_m1_u, p_y1_m1, p_gm_m2, p_gm_m1, p_gy_y2, p_gy_y1)

grid1 <- as_tibble(grid1)
setDT(grid1)


# Convert to a tibble for better display if desired

grid1[,`:=`(eps_m1 = sqrt(1-p_m1_u^2),
            eps_y1 = sqrt(1-p_y1_u^2-
                            p_y1_m1^2-
                            2*p_y1_u*p_m1_u*p_y1_m1),
            eps_m2 = sqrt(1-p_m2_m1^2-
                            p_m2_u^2-
                            p_m2_x^2-
                            2*p_m2_m1*p_m2_u*p_m1_u),
            eps_y2 = sqrt(1-p_y2_y1^2-
                            p_y2_m2^2-
                            p_y2_u^2-
                            p_y2_x^2-
                            2*p_y2_x*p_m2_x*p_y2_x-
                            2*p_y2_m2*p_m2_m1*(p_y1_m1*p_y2_y1+p_m1_u*p_y2_u+p_m1_u*p_y1_u*p_y2_y1)-
                            2*p_y2_u*p_y1_u*p_y2_y1))]
# use na.omit code to eliminate imaginary number data
grid1 <- na.omit(grid1)

# computing path correlation
grid1[,`:=`(rm2y2 = p_y2_m2 + 
              p_m2_u*p_y2_u + 
              p_m2_m1*p_y1_m1*p_y2_y1 + 
              p_m2_m1*p_m1_u*p_y1_u*p_y2_y1 + 
              p_m2_m1*p_m1_u*p_y2_u +
              p_m2_u*p_y1_u*p_y2_y1 +
              p_m2_u*p_m1_u*p_y1_m1*p_y2_y1 ,
            rm1m2 = p_m2_m1 + p_m1_u*p_m2_u,
            rm1y2 = p_m2_m1*p_y2_m2 + p_m1_u*p_y2_u + p_y1_m1*p_y2_y1 + p_m1_u*p_y1_u*p_y2_y1 +
              p_m1_u*p_m2_u*p_y2_m2,
            ry1m2 = p_m2_m1*p_y1_m1 + p_m2_m1*p_m1_u*p_y1_u +
              p_m2_u*p_m1_u*p_y1_m1 +
              p_m2_u*p_y1_u,
            ry1y2 = p_y2_y1 + 
              p_y1_m1*p_m2_m1*p_y2_m2 + 
              p_y1_u*p_m1_u*p_m2_m1*p_y2_m2 + 
              p_y1_u*p_y2_u +
              p_y1_u*p_m2_u*p_y2_m2 +
              p_y1_m1*p_m1_u*p_m2_u*p_y2_m2 +
              p_y1_m1*p_m1_u*p_y2_u,
            rm1y1 = p_y1_m1 + p_m1_u*p_y1_u)]

grid1[,`:=`(ry1m2.m1 = (ry1m2-rm1y1*rm1m2)/sqrt((1-rm1y1^2)*(1-rm1m2^2)),
            ry1y2.m1 = (ry1y2-rm1y1*rm1y2)/sqrt((1-rm1y1^2)*(1-rm1y2^2)),
            rm2y2.m1 = (rm2y2-rm1m2*rm1y2)/sqrt((1-rm1m2^2)*(1-rm1y2^2)),
            sd_m2.m1 = sqrt(1-rm1m2^2),
            sd_y2.m1 = sqrt(1-rm1y2^2),
            sd_y1.m1 = sqrt(1-rm1y1^2))]

# use na.omit code to eliminate imaginary number data
grid1 <- na.omit(grid1)

# computing estimators
grid1[,`:=`(change = p_m2_x*(p_gm_m1*p_m2_m1*p_y2_m2*p_gy_y2 +
                               p_gm_m1*p_y1_m1*p_gy_y1+
                               p_gm_m1*p_y1_m1*p_y2_y1*p_gy_y2+
                               p_gm_m1*p_m1_u*p_m2_u*p_y2_m2*p_gy_y2+
                               p_gm_m1*p_m1_u*p_y1_u*p_gy_y1+
                               p_gm_m1*p_m1_u*p_y1_u*p_y2_y1*p_gy_y2+
                               p_gm_m1*p_m1_u*p_y2_u*p_gy_y2+
                               p_gm_m2*p_y2_m2*p_gy_y2+
                               p_gm_m2*p_m2_m1*p_y1_m1*p_gy_y1+
                               p_gm_m2*p_m2_m1*p_y1_m1*p_y2_y1*p_gy_y2+
                               p_gm_m2*p_m2_m1*p_m1_u*p_y1_u*p_gy_y1+
                               p_gm_m2*p_m2_m1*p_m1_u*p_y1_u*p_y2_y1*p_gy_y2+
                               p_gm_m2*p_m2_m1*p_m1_u*p_y2_u*p_gy_y2+
                               p_gm_m2*p_m2_u*p_m1_u*p_y1_m1*p_gy_y1+
                               p_gm_m2*p_m2_u*p_m1_u*p_y1_m1*p_y2_y1*p_gy_y2+
                               p_gm_m2*p_m2_u*p_y1_u*p_gy_y1+
                               p_gm_m2*p_m2_u*p_y1_u*p_y2_y1*p_gy_y2+
                               p_gm_m2*p_m2_u*p_y2_u*p_gy_y2)/(2-2*(p_m2_m1+p_m1_u*p_m2_u)),
            change_2 = p_m2_x*(p_gm_m1*p_m2_m1*p_y2_m2*p_gy_y2 +
                                 p_gm_m1*p_y1_m1*p_gy_y1+
                                 p_gm_m1*p_y1_m1*p_y2_y1*p_gy_y2+
                                 p_gm_m1*p_m1_u*p_m2_u*p_y2_m2*p_gy_y2+
                                 p_gm_m1*p_m1_u*p_y1_u*p_gy_y1+
                                 p_gm_m1*p_m1_u*p_y1_u*p_y2_y1*p_gy_y2+
                                 p_gm_m1*p_m1_u*p_y2_u*p_gy_y2+
                                 p_gm_m2*p_y2_m2*p_gy_y2+
                                 p_gm_m2*p_m2_m1*p_y1_m1*p_gy_y1+
                                 p_gm_m2*p_m2_m1*p_y1_m1*p_y2_y1*p_gy_y2+
                                 p_gm_m2*p_m2_m1*p_m1_u*p_y1_u*p_gy_y1+
                                 p_gm_m2*p_m2_m1*p_m1_u*p_y1_u*p_y2_y1*p_gy_y2+
                                 p_gm_m2*p_m2_m1*p_m1_u*p_y2_u*p_gy_y2+
                                 p_gm_m2*p_m2_u*p_m1_u*p_y1_m1*p_gy_y1+
                                 p_gm_m2*p_m2_u*p_m1_u*p_y1_m1*p_y2_y1*p_gy_y2+
                                 p_gm_m2*p_m2_u*p_y1_u*p_gy_y1+
                                 p_gm_m2*p_m2_u*p_y1_u*p_y2_y1*p_gy_y2+
                                 p_gm_m2*p_m2_u*p_y2_u*p_gy_y2)/(1-1*(p_m2_m1+p_m1_u*p_m2_u)),
            change_1 = p_m2_x*(p_y2_m2*p_gy_y2 +
                                 p_m2_m1*p_y1_m1*p_gy_y1 +
                                 p_m2_m1*p_y1_m1*p_y2_y1*p_gy_y2 +
                                 p_m2_m1*p_m1_u*p_y1_u*p_gy_y1 +
                                 p_m2_m1*p_m1_u*p_y1_u*p_y2_y1*p_gy_y2 +
                                 p_m2_m1*p_m1_u*p_y2_u*p_gy_y2 +
                                 p_m2_u*p_y1_u*p_gy_y1 +
                                 p_m2_u*p_y1_u*p_y2_y1*p_gy_y2 +
                                 p_m2_u*p_y2_u*p_gy_y2),
            ancova = p_m2_x*(((rm2y2.m1*sd_y2.m1/sd_m2.m1)-(ry1m2.m1*sd_y1.m1/sd_m2.m1)*(ry1y2.m1*sd_y2.m1/sd_y1.m1))/
                               (1-((ry1m2.m1*sd_y1.m1/sd_m2.m1)^2*sd_m2.m1^2/sd_y1.m1^2))),
            truth = p_m2_x*p_y2_m2,
            naive =p_m2_x*(p_y2_m2+
                             p_m2_m1*p_y1_m1*p_y2_y1+
                             p_m2_m1*p_m1_u*p_y1_u*p_y2_y1+
                             p_m2_m1*p_m1_u*p_y2_u+
                             p_m2_u*p_m1_u*p_y1_m1*p_y2_y1+
                             p_m2_u*p_y1_u*p_y2_y1+
                             p_m2_u*p_y2_u))]


# computing bias
grid1[,`:=`(bias_change = change - truth, ### 2 wave change score that people usually use
            bias_change_2 = change_2 - truth, ### 2 wave change score multiplied by 2
            bias_ancova = ancova - truth,
            bias_naive = naive - truth,
            bias_change_1 = change_1 - truth)] ### change score with one pre-test
grid1$diff <- abs(grid1$p_m2_u*grid1$p_y2_u)+abs(grid1$p_m1_u*grid1$p_y1_u)
grid1$diff2 <- abs(grid1$p_m1_u- grid1$p_m2_u - grid1$p_m1_u*grid1$p_m2_m1)+
                     abs(grid1$p_y1_u - grid1$p_y2_u - grid1$p_y1_u*grid1$p_y2_y1)
grid1$y2y1 <- grid1$p_y2_y1
grid1 <- grid1[,-(15:36)]
colnames(grid1)


grid1[15:19] <- round_df(grid1[15:19], 3)
grid1$Min_1 <- colnames(grid1[,15:19])[minCol(abs(grid1[,15:19]), ties.method = "first")]
grid1$Min_2 <- colnames(grid1[,15:19])[minCol(abs(grid1[,15:19]), ties.method = "last")]
grid1<- grid1 |> mutate(Min = ifelse(Min_1 != Min_2, "Tie", Min_1))


# RESULTS---------------------------------------------------------------------------------------------
# Final number of |Bias_naive| >= |Bias_fd|
# no constraint
(f.tab <- with(grid1, table(Min)))
prop.table(f.tab)

# p_m1_u == p_m2_u
(f.tab <- with(grid1[grid1$p_m1_u == grid1$p_m2_u+grid1$p_m1_u*grid1$p_m2_m1,], 
               table(Min)))
prop.table(f.tab)

# p_y1_u == p_y2_u
(f.tab <- with(grid1[grid1$p_y1_u == grid1$p_y2_u+ grid1$p_y1_u*grid1$p_y2_y1,], 
               table(Min)))
prop.table(f.tab)

# p_y1_m1 == 0
(f.tab <- with(grid1[grid1$p_y1_m1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_m2_m1 == 0
(f.tab <- with(grid1[grid1$p_m2_m1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_y2_y1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_u == p_y2_u
(f.tab <- with(grid1[grid1$p_m1_u == grid1$p_m2_u+grid1$p_m1_u*grid1$p_m2_m1 &
                       grid1$p_y1_u == grid1$p_y2_u + grid1$p_y1_u*grid1$p_y2_y1,], 
               table(Min)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_m1==0
(f.tab <- with(grid1[grid1$p_m1_u == grid1$p_m2_u+grid1$p_m1_u*grid1$p_m2_m1 &
                       grid1$p_y1_m1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_m2_m1==0
(f.tab <- with(grid1[grid1$p_m1_u == grid1$p_m2_u+grid1$p_m1_u*grid1$p_m2_m1 &
                       grid1$p_m2_m1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u == grid1$p_m2_u+grid1$p_m1_u*grid1$p_m2_m1 &
                       grid1$p_y2_y1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_y1_m1 == 0
(f.tab <- with(grid1[grid1$p_y1_u == grid1$p_y2_u+ grid1$p_y1_u*grid1$p_y2_y1 &
                       grid1$p_y1_m1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_m2_m1 == 0
(f.tab <- with(grid1[grid1$p_y1_u == grid1$p_y2_u+ grid1$p_y1_u*grid1$p_y2_y1 &
                       grid1$p_m2_m1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_y1_u == grid1$p_y2_u+ grid1$p_y1_u*grid1$p_y2_y1 &
                       grid1$p_y2_y1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_y1_m1 == 0 & p_m2_m1 == 0
(f.tab <- with(grid1[grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_u == p_y2_u & p_y1_m1 == 0
(f.tab <- with(grid1[grid1$p_m1_u == grid1$p_m2_u+grid1$p_m1_u*grid1$p_m2_m1 &
                       grid1$p_y1_u == grid1$p_y2_u+ grid1$p_y1_u*grid1$p_y2_y1 &
                       grid1$p_y1_m1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_u == p_y2_u & p_m2_m1 == 0
(f.tab <- with(grid1[grid1$p_m1_u == grid1$p_m2_u+grid1$p_m1_u*grid1$p_m2_m1 &
                       grid1$p_y1_u == grid1$p_y2_u+ grid1$p_y1_u*grid1$p_y2_y1 &
                       grid1$p_m2_m1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_u == p_y2_u & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u == grid1$p_m2_u+grid1$p_m1_u*grid1$p_m2_m1 &
                       grid1$p_y1_u == grid1$p_y2_u &
                       grid1$p_y2_y1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_m1 == 0 & p_m2_m1 == 0
(f.tab <- with(grid1[grid1$p_m1_u == grid1$p_m2_u+grid1$p_m1_u*grid1$p_m2_m1 &
                       grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u == grid1$p_m2_u+grid1$p_m1_u*grid1$p_m2_m1 &
                       grid1$p_y1_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(Min)))
prop.table(f.tab)


# p_m1_u == p_m2_u & p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u == grid1$p_m2_u+grid1$p_m1_u*grid1$p_m2_m1 &
                       grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_y1_m1 == 0 & p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_y1_m1 == 0 & p_m2_m1 == 0
(f.tab <- with(grid1[grid1$p_y1_u == grid1$p_y2_u+ grid1$p_y1_u*grid1$p_y2_y1 &
                       grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_y1_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_y1_u == grid1$p_y2_u+ grid1$p_y1_u*grid1$p_y2_y1 &
                       grid1$p_y1_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_y1_u == grid1$p_y2_u+ grid1$p_y1_u*grid1$p_y2_y1 &
                       grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_m1 == 0 & p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u == grid1$p_m2_u+grid1$p_m1_u*grid1$p_m2_m1 &
                       grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_y1_m1 == 0 & p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_y1_u == grid1$p_y2_u+ grid1$p_y1_u*grid1$p_y2_y1 &
                       grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_m1_u == p_m2_u & p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u == grid1$p_m2_u+grid1$p_m1_u*grid1$p_m2_m1 &
                       grid1$p_y1_u == grid1$p_y2_u+ grid1$p_y1_u*grid1$p_y2_y1 &
                       grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_m1_u == p_m2_u & p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u == grid1$p_m2_u+grid1$p_m1_u*grid1$p_m2_m1 &
                       grid1$p_y1_u == grid1$p_y2_u+ grid1$p_y1_u*grid1$p_y2_y1 &
                       grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(Min)))
prop.table(f.tab)

####### graphs
#overlaid densities of both biases
graph1 <- grid1 |>
  select(bias_ancova, diff, diff2, y2y1) |>
  as_tibble() |>
  round_df(4)
graph1[1:4,]
graph1_ancova <- graph1 |>  # Group by all columns
  group_by(bias_ancova, diff, diff2, y2y1) |>
  count()
graph1_ancova <- graph1_ancova |>
  mutate(perc = n/nrow(graph1))

graph1 <- grid1 |>
  select(bias_change, bias_naive, diff, diff2, y2y1) |>
  as_tibble() |>
  round_df(4)

graph1_change <- graph1 |>  # Group by all columns
  group_by(bias_change, diff, diff2, y2y1) |>
  count()
graph1_change <- graph1_change|>
  mutate(perc = n/nrow(graph1))

graph1 <- grid1 |>
  select(bias_change_1,  diff, diff2, y2y1) |>
  as_tibble() |>
  round_df(4)

graph1_change1 <- graph1 |>  # Group by all columns
  group_by(bias_change_1, diff, diff2, y2y1) |>
  count()
graph1_change1 <- graph1_change1|>
  mutate(perc = n/nrow(graph1))


graph1 <- grid1 |>
  select(bias_change_2, diff, diff2, y2y1) |>
  as_tibble() |>
  round_df(4)

graph1_change_2 <- graph1 |>  # Group by all columns
  group_by(bias_change_2, diff, diff2, y2y1) |>
  count()
graph1_change_2 <- graph1_change_2|>
  mutate(perc = n/nrow(graph1))
graph1 <- grid1 |>
  select(bias_change_1, diff, diff2, y2y1) |>
  as_tibble() |>
  round_df(4)
graph1_change_1 <- graph1 |>  # Group by all columns
  group_by(bias_change_1, diff, diff2, y2y1) |>
  count()
graph1_change_1 <- graph1_change_1|>
  mutate(perc = n/nrow(graph1))

graph1 <- grid1 |>
  select(bias_naive, diff, diff2, y2y1) |>
  as_tibble() |>
  round_df(4)

graph1_naive <- graph1 |>  # Group by all columns
  group_by(bias_naive, diff, diff2, y2y1) |>
  count()
graph1_naive<- graph1_naive |>
  mutate(perc = n/nrow(graph1))
write.csv(graph1_ancova,"ancova_sim1_1.csv")
write.csv(graph1_change,"change_sim1_1.csv")
write.csv(graph1_change1,"change1_sim1_1.csv")
write.csv(graph1_change_2,"change2_sim1_1.csv")
write.csv(graph1_naive,"naive_sim1_1.csv")
write.csv(graph1_change_1,"change1_sim1_1.csv")

p11 <- ggplot() +
  geom_density(data = graph1_ancova,aes( x=bias_ancova,weight = perc), linetype = "solid") +
  geom_density(data = graph1_change, aes(x=bias_change, weight = perc), linetype = "dashed", color = "red") +
  geom_density(data = graph1_change_2, aes(x=bias_change_2, weight = perc), linetype = "dashed", color = "blue") +
  geom_density(data = graph1_change_1, aes(x=bias_change_1, weight = perc), linetype = "dashed", color = "purple") +
  geom_density(data = graph1_naive, aes(x=bias_naive, weight = perc),linetype = "dotted") +
  labs(x="Raw bias", y="Percent") + 
  xlim(-1.7,1.7) + 
  theme_apa() 

tikz("biasdensity.tex",width = 4,height = 3, standAlone = TRUE)
pdf("figure_p11.pdf", height=6, width=6)
p11
dev.off()

grid1small1 <- grid1 |> filter(grid1$p_m1_u == grid1$p_m2_u+grid1$p_m1_u*grid1$p_m2_m1 &
                                 grid1$p_y1_u == grid1$p_y2_u+ grid1$p_y1_u*grid1$p_y2_y1 &
                                 grid1$p_y1_m1 == 0)
range(grid1small1$diff2)

graph1 <- grid1small1 |>
  select(bias_ancova, diff, diff2) |>
  as_tibble() |>
  round_df(4)

graph1_ancova <- graph1 |>  # Group by all columns
  group_by(bias_ancova, diff, diff2) |>
  count()
graph1_ancova <- graph1_ancova |>
  mutate(perc = n/nrow(graph1))

graph1 <- grid1small1 |>
  select(bias_change, bias_naive, diff, diff2) |>
  as_tibble() |>
  round_df(4)

graph1_change <- graph1 |>  # Group by all columns
  group_by(bias_change, diff, diff2) |>
  count()
graph1_change <- graph1_change|>
  mutate(perc = n/nrow(graph1))


graph1 <- grid1small1 |>
  select(bias_change_1, diff, diff2) |>
  as_tibble() |>
  round_df(4)

graph1_change_1 <- graph1|> # Group by all columns
  group_by(bias_change_1, diff, diff2) |>
  count()
graph1_change_1 <- graph1_change_1|>
  mutate(perc = n/nrow(graph1))

graph1 <- grid1small1 |>
  select(bias_change_2, diff, diff2) |>
  as_tibble() |>
  round_df(4)

graph1_change_2 <- graph1 |>  # Group by all columns
  group_by(bias_change_2, diff, diff2) |>
  count()
graph1_change_2 <- graph1_change_2|>
  mutate(perc = n/nrow(graph1))


graph1 <- grid1small1 |>
  select(bias_naive, diff, diff2) |>
  as_tibble() |>
  round_df(4)

graph1_naive <- graph1 |>  # Group by all columns
  group_by(bias_naive, diff, diff2) |>
  count()
graph1_naive<- graph1_naive |>
  mutate(perc = n/nrow(graph1))
write.csv(graph1_ancova,"ancova_sim1_2.csv")
write.csv(graph1_change,"change_sim1_2.csv")
write.csv(graph1_change_2,"change2_sim1_2.csv")
write.csv(graph1_change_1,"change1_sim1_2.csv")
write.csv(graph1_naive,"naive_sim1_2.csv")

grid1small1 <- grid1 |> filter(grid1$p_m1_u * grid1$p_m2_u +
                                 grid1$p_y1_u * grid1$p_y2_u == 0)
range(grid1small1$diff2)

graph1 <- grid1small1 |>
  select(bias_ancova, diff, diff2) |>
  as_tibble() |>
  round_df(4)

graph1_ancova <- graph1 |>  # Group by all columns
  group_by(bias_ancova, diff, diff2) |>
  count()
graph1_ancova <- graph1_ancova |>
  mutate(perc = n/nrow(graph1))

graph1 <- grid1small1 |>
  select(bias_change, bias_naive, diff, diff2) |>
  as_tibble() |>
  round_df(4)

graph1_change <- graph1 |>  # Group by all columns
  group_by(bias_change, diff, diff2) |>
  count()
graph1_change <- graph1_change|>
  mutate(perc = n/nrow(graph1))


graph1 <- grid1small1 |>
  select(bias_change_1, diff, diff2) |>
  as_tibble() |>
  round_df(4)

graph1_change_1 <- graph1|> # Group by all columns
  group_by(bias_change_1, diff, diff2) |>
  count()
graph1_change_1 <- graph1_change_1|>
  mutate(perc = n/nrow(graph1))

graph1 <- grid1small1 |>
  select(bias_change_2, diff, diff2) |>
  as_tibble() |>
  round_df(4)

graph1_change_2 <- graph1 |>  # Group by all columns
  group_by(bias_change_2, diff, diff2) |>
  count()
graph1_change_2 <- graph1_change_2|>
  mutate(perc = n/nrow(graph1))


graph1 <- grid1small1 |>
  select(bias_naive, diff, diff2) |>
  as_tibble() |>
  round_df(4)

graph1_naive <- graph1 |>  # Group by all columns
  group_by(bias_naive, diff, diff2) |>
  count()
graph1_naive<- graph1_naive |>
  mutate(perc = n/nrow(graph1))
write.csv(graph1_ancova,"ancova_sim1_3.csv")
write.csv(graph1_change,"change_sim1_3.csv")
write.csv(graph1_change_2,"change2_sim1_3.csv")
write.csv(graph1_change_1,"change1_sim1_3.csv")
write.csv(graph1_naive,"naive_sim1_3.csv")
grid1[grid1$diff2>2]
graph1 <- grid1small1$bias_ancova |>
  as_tibble() |>
  round_df(4)

graph1_ancova <- graph1 |>  # Group by all columns
  group_by(value) |>
  count()
graph1_ancova <- graph1_ancova |>
  mutate(perc = n/nrow(graph1))

graph1 <- grid1small1$bias_change |>
  as_tibble() |>
  round_df(4)

graph1_change <- graph1 |>  # Group by all columns
  group_by(value) |>
  count()
graph1_change <- graph1_change|>
  mutate(perc = n/nrow(graph1))


graph1 <- grid1small1$bias_change_2 |>
  as_tibble() |>
  round_df(4)

graph1_change_2 <- graph1 |>  # Group by all columns
  group_by(value) |>
  count()
graph1_change_2 <- graph1_change_2|>
  mutate(perc = n/nrow(graph1))

graph1 <- grid1small1$bias_naive |>
  as_tibble() |>
  round_df(4)

graph1_naive <- graph1 |>  # Group by all columns
  group_by(value) |>
  count()
graph1_naive<- graph1_naive |>
  mutate(perc = n/nrow(graph1))
write.csv(graph1_ancova,"ancova_sim1_2.csv")
write.csv(graph1_change,"change_sim1_2.csv")
write.csv(graph1_change_2,"change2_sim1_2.csv")
write.csv(graph1_naive,"naive_sim1_2.csv")
write.csv(grid1small1,"grid1small1.csv")


pdf("figure_p12.pdf", height=6, width=6)
p12
dev.off()

#numerical values from table
#both strictly larger than 0 (bias_diff_ind)
#and larger or equal than 0 (bias_diff_ind_eq)
ptable1 <- prop.table(table(grid1$bias_ind))
ptable1b <- prop.table(table(grid1$bias_ind_eq))
ptable1_2 <- prop.table(table(grid1$bias_ind_2))
ptable1b_2 <- prop.table(table(grid1$bias_ind_eq_2))

grid1small2 <- grid1 |> filter(p_m1_u == p_m2_u+p_m1_u*p_m2_m1 &
                                 p_y1_u == p_y2_u+p_y1_u*p_y2_y1 & 
                                 p_y1_m1 == 0 &
                                 p_m2_m1 == 0 &
                                 p_y2_y1 == 0)

write.csv(grid1small2,"grid1smal21.csv")

#constraining direct effect gamma = 0, no direct effect violation
graph3 <- pivot_longer(grid1small1, cols = c(bias_change_2, bias_naive, bias_ancova),
                       names_to = "bias_type", values_to = "bias_result")
p2 <- ggplot(graph3,aes(x=bias_result,y=..count../max(count), linetype = bias_type)) + 
  geom_density(bw=.1) + 
  labs(x="Raw bias", y="Percent") + 
  xlim(-1.7,1.7) + 
  theme_apa()
tikz("biasdensitynogamma.tex",width = 4,height = 3, standAlone = TRUE)
pdf("figure_p2.pdf", height=6, width=6)

p2
dev.off()
grid1small2 <- grid1 |> filter(p_m1_u == p_m2_u+p_m1_u*p_m2_m1 &
                                 p_y1_u == p_y2_u+p_y1_u*p_y2_y1 & 
                                 p_m2_m1 == 0 &
                                 p_y1_m1 == 0)
# bias_change_2 shows unbiased result only when the variance of gamma M is 2-2alpha^2
graph4 <- pivot_longer(grid1small2, cols = c(bias_change_2, bias_naive, bias_ancova, bias_change_1),
                       names_to = "bias_type", values_to = "bias_result")
p22 <- ggplot(graph4,aes(x=bias_result,y=..count../max(count), color = bias_type)) + 
  geom_density(bw=.1) + 
  labs(x="Raw bias", y="Percent") + 
  xlim(-1.7,1.7) + 
  theme_apa() + theme(legend.position = c(.2,.8))
tikz("biasdensitynogamma_2.tex",width = 4,height = 3, standAlone = TRUE)
pdf("figure_p22.pdf", height=6, width=6)

p22
dev.off()

table(grid1small1$bias_change_1)
ptable2 <- prop.table(table(grid1small1$bias_ind_2))
ptable2b <- prop.table(table(grid1small1$bias_ind_eq_2))

graph5 <- pivot_longer(grid1, cols = c(bias_change_2, bias_naive, bias_ancova, bias_change_1),
                       names_to = "bias_type", values_to = "bias_result")
ggplot(graph5,aes(x=bias_result,y=..count../max(count), color = bias_type)) + 
  geom_density(bw=.1) + 
  labs(x="Raw bias", y="Percent") + 
  xlim(-1.7,1.7) + 
  theme_apa() + theme(legend.position = c(.2,.8))



stopCluster(cl)
