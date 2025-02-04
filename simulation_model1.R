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


#### analytic bias from u
#### 7 grid points will only run on supercomputer
#### 9 grid points yield similar results
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

p_m2_x <- .2 # a path-fixed # these paths are fixed to reduce the running time. 
p_y2_x <- .2 # a path-fixed 
# As the x variable is considered to be controlled for in the product method when calculating the B path, it does not affect the outcome.
p_y2_m2 <- ng
p_y2_y1 <- ng 
p_m2_m1 <- ng 
p_y2_u <- ng 
p_y2_m1 <- ng 
p_m2_u <- ng 
p_m2_y1 <- ng 
p_y1_u <- ng 
p_m1_u <- ng 
p_y1_m1 <- ng 
p_gm_m2 <- 1
p_gm_m1 <- -1
p_gy_y2 <- 1
p_gy_y1 <- -1

grid1 <- data_frame()
#using data.table here, in the hope that it is faster

grid1 <- CJ(p_m2_x, p_y2_x, p_y2_m2, p_m2_m1, p_y2_y1, p_y2_u, p_m2_u,p_m2_y1 ,
            p_y1_u, p_y2_m1,
            p_m1_u, p_y1_m1, p_gm_m2, p_gm_m1, p_gy_y2, p_gy_y1)


grid1 <- as_tibble(grid)
setDT(grid1)

# Convert to a tibble for better display if desired

grid1[,`:=`(eps_m1 = sqrt(1-p_m1_u^2),
            eps_y1 = sqrt(1-p_y1_u^2-
                            p_y1_m1^2-
                            2*p_y1_u*p_m1_u*p_y1_m1),
            eps_m2 = sqrt(1-p_m2_m1^2-
                            p_m2_u^2-
                            p_m2_x^2-
                            p_m2_y1^2-
                            p_m2_m1*(2*p_m2_u*p_m1_u+
                                       2*p_m2_y1*p_y1_m1+
                                       2*p_m2_y1*p_y1_u*p_m1_u)-
                            p_m2_u*(2*p_m2_m1*p_m1_u+
                                      2*p_m2_y1*p_y1_m1*p_m1_u+
                                      2*p_m2_y1*p_y1_u)-
                            p_m2_y1*(2*p_m2_m1*p_y1_m1 +
                                       2*p_m2_m1*p_m1_u*p_y1_u+
                                       2*p_m2_u*p_m1_u*p_y1_m1+
                                       2*p_m2_u*p_y1_u)),
            eps_y2 = sqrt(1-p_y2_y1^2-
                            p_y2_m2^2-
                            p_y2_u^2-
                            p_y2_x^2-
                            p_y2_m1^2-
                            p_y2_m2*(2*p_y2_m1*p_m2_m1+
                                       2*p_y2_m1*p_y1_m1*p_m2_y1+
                                       2*p_y2_m1*p_m1_u*p_m2_u+
                                       2*p_y2_m1*p_m1_u*p_y1_u*p_m2_y1+
                                       2*p_y2_u*p_m1_u*p_m2_m1+
                                       2*p_y2_u*p_m1_u*p_y1_m1*p_m2_y1+
                                       2*p_y2_u*p_m2_u+
                                       2*p_y2_u*p_y1_u*p_m2_y1+
                                       2*p_y2_x*p_m2_x+
                                       2*p_y2_y1*p_m2_y1+
                                       2*p_y2_y1*p_y1_m1*p_m2_m1+
                                       2*p_y2_y1*p_y1_m1*p_m1_u*p_m2_u+
                                       2*p_y2_y1*p_y1_u*p_m1_u*p_m2_m1+
                                       2*p_y2_y1*p_y1_u*p_m2_u)-
                            p_y2_m1*(2*p_y2_m2*p_m2_m1+
                                       2*p_y2_m2*p_m2_u*p_m1_u+
                                       2*p_y2_m2*p_m2_y1*p_y1_m1+
                                       2*p_y2_m2*p_m2_y1*p_y1_u*p_m1_u+
                                       2*p_y2_u*p_m1_u+
                                       2*p_y2_y1*p_y1_m1+
                                       2*p_y2_y1*p_y1_u*p_m1_u)-
                            p_y2_y1*(2*p_y2_m1*p_y1_m1+
                                       2*p_y2_m1*p_m1_u*p_y1_u+
                                       2*p_y2_m2*p_m2_m1*p_y1_m1+
                                       2*p_y2_m2*p_m2_m1*p_m1_u*p_y1_u+
                                       2*p_y2_m2*p_m2_u*p_m1_u*p_y1_m1+
                                       2*p_y2_m2*p_m2_u*p_y1_u+
                                       2*p_y2_m2*p_m2_y1+
                                       2*p_y2_u*p_m1_u*p_y1_m1+
                                       2*p_y2_u*p_y1_u)-
                            p_y2_u*(2*p_y2_m1*p_m1_u+
                                      2*p_y2_m2*p_m2_m1*p_m1_u+
                                      2*p_y2_m2*p_m2_u+
                                      2*p_y2_m2*p_m2_y1*p_y1_m1*p_m1_u+
                                      2*p_y2_m2*p_m2_y1*p_y1_u+
                                      2*p_y2_y1*p_y1_m1*p_m1_u+
                                      2*p_y2_y1*p_y1_u)-
                            2*p_y2_m2*p_m2_x*p_y2_x))]


grid1 <- grid[eps_m1>=0 & eps_y1>=0 & eps_m2>=0 & eps_y2>=0]


# computing estimators
grid1[,`:=`(rm2y2 = p_y2_m2 +
              p_m2_m1*p_y1_m1*p_y2_y1 +
              p_m2_m1*p_y2_m1 +
              p_m2_m1*p_m1_u*p_y1_u*p_y2_y1 +
              p_m2_m1*p_m1_u*p_y2_u +
              p_m2_u*p_m1_u*p_y1_m1*p_y2_y1 +
              p_m2_u*p_m1_u*p_y2_m1 +
              p_m2_u*p_y1_u*p_y2_y1 +
              p_m2_u*p_y2_u +
              p_m2_x*p_y2_x +
              p_m2_y1*p_y2_y1 +
              p_m2_y1*p_y1_m1*p_y2_m1 +
              p_m2_y1*p_y1_m1*p_m1_u*p_y2_u +
              p_m2_y1*p_y1_u*p_m1_u*p_y2_m1 +
              p_m2_y1*p_y1_u*p_y2_u ,
            rm1m2 = p_m2_m1 + 
              p_m1_u*p_m2_u +
              p_y1_m1*p_m2_y1 +
              p_m1_u*p_y1_u*p_m2_y1,
            rm1y2 = p_m2_m1*p_y2_m2 +
              p_y1_m1*p_m2_y1*p_y2_m2 +
              p_y1_m1*p_y2_y1 +
              p_y2_m1 +
              p_m1_u*p_m2_u*p_y2_m2 +
              p_m1_u*p_y1_u*p_m2_y1*p_y2_m2 +
              p_m1_u*p_y1_u*p_y2_y1 +
              p_m1_u*p_y2_u,
            ry1m2 = p_m2_y1 +
              p_m2_m1*p_y1_m1 +
              p_m2_m1*p_m1_u*p_y1_u +
              p_m2_u*p_m1_u*p_y1_m1 +
              p_m2_u*p_y1_u,
            ry1y2 = p_y2_m2*p_m2_y1 + #Y1 -> M2 -> Y2"
              p_y2_y1 +#"Y1 -> Y2"  
              p_y2_m2*p_m2_m1*p_y1_m1 +#"Y1 <- M1 -> M2 -> Y2"   
              p_y2_m1*p_y1_m1 +# "Y1 <- M1 -> Y2"    
              p_y2_m2*p_m2_m1*p_m1_u*p_y1_u +#"Y1 <- U -> M1 -> M2 -> Y2" 
              p_y2_m1*p_m1_u*p_y1_u +#"Y1 <- U -> M1 -> Y2" 
              p_y2_m2*p_m2_u*p_m1_u*p_y1_m1 +#     "Y1 <- M1 <- U -> M2 -> Y2" 
              p_y2_u*p_m1_u*p_y1_m1 +#"Y1 <- M1 <- U -> Y2"    
              p_y2_m2*p_m2_u*p_y1_u +# "Y1 <- U -> M2 -> Y2"      
              p_y2_u*p_y1_u, #"Y1 <- U -> Y2",
            rm1y1 = p_y1_m1 + p_m1_u*p_y1_u,
            rm1x = 0,
            ry1x = 0,
            rm2x = p_m2_x,
            ry2x = p_m2_x*p_y2_m2 +
              p_y2_x,
            rgygm = 
              p_gm_m1 * p_m2_m1 * p_y2_m2 * p_gy_y2 +                   # "GM <- M1 -> M2 -> Y2 -> GY"
              p_gm_m1 * p_y1_m1 * p_gy_y1 +                             # "GM <- M1 -> Y1 -> GY"
              p_gm_m1 * p_y1_m1 * p_m2_y1 * p_y2_m2 * p_gy_y2 +         # "GM <- M1 -> Y1 -> M2 -> Y2 -> GY"
              p_gm_m1 * p_y1_m1 * p_y2_y1 * p_gy_y2 +                   # "GM <- M1 -> Y1 -> Y2 -> GY"
              p_gm_m1 * p_y2_m1 * p_gy_y2 +                             # "GM <- M1 -> Y2 -> GY"
              p_gm_m1 * p_m1_u * p_m2_u * p_y2_m2 * p_gy_y2 +           # "GM <- M1 <- U -> M2 -> Y2 -> GY"
              p_gm_m1 * p_m1_u * p_y1_u * p_gy_y1 +                     # "GM <- M1 <- U -> Y1 -> GY"
              p_gm_m1 * p_m1_u * p_y1_u * p_m2_y1 * p_y2_m2 * p_gy_y2 + # "GM <- M1 <- U -> Y1 -> M2 -> Y2 -> GY"
              p_gm_m1 * p_m1_u * p_y1_u * p_y2_y1 * p_gy_y2 +           # "GM <- M1 <- U -> Y1 -> Y2 -> GY"
              p_gm_m1 * p_m1_u * p_y2_u * p_gy_y2 +                     # "GM <- M1 <- U -> Y2 -> GY"
              p_gm_m2 * p_y2_m2 * p_gy_y2 +                             # "GM <- M2 -> Y2 -> GY"
              p_gm_m2 * p_m2_m1 * p_y1_m1 * p_gy_y1 +                   # "GM <- M2 <- M1 -> Y1 -> GY"
              p_gm_m2 * p_m2_m1 * p_y1_m1 * p_y2_y1 * p_gy_y2 +         # "GM <- M2 <- M1 -> Y1 -> Y2 -> GY"
              p_gm_m2 * p_m2_m1 * p_y2_m1 * p_gy_y2 +                   # "GM <- M2 <- M1 -> Y2 -> GY"
              p_gm_m2 * p_m2_m1 * p_m1_u * p_y1_u * p_gy_y1 +           # "GM <- M2 <- M1 <- U -> Y1 -> GY"
              p_gm_m2 * p_m2_m1 * p_m1_u * p_y1_u * p_y2_y1 * p_gy_y2 + # "GM <- M2 <- M1 <- U -> Y1 -> Y2 -> GY"
              p_gm_m2 * p_m2_m1 * p_m1_u * p_y2_u * p_gy_y2 +           # "GM <- M2 <- M1 <- U -> Y2 -> GY"
              p_gm_m2 * p_m2_u * p_m1_u * p_y1_m1 * p_gy_y1 +           # "GM <- M2 <- U -> M1 -> Y1 -> GY"
              p_gm_m2 * p_m2_u * p_m1_u * p_y1_m1 * p_y2_y1 * p_gy_y2 + # "GM <- M2 <- U -> M1 -> Y1 -> Y2 -> GY"
              p_gm_m2 * p_m2_u * p_m1_u * p_y2_m1 * p_gy_y2 +           # "GM <- M2 <- U -> M1 -> Y2 -> GY"
              p_gm_m2 * p_m2_u * p_y1_u * p_gy_y1 +                     # "GM <- M2 <- U -> Y1 -> GY"
              p_gm_m2 * p_m2_u * p_y1_u * p_y2_y1 * p_gy_y2 +           # "GM <- M2 <- U -> Y1 -> Y2 -> GY"
              p_gm_m2 * p_m2_u * p_y2_u * p_gy_y2 +                     # "GM <- M2 <- U -> Y2 -> GY"
              p_gm_m2 * p_m2_y1 * p_gy_y1 +                             # "GM <- M2 <- Y1 -> GY"
              p_gm_m2 *  p_m2_y1 * p_y2_y1 * p_gy_y2 +                   # "GM <- M2 <- Y1 -> Y2 -> GY"
              p_gm_m2 *  p_m2_y1 * p_y1_m1 * p_y2_m1 * p_gy_y2 +         # "GM <- M2 <- Y1 <- M1 -> Y2 -> GY"
              p_gm_m2 *  p_m2_y1 * p_m1_u * p_y2_u * p_gy_y2 +           # "GM <- M2 <- Y1 <- M1 <- U -> Y2 -> GY"
              p_gm_m2 *  p_m2_y1 * p_y1_u * p_m1_u * p_y2_m1 * p_gy_y2 + # "GM <- M2 <- Y1 <- U -> M1 -> Y2 -> GY"
              p_gm_m2 *  p_m2_y1 * p_y1_u * p_y2_u * p_gy_y2 + # "GM <- M2 <- Y1 <- U -> Y2 -> GY"
              p_gm_m2 *  p_m2_x * p_y2_x * p_gy_y2
            ,
            rgym2 = 
              p_y2_m2 * p_gy_y2 +                                         # "M2 -> Y2 -> GY"
              p_m2_m1 * p_y1_m1 * p_gy_y1 +                               # "M2 <- M1 -> Y1 -> GY"
              p_m2_m1 * p_y1_m1 * p_y2_y1 * p_gy_y2 +                     # "M2 <- M1 -> Y1 -> Y2 -> GY"
              p_m2_m1 * p_y2_m1 * p_gy_y2 +                               # "M2 <- M1 -> Y2 -> GY"
              p_m2_m1 * p_m1_u * p_y1_u * p_gy_y1 +                       # "M2 <- M1 <- U -> Y1 -> GY"
              p_m2_m1 * p_m1_u * p_y1_u * p_y2_y1 * p_gy_y2 +             # "M2 <- M1 <- U -> Y1 -> Y2 -> GY"
              p_m2_m1 * p_m1_u * p_y2_u * p_gy_y2 +                       # "M2 <- M1 <- U -> Y2 -> GY"
              p_m2_u * p_m1_u * p_y1_m1 * p_gy_y1 +                       # "M2 <- U -> M1 -> Y1 -> GY"
              p_m2_u * p_m1_u * p_y1_m1 * p_y2_y1 * p_gy_y2 +             # "M2 <- U -> M1 -> Y1 -> Y2 -> GY"
              p_m2_u * p_m1_u * p_y2_m1 * p_gy_y2 +                       # "M2 <- U -> M1 -> Y2 -> GY"
              p_m2_u * p_y1_u * p_gy_y1 +                                 # "M2 <- U -> Y1 -> GY"
              p_m2_u * p_y1_u * p_y2_y1 * p_gy_y2 +                       # "M2 <- U -> Y1 -> Y2 -> GY"
              p_m2_u * p_y2_u * p_gy_y2 +                                 # "M2 <- U -> Y2 -> GY"
              p_m2_y1 * p_gy_y1 +                                         # "M2 <- Y1 -> GY"
              p_m2_y1 * p_y2_y1 * p_gy_y2 +                               # "M2 <- Y1 -> Y2 -> GY"
              p_m2_y1 * p_y1_m1 * p_y2_m1 * p_gy_y2 +                     # "M2 <- Y1 <- M1 -> Y2 -> GY"
              p_m2_y1 * p_m1_u * p_y2_u * p_gy_y2 +                       # "M2 <- Y1 <- M1 <- U -> Y2 -> GY"
              p_m2_y1 * p_y1_u * p_m1_u * p_y2_m1 * p_gy_y2 +             # "M2 <- Y1 <- U -> M1 -> Y2 -> GY"
              p_m2_y1 * p_y1_u * p_y2_u * p_gy_y2 +
              p_m2_x * p_y2_x * p_gy_y2,
            rgyx = p_gy_y2*p_y2_m2*p_m2_x+
              p_gy_y2*p_y2_x,
            rgmx = p_gm_m2*p_m2_x )]


grid1[,`:=`( ry1m2_m1 = (ry1m2 - rm1y1 * rm1m2) / sqrt((1 - rm1y1^2) * (1 - rm1m2^2)),
             ry1y2_m1 = (ry1y2 - rm1y1 * rm1y2) / sqrt((1 - rm1y1^2) * (1 - rm1y2^2)),
             rm2y2_m1 = (rm2y2 - rm1m2 * rm1y2) / sqrt((1 - rm1m2^2) * (1 - rm1y2^2)),
             sd_m2_m1 = sqrt(1 - rm1m2^2),
             sd_y2_m1 = sqrt(1 - rm1y2^2),
             sd_y1_m1 = sqrt(1 - rm1y1^2),
             ry1x_m1 = (ry1x - rm1y1 * rm1x) / sqrt((1 - rm1y1^2) * (1 - rm1x^2)),
             rxy2_m1 = (ry2x - rm1x * rm1y2) / sqrt((1 - rm1x^2) * (1 - rm1y2^2)),
             sd_x_m1 = sqrt(1 - rm1x^2),
             rxm2_m1 = (rm2x - rm1x * rm1m2) / sqrt((1 - rm1x^2) * (1 - rm1m2^2)),
             sd_m2_m1y1 = sqrt(1 - (ry1m2^2 + rm1m2^2 - 2 * ry1m2 * rm1m2 * rm1y1) / (1 - rm1y1^2)),
             sd_y2_m1y1 = sqrt(1 - (ry1y2^2 + rm1y2^2 - 2 * ry1y2 * rm1y2 * rm1y1) / (1 - rm1y1^2)),
             sd_x_m1y1 = sqrt(1 - (ry1x^2 + rm1x^2 - 2 * ry1x * rm1x * rm1y1) / (1 - rm1y1^2)),
             bgygm_x = ((rgygm - rgyx * rgmx) / (1 - rgmx^2)),
             bgym2_x = ((rgym2 - rgyx * rm2x) / (1 - rm2x^2)),
             rgygm_x = (rgygm - rgyx * rgmx) / sqrt((1 - rgyx^2) * (1 - rgmx^2)),
             rgym2_x = (rgym2 - rgyx * rm2x) / sqrt((1 - rgyx^2) * (1 - rm2x^2)),
             by2m2_x = (rm2y2 - ry2x * rm2x) / ((1 -  rm2x^2)),
             rm1m2_x = (rm1m2 - rm1x * rm2x) / sqrt((1 - rm1x^2) * (1 - rm2x^2))
)]


grid1[,`:=`( rm2y2_m1y1 = (rm2y2_m1 - ry1m2_m1 * ry1y2_m1) / sqrt((1 - ry1m2_m1^2) * (1 - ry1y2_m1^2)),
             rxy2_m1y1 = (rxy2_m1 - ry1x_m1 * ry1y2_m1) / sqrt((1 - ry1x_m1^2) * (1 - ry1y2_m1^2)),
             rxm2_m1y1 = (rxm2_m1 - ry1x_m1 * ry1m2_m1) / sqrt((1 - ry1x_m1^2) * (1 - ry1m2_m1^2)))]

grid1[,`:=`(rm2y2_xm1y1 = (rm2y2_m1y1- rxm2_m1y1 *rxy2_m1y1) /
                             sqrt((1 - rxm2_m1y1^2) * (1 - rxy2_m1y1^2))
)]
grid1[,`:=`(bm2y2_xm1y1 = ((rm2y2_m1y1- rxm2_m1y1 *rxy2_m1y1) /
                             (1 - rxm2_m1y1 ^2 ))* sd_y2_m1y1 / sd_m2_m1y1
)]




grid1 <- na.omit(grid1)

# Remove cases where the absolute correlation is greater than 1, as this is impossible.
grid1 <- grid1[!(abs(rm2y2) > 1 |
                   abs(rm1m2) > 1 |
                   abs(rm1y2) > 1 |
                   abs(ry1m2) > 1 |
                   abs(ry1y2) > 1 |
                   abs(rm1y1) > 1 |
                   abs(rm2x) > 1 |
                   abs(ry2x) > 1|
                   abs(rgygm_x) > 1|
                   abs(rgym2_x) > 1|
                   abs(rxm2_m1) > 1|
                   abs(rm2y2_m1y1) > 1|
                   abs(rxy2_m1y1) > 1|
                   abs(rxm2_m1y1) > 1|
                   abs(rm2y2_xm1y1)>1)]

grid1[,`:=`(change = p_m2_x*bgygm_x/(2-2*rm1m2_x),
            change_2 = p_m2_x*bgygm_x/(1-rm1m2_x),
            change_1 = p_m2_x*bgym2_x,
            ancova = p_m2_x*bm2y2_xm1y1 ,
            truth = p_m2_x*p_y2_m2,
            naive =p_m2_x*by2m2_x)]

range(grid1$ancova)

# Convert to a tibble for better display if desired

# computing bias
grid1[,`:=`(bias_change = change - truth, ### 2 wave change score that people usually use
            bias_change_2 = change_2 - truth, ### 2 wave change score multiplied by 2
            bias_ancova = ancova - truth,
            bias_naive = naive - truth,
            bias_change_1 = change_1 - truth)] ### change score with one pre-test (outcome)
# add the column showing the degree of violation of the ignorability assumption
grid1$diff <- abs(grid1$p_m2_u*grid1$p_y2_u)

grid1 <- grid1[,-(13:42)] # Remove unnecessary columns for faster simulation

grid1[,27:39] <- round_df(grid1[,27:39], 4)
grid1$Min_1 <- colnames(grid1[,34:38])[minCol(abs(grid1[,34:38]), ties.method = "first")]
grid1$Min_2 <- colnames(grid1[,34:38])[minCol(abs(grid1[,34:38]), ties.method = "last")]
grid1<- grid1 |> mutate(Min = ifelse(Min_1 != Min_2, "Tie", Min_1))


####### Save data

# Graph under the no time varying confounder assumption only

# save dataset
graph1 <- grid1 |>
  select(bias_ancova, diff) |>
  as_tibble() |>
  round_df(4)

graph1_ancova <- graph1 |>  # Group by all columns
  group_by(bias_ancova, diff) |>
  count()


graph1 <- grid1 |>
  select(bias_change, diff) |>
  as_tibble() |>
  round_df(4)

graph1_change <- graph1 |>  # Group by all columns
  group_by(bias_change, diff) |>
  count()


graph1 <- grid1 |>
  select(bias_change_1, diff) |>
  as_tibble() |>
  round_df(4)

graph1_change1 <- graph1 |>  # Group by all columns
  group_by(bias_change_1, diff) |>
  count()



graph1 <- grid1 |>
  select(bias_change_2, diff) |>
  as_tibble() |>
  round_df(4)

graph1_change_2 <- graph1 |>  # Group by all columns
  group_by(bias_change_2, diff) |>
  count()


graph1 <- grid1 |>
  select(bias_change_1, diff) |>
  as_tibble() |>
  round_df(4)
graph1_change_1 <- graph1 |>  # Group by all columns
  group_by(bias_change_1, diff) |>
  count()


graph1 <- grid1 |>
  select(bias_naive,diff) |>
  as_tibble() |>
  round_df(4)

graph1_naive <- graph1 |>  # Group by all columns
  group_by(bias_naive, diff) |>
  count()
graph1_naive<- graph1_naive |>
  mutate(perc = n/nrow(graph1))

write.csv(graph1_ancova,"ancova_sim1_r1.csv")
write.csv(graph1_change,"change_sim1_r1.csv")
write.csv(graph1_change1,"change1_sim1_r1.csv")
write.csv(graph1_change_2,"change2_sim1_r1.csv")
write.csv(graph1_naive,"naive_sim1_r1.csv")
write.csv(graph1_change_1,"change1_sim1_r1.csv")



grid1small1 <- grid1 |> filter(p_y1_u + p_m1_u*p_y1_m1== p_y2_u+p_y1_u*p_y2_y1+
                                 p_m1_u*p_y2_m1+p_m1_u*p_y1_m1*p_y2_y1&
                                 p_y1_m1 + p_y1_u*p_m1_u ==
                                 p_y1_m1*p_y2_y1 + p_y2_m1 + p_m1_u*p_y1_u*p_y2_y1+
                                 p_m1_u*p_y2_u &
                                 p_m2_y1 == 0
)


graph1 <- grid1small1 |>
  select(bias_ancova) |>
  as_tibble() |>
  round_df(4)

graph1_ancova <- graph1 |>  # Group by all columns
  group_by(bias_ancova) |>
  count()


graph1 <- grid1small1 |>
  select(bias_change) |>
  as_tibble() |>
  round_df(4)

graph1_change <- graph1 |>  # Group by all columns
  group_by(bias_change) |>
  count()



graph1 <- grid1small1 |>
  select(bias_change_1) |>
  as_tibble() |>
  round_df(4)

graph1_change_1 <- graph1|> # Group by all columns
  group_by(bias_change_1) |>
  count()


graph1 <- grid1small1 |>
  select(bias_change_2) |>
  as_tibble() |>
  round_df(4)

graph1_change_2 <- graph1 |>  # Group by all columns
  group_by(bias_change_2) |>
  count()



graph1 <- grid1small1 |>
  select(bias_naive) |>
  as_tibble() |>
  round_df(4)

graph1_naive <- graph1 |>  # Group by all columns
  group_by(bias_naive) |>
  count()

write.csv(graph1_ancova,"ancova_sim1_rc.csv")
write.csv(graph1_change,"change_sim1_rc.csv")
write.csv(graph1_change_2,"change2_sim1_rc.csv")
write.csv(graph1_change_1,"change1_sim1_rc.csv")
write.csv(graph1_naive,"naive_sim1_rc.csv")

grid1small1 <- grid1 |> filter(grid1$p_m2_u*
                                 grid1$p_y2_u==0)
range(grid1small1$bias_ancova)
range(grid1small1$diff2)
graph1 <- grid1small1 |>
  select(bias_ancova) |>
  as_tibble() |>
  round_df(4)

graph1_ancova <- graph1 |>  # Group by all columns
  group_by(bias_ancova) |>
  count()


graph1 <- grid1small1 |>
  select(bias_change) |>
  as_tibble() |>
  round_df(4)

graph1_change <- graph1 |>  # Group by all columns
  group_by(bias_change) |>
  count()



graph1 <- grid1small1 |>
  select(bias_change_1) |>
  as_tibble() |>
  round_df(4)

graph1_change_1 <- graph1|> # Group by all columns
  group_by(bias_change_1) |>
  count()


graph1 <- grid1small1 |>
  select(bias_change_2) |>
  as_tibble() |>
  round_df(4)

graph1_change_2 <- graph1 |>  # Group by all columns
  group_by(bias_change_2) |>
  count()



graph1 <- grid1small1 |>
  select(bias_naive) |>
  as_tibble() |>
  round_df(4)

graph1_naive <- graph1 |>  # Group by all columns
  group_by(bias_naive) |>
  count()

write.csv(graph1_ancova,"ancova_sim1_ri.csv")
write.csv(graph1_change,"change_sim1_ri.csv")
write.csv(graph1_change_2,"change2_sim1_ri.csv")
write.csv(graph1_change_1,"change1_sim1_ri.csv")
write.csv(graph1_naive,"naive_sim1_ri.csv")

stopCluster(cl)
