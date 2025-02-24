library(profvis)
library(tictoc)
library(plyr)
library(doParallel)
cl<-49
registerDoParallel(cl)
library(dagitty) # to draw DAGs
library(raster)
library(ggplot2)
library(ggthemes)
library(patchwork)
library(tikzDevice)
library(matrixcalc)
library(data.table)
library(dplyr)
library(tidyr)
library(jtools)
library(cgwtools)
library(forestmangr)
library(dendroTools)
#setDTthreads(percent = 90)

ng <- list()
ng <- append(.1, ng)
ng <- append(sqrt(ng[[1]]),ng)
ng <- append(sqrt(ng[[1]]),ng)
#ng <- append(sqrt(ng[[1]]),ng)
ng <- append(0,ng)
ng <- append(-1*ng[[2]],ng)
ng <- append(-1*ng[[4]],ng)
ng <- append(-1*ng[[6]],ng)
#ng <- append(-1*ng[[8]],ng)
ng <- as.numeric(ng)

p_m2_x <- .2 # a path-fixed
p_y2_x <- .2
p_m2_m1 <- ng
p_y2_y1 <- ng
p_y2_u2 <- ng
p_m2_u2 <- ng
p_u2_u1 <- ng
p_m1_u1 <- ng
p_y1_u1 <- ng
p_y1_m1 <- ng
p_m2_y1 <- ng
p_y2_m1 <- ng
p_gm_m2 <- 1
p_gm_m1 <- -1
p_gy_y2 <- 1
p_gy_y1 <- -1

grid1 <- data_frame()
grid <- data.table()

for (i in 1:2){
  p_y2_m2 <- ng[[i]]
  grid <- CJ(p_m2_x, p_y2_x, p_y2_m2, p_m2_m1, p_y2_y1, p_y2_u2, p_m2_u2, p_y1_u1,
             p_m1_u1, p_u2_u1, p_y1_m1, p_gm_m2, p_gm_m1, p_gy_y2, p_gy_y1, p_y2_m1, p_m2_y1)
  grid[,`:=`(eps_m1 = sqrt(1-p_m1_u1^2),
             eps_y1 = sqrt(1-p_y1_u1^2-
                             p_y1_m1^2-
                             2*p_y1_u1*p_m1_u1*p_y1_m1),
             eps_m2 = sqrt(1-p_m2_m1^2-
                             p_m2_u2^2-
                             p_m2_x^2-
                             p_m2_y1^2-
                             p_m2_m1*(2*p_m2_u2*p_u2_u1*p_m1_u1+
                                        2*p_m2_y1*p_y1_m1+
                                        2*p_m2_y1*p_y1_u1*p_m1_u1)-
                             p_m2_u2*(2*p_m2_m1*p_m1_u1*p_u2_u1+
                                        2*p_m2_y1*p_y1_m1*p_m1_u1*p_u2_u1+
                                        2*p_m2_y1*p_y1_u1*p_u2_u1)-
                             p_m2_y1*(2*p_m2_m1*p_y1_m1 +
                                        2*p_m2_m1*p_m1_u1*p_y1_u1+
                                        2*p_m2_u2*p_u2_u1*p_m1_u1*p_y1_m1+
                                        2*p_m2_u2*p_u2_u1*p_y1_u1)),
             eps_y2 = sqrt(1-p_y2_y1^2-
                             p_y2_m2^2-
                             p_y2_u2^2-
                             p_y2_x^2-
                             p_y2_m1^2-
                             p_y2_m2*(2*p_y2_m1*p_m2_m1+
                                        2*p_y2_m1*p_y1_m1*p_m2_y1+
                                        2*p_y2_m1*p_m1_u1*p_u2_u1*p_m2_u2+
                                        2*p_y2_m1*p_m1_u1*p_y1_u1*p_m2_y1+
                                        2*p_y2_u2*p_u2_u1*p_m1_u1*p_m2_m1+
                                        2*p_y2_u2*p_u2_u1*p_m1_u1*p_y1_m1*p_m2_y1+
                                        2*p_y2_u2*p_m2_u2+
                                        2*p_y2_u2*p_u2_u1*p_y1_u1*p_m2_y1+
                                        2*p_y2_x*p_m2_x+
                                        2*p_y2_y1*p_m2_y1+
                                        2*p_y2_y1*p_y1_m1*p_m2_m1+
                                        2*p_y2_y1*p_y1_m1*p_m1_u1*p_u2_u1*p_m2_u2+
                                        2*p_y2_y1*p_y1_u1*p_m1_u1*p_m2_m1+
                                        2*p_y2_y1*p_y1_u1*p_u2_u1*p_m2_u2)-
                             p_y2_m1*(2*p_y2_m2*p_m2_m1+
                                        2*p_y2_m2*p_m2_u2*p_u2_u1*p_m1_u1+
                                        2*p_y2_m2*p_m2_y1*p_y1_m1+
                                        2*p_y2_m2*p_m2_y1*p_y1_u1*p_m1_u1+
                                        2*p_y2_u2*p_u2_u1*p_m1_u1+
                                        2*p_y2_y1*p_y1_m1+
                                        2*p_y2_y1*p_y1_u1*p_m1_u1)-
                             p_y2_y1*(2*p_y2_m1*p_y1_m1+
                                        2*p_y2_m1*p_m1_u1*p_y1_u1+
                                        2*p_y2_m2*p_m2_m1*p_y1_m1+
                                        2*p_y2_m2*p_m2_m1*p_m1_u1*p_y1_u1+
                                        2*p_y2_m2*p_m2_u2*p_u2_u1*p_m1_u1*p_y1_m1+
                                        2*p_y2_m2*p_m2_u2*p_u2_u1*p_y1_u1+
                                        2*p_y2_m2*p_m2_y1+
                                        2*p_y2_u2*p_u2_u1*p_m1_u1*p_y1_m1+
                                        2*p_y2_u2*p_u2_u1*p_y1_u1)-
                             p_y2_u2*(2*p_y2_m1*p_m1_u1*p_u2_u1+
                                        2*p_y2_m2*p_m2_m1*p_m1_u1*p_u2_u1+
                                        2*p_y2_m2*p_m2_u2+
                                        2*p_y2_m2*p_m2_y1*p_y1_m1*p_m1_u1*p_u2_u1+
                                        2*p_y2_m2*p_m2_y1*p_y1_u1*p_u2_u1+
                                        2*p_y2_y1*p_y1_m1*p_m1_u1*p_u2_u1+
                                        2*p_y2_y1*p_y1_u1*p_u2_u1)-
                             2*p_y2_m2*p_m2_x*p_y2_x))]
  grid <- grid[eps_m1>=0 & eps_y1>=0 & eps_m2>=0 & eps_y2>=0]
  grid1 <- rbind(grid1,grid)
}
grid1[1:10]
rm(grid)
grid1 <- na.omit(grid1)
# computing estimators
grid1[,`:=`(rm2y2 = p_y2_m2 +  # "M2 -> Y2"
              p_m2_m1 * p_y1_m1 * p_y2_y1 +  # "M2 <- M1 -> Y1 -> Y2"
              p_m2_m1 * p_y2_m1 +  # "M2 <- M1 -> Y2"
              p_m2_m1 * p_m1_u1 * p_u2_u1 * p_y2_u2 +  # "M2 <- M1 <- U1 -> U2 -> Y2"
              p_m2_m1 * p_m1_u1 * p_y1_u1 * p_y2_y1 +  # "M2 <- M1 <- U1 -> Y1 -> Y2"
              p_m2_u2 * p_y2_u2 +  # "M2 <- U2 -> Y2"
              p_m2_u2 * p_u2_u1 * p_m1_u1 * p_y1_m1 * p_y2_y1 +  # "M2 <- U2 <- U1 -> M1 -> Y1 -> Y2"
              p_m2_u2 * p_u2_u1 * p_m1_u1 * p_y2_m1 +  # "M2 <- U2 <- U1 -> M1 -> Y2"
              p_m2_u2 * p_u2_u1 * p_y1_u1 * p_y2_y1 +  # "M2 <- U2 <- U1 -> Y1 -> Y2"
              p_m2_x * p_y2_x +  # "M2 <- X -> Y2"
              p_m2_y1 * p_y2_y1 +  # "M2 <- Y1 -> Y2"
              p_m2_y1 * p_y1_m1 * p_y2_m1 +  # "M2 <- Y1 <- M1 -> Y2"
              p_m2_y1 * p_y1_m1 * p_m1_u1 * p_u2_u1 * p_y2_u2 +  # "M2 <- Y1 <- M1 <- U1 -> U2 -> Y2"
              p_m2_y1 * p_y1_u1 * p_m1_u1 * p_y2_m1 +  # "M2 <- Y1 <- U1 -> M1 -> Y2"
              p_m2_y1 * p_y1_u1 * p_u2_u1 * p_y2_u2 , # "M2 <- Y1 <- U1 -> U2 -> Y2" 
            rm1m2 = p_m2_m1 + #"M1 -> M2" 
              p_m1_u1 * p_u2_u1 * p_m2_u2 +#  "M1 <- U1 -> U2 -> M2" 
              p_y1_m1 * p_m2_y1 + # "M1 -> Y1 -> M2"   
              p_m1_u1 * p_y1_u1 * p_m2_y1,#"M1 <- U1 -> Y1 -> M2"
            rm1y2 = p_m2_m1*p_y2_m2 +  # "M1 -> M2 -> Y2"
              p_y1_m1 * p_y2_y1 +  # "M1 -> Y1 -> Y2"
              p_m1_u1 * p_u2_u1 * p_y2_u2 +  # "M1 <- U1 -> U2 -> Y2"
              p_m1_u1 * p_y1_u1 * p_y2_y1 +  # "M1 <- U1 -> Y1 -> Y2"
              p_m1_u1 * p_u2_u1 * p_m2_u2 * p_y2_m2 +  # "M1 <- U1 -> U2 -> M2 -> Y2"
              p_y1_m1 * p_m2_y1 * p_y2_m2 +  # "M1 -> Y1 -> M2 -> Y2"
              p_y2_m1 +  # "M1 -> Y2"
              p_m1_u1 * p_y1_u1 * p_m2_y1 * p_y2_m2,  # "M1 <- U1 -> Y1 -> M2 -> Y2"
            ry1m2 = p_m2_m1*p_y1_m1 +  # "Y1 -> M2"
              p_m2_m1*p_m1_u1*p_y1_u1 +  # "Y1 <- M1 -> M2"
              p_y1_m1*p_m1_u1*p_u2_u1*p_m2_u2 +  # "Y1 <- M1 <- U1 -> U2 -> M2"
              p_y1_u1*p_u2_u1*p_m2_u2 +  # "Y1 <- U1 -> U2 -> M2"
              p_m2_y1,  # "Y1 <- U1 -> M1 -> M2"
            ry1y2 = p_y2_m2*p_m2_y1 + #Y1 -> M2 -> Y2"
              p_y2_y1 +#"Y1 -> Y2"  
              p_y2_m2*p_m2_m1*p_y1_m1 +#"Y1 <- M1 -> M2 -> Y2"   
              p_y2_m1*p_y1_m1 +# "Y1 <- M1 -> Y2"    
              p_y2_m2*p_m2_m1*p_m1_u1*p_y1_u1 +#"Y1 <- U -> M1 -> M2 -> Y2" 
              p_y2_m1*p_m1_u1*p_y1_u1 +#"Y1 <- U -> M1 -> Y2" 
              p_y2_m2*p_m2_u2 *p_u2_u1 *p_m1_u1*p_y1_m1 +#     "Y1 <- M1 <- U -> M2 -> Y2" 
              p_y2_u2 *p_u2_u1 *p_m1_u1*p_y1_m1 +#"Y1 <- M1 <- U -> Y2"    
              p_y2_m2*p_m2_u2 *p_u2_u1 *p_y1_u1 +# "Y1 <- U -> M2 -> Y2"      
              p_y2_u2 *p_u2_u1 *p_y1_u1, #"Y1 <- U -> Y2",
            rm1y1 = p_y1_m1 + 
              p_m1_u1*p_y1_u1,
            rm1x = 0,
            ry1x = 0,
            rm2x = p_m2_x,
            ry2x = p_m2_x*p_y2_m2 +
              p_y2_x,
            rgygm = 
              p_gy_y2 * p_y2_m2 * p_m2_m1 * p_gm_m1 +  # "GM <- M1 -> M2 -> Y2 -> GY"
              p_gy_y1 * p_y1_m1 * p_gm_m1 +  # "GM <- M1 -> Y1 -> GY"
              p_gy_y2 * p_y2_m2 * p_m2_y1 * p_y1_m1 * p_gm_m1 +  # "GM <- M1 -> Y1 -> M2 -> Y2 -> GY"
              p_gy_y2 * p_y2_y1 * p_y1_m1 * p_gm_m1 +  # "GM <- M1 -> Y1 -> Y2 -> GY"
              p_gy_y2 * p_y2_m1 * p_gm_m1 +  # "GM <- M1 -> Y2 -> GY"
              p_gy_y2 * p_y2_m2 * p_m2_u2 * p_u2_u1 * p_m1_u1 * p_gm_m1 +  # "GM <- M1 <- U1 -> U2 -> M2 -> Y2 -> GY"
              p_gy_y2 * p_y2_u2 * p_u2_u1 * p_m1_u1 * p_gm_m1 +  # "GM <- M1 <- U1 -> U2 -> Y2 -> GY"
              p_gy_y1 * p_y1_u1 * p_m1_u1 * p_gm_m1 +  # "GM <- M1 <- U1 -> Y1 -> GY"
              p_gy_y2 * p_y2_m2 * p_m2_y1 * p_y1_u1 * p_m1_u1 * p_gm_m1 +  # "GM <- M1 <- U1 -> Y1 -> M2 -> Y2 -> GY"
              p_gy_y2 * p_y2_y1 * p_y1_u1 * p_m1_u1 * p_gm_m1 +  # "GM <- M1 <- U1 -> Y1 -> Y2 -> GY"
              p_gy_y2 * p_y2_m2 * p_gm_m2 +  # "GM <- M2 -> Y2 -> GY"
              p_gy_y1 * p_y1_m1 * p_m2_m1 * p_gm_m2 +  # "GM <- M2 <- M1 -> Y1 -> GY"
              p_gy_y2 * p_y2_y1 * p_y1_m1 * p_m2_m1 * p_gm_m2 +  # "GM <- M2 <- M1 -> Y1 -> Y2 -> GY"
              p_gy_y2 * p_y2_m1 * p_m2_m1 * p_gm_m2 +  # "GM <- M2 <- M1 -> Y2 -> GY"
              p_gy_y2 * p_y2_u2 * p_u2_u1 * p_m1_u1 * p_m2_m1 * p_gm_m2 +  # "GM <- M2 <- M1 <- U1 -> U2 -> Y2 -> GY"
              p_gy_y1 * p_y1_u1 * p_m1_u1 * p_m2_m1 * p_gm_m2 +  # "GM <- M2 <- M1 <- U1 -> Y1 -> GY"
              p_gy_y2 * p_y2_y1 * p_y1_u1 * p_m1_u1 * p_m2_m1 * p_gm_m2 +  # "GM <- M2 <- M1 <- U1 -> Y1 -> Y2 -> GY"
              p_gy_y2 * p_y2_u2 * p_m2_u2 * p_gm_m2 +  # "GM <- M2 <- U2 -> Y2 -> GY"
              p_gy_y1 * p_y1_m1 * p_m1_u1 * p_u2_u1 * p_m2_u2 * p_gm_m2 +  # "GM <- M2 <- U2 <- U1 -> M1 -> Y1 -> GY"
              p_gy_y2 * p_y2_y1 * p_y1_m1 * p_m1_u1 * p_u2_u1 * p_m2_u2 * p_gm_m2 +  # "GM <- M2 <- U2 <- U1 -> M1 -> Y1 -> Y2 -> GY"
              p_gy_y2 * p_y2_m1 * p_m1_u1 * p_u2_u1 * p_m2_u2 * p_gm_m2 +  # "GM <- M2 <- U2 <- U1 -> M1 -> Y2 -> GY"
              p_gy_y1 * p_y1_u1 * p_u2_u1 * p_m2_u2 * p_gm_m2 +  # "GM <- M2 <- U2 <- U1 -> Y1 -> GY"
              p_gy_y2 * p_y2_y1 * p_y1_u1 * p_u2_u1 * p_m2_u2 * p_gm_m2 +  # "GM <- M2 <- U2 <- U1 -> Y1 -> Y2 -> GY"
              p_gy_y2 * p_y2_x * p_m2_x * p_gm_m2 +  # "GM <- M2 <- X -> Y2 -> GY"
              p_gy_y1 * p_m2_y1 * p_gm_m2 +  # "GM <- M2 <- Y1 -> GY"
              p_gy_y2 * p_y2_y1 * p_m2_y1 * p_gm_m2 +  # "GM <- M2 <- Y1 -> Y2 -> GY"
              p_gy_y2 * p_y2_m1 * p_m2_m1 * p_m2_y1 * p_gm_m2 +  # "GM <- M2 <- Y1 <- M1 -> Y2 -> GY"
              p_gy_y2 * p_y2_u2 * p_u2_u1 * p_m1_u1 * p_m2_m1 * p_m2_y1 * p_gm_m2 +  # "GM <- M2 <- Y1 <- M1 <- U1 -> U2 -> Y2 -> GY"
              p_gy_y2 * p_y2_m1 * p_m1_u1 * p_y1_u1 * p_m2_y1 * p_gm_m2 +  # "GM <- M2 <- Y1 <- U1 -> M1 -> Y2 -> GY"
              p_gy_y2 * p_y2_u2 * p_u2_u1 * p_y1_u1 * p_m2_y1 * p_gm_m2,  # "GM <- M2 <- Y1 <- U1 -> U2 -> Y2 -> GY"
            rgym2 = 
              p_gy_y2 * p_y2_m2 +  # "M2 -> Y2 -> GY"
              p_gy_y1 * p_y1_m1 * p_m2_m1 +  # "M2 <- M1 -> Y1 -> GY"
              p_gy_y2 * p_y2_y1 * p_y1_m1 * p_m2_m1 +  # "M2 <- M1 -> Y1 -> Y2 -> GY"
              p_gy_y2 * p_y2_m1 * p_m2_m1 +  # "M2 <- M1 -> Y2 -> GY"
              p_gy_y2 * p_y2_u2 * p_u2_u1 * p_m1_u1 * p_m2_m1 +  # "M2 <- M1 <- U1 -> U2 -> Y2 -> GY"
              p_gy_y1 * p_y1_u1 * p_m1_u1 * p_m2_m1 +  # "M2 <- M1 <- U1 -> Y1 -> GY"
              p_gy_y2 * p_y2_y1 * p_y1_u1 * p_m1_u1 * p_m2_m1 +  # "M2 <- M1 <- U1 -> Y1 -> Y2 -> GY"
              p_gy_y2 * p_y2_u2 * p_m2_u2 +  # "M2 <- U2 -> Y2 -> GY"
              p_gy_y1 * p_y1_m1 * p_m1_u1 * p_u2_u1 * p_m2_u2 +  # "M2 <- U2 <- U1 -> M1 -> Y1 -> GY"
              p_gy_y2 * p_y2_y1 * p_y1_m1 * p_m1_u1 * p_u2_u1 * p_m2_u2 +  # "M2 <- U2 <- U1 -> M1 -> Y1 -> Y2 -> GY"
              p_gy_y2 * p_y2_m1 * p_m1_u1 * p_u2_u1 * p_m2_u2 +  # "M2 <- U2 <- U1 -> M1 -> Y2 -> GY"
              p_gy_y1 * p_y1_u1 * p_u2_u1 * p_m2_u2 +  # "M2 <- U2 <- U1 -> Y1 -> GY"
              p_gy_y2 * p_y2_y1 * p_y1_u1 * p_u2_u1 * p_m2_u2 +  # "M2 <- U2 <- U1 -> Y1 -> Y2 -> GY"
              p_gy_y2 * p_y2_x * p_m2_x +  # "M2 <- X -> Y2 -> GY"
              p_gy_y1 * p_m2_y1 +  # "M2 <- Y1 -> GY"
              p_gy_y2 * p_y2_y1 * p_m2_y1 +  # "M2 <- Y1 -> Y2 -> GY"
              p_gy_y2 * p_y2_m1 * p_y1_m1 * p_m2_y1 +  # "M2 <- Y1 <- M1 -> Y2 -> GY"
              p_gy_y2 * p_y2_u2 * p_u2_u1 * p_m1_u1 * p_y1_m1 * p_m2_y1 +  # "M2 <- Y1 <- M1 <- U1 -> U2 -> Y2 -> GY"
              p_gy_y2 * p_y2_m1 * p_m1_u1 * p_y1_u1 * p_m2_y1 +  # "M2 <- Y1 <- U1 -> M1 -> Y2 -> GY"
              p_gy_y2 * p_y2_u2 * p_u2_u1 * p_y1_u1 * p_m2_y1  # "M2 <- Y1 <- U1 -> U2 -> Y2 -> GY"
            ,
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
grid1 <- na.omit(grid1)

grid1[,`:=`( rm2y2_m1y1 = (rm2y2_m1 - ry1m2_m1 * ry1y2_m1) / sqrt((1 - ry1m2_m1^2) * (1 - ry1y2_m1^2)),
             rxy2_m1y1 = (rxy2_m1 - ry1x_m1 * ry1y2_m1) / sqrt((1 - ry1x_m1^2) * (1 - ry1y2_m1^2)),
             rxm2_m1y1 = (rxm2_m1 - ry1x_m1 * ry1m2_m1) / sqrt((1 - ry1x_m1^2) * (1 - ry1m2_m1^2)))]

grid1 <- na.omit(grid1)
grid1[,`:=`(bm2y2_xm1y1 = ((rm2y2_m1y1- rxm2_m1y1 *rxy2_m1y1) /
                             (1 - rxm2_m1y1 ^2 ))* sd_y2_m1y1 / sd_m2_m1y1
)]



grid1 <- na.omit(grid1)

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
                   abs(rxm2_m1y1) > 1)]

grid1[,`:=`(change = p_m2_x*bgygm_x/(2-2*rm1m2_x),
            change_2 = p_m2_x*bgygm_x/(1-rm1m2_x),
            change_1 = p_m2_x*bgym2_x,
            ancova = p_m2_x*bm2y2_xm1y1 ,
            truth = p_m2_x*p_y2_m2,
            naive =p_m2_x*by2m2_x)]
# computing bias
grid1[,`:=`(bias_change = change - truth,
            bias_change_2 = change_2 - truth,
            bias_ancova = ancova - truth,
            bias_naive = naive - truth,
            bias_change_1 = change_1 - truth)]
grid1[1,]

grid1$diff <- abs(grid1$p_m2_u2*grid1$p_y2_u2)

grid1 <- grid1[,-(18:58)]
colnames(grid1)
grid1[,17:28] <- round_df(grid1[,17:28],4)
grid1$Min_1 <- colnames(grid1[,24:28])[minCol(abs(grid1[,24:28]), ties.method = "first")]
grid1$Min_2 <- colnames(grid1[,24:28])[minCol(abs(grid1[,24:28]), ties.method = "last")]
grid1<- grid1 |> mutate(Min = ifelse(Min_1 != Min_2, "Tie", Min_1))
# RESULTS---------------------------------------------------------------------------------------------

#overlaid densities of both biases
graph1 <- grid1$bias_ancova |>
  as_tibble()

graph1_ancova <- graph1 |>  # Group by all columns
  group_by(value, diff) |>
  count()
graph1_ancova <- graph1_ancova |>
  mutate(perc = n/nrow(graph1))
write.csv(graph1_ancova,"ancova_2_1.csv")

graph1 <- grid1$bias_change |>
  as_tibble()

graph1_change <- graph1 |>  # Group by all columns
  group_by(value, diff) |>
  count()
graph1_change <- graph1_change|>
  mutate(perc = n/nrow(graph1))

write.csv(graph1_change,"change_2_1.csv")

graph1 <- grid1$bias_change_2 |>
  as_tibble()

graph1_change2 <- graph1 |>  # Group by all columns
  group_by(value, diff) |>
  count()
graph1_change2 <- graph1_change2|>
  mutate(perc = n/nrow(graph1))

write.csv(graph1_change2,"change2_2_1.csv")

graph1 <- grid1$bias_naive |>
  as_tibble()
graph1_naive <- graph1 |>  # Group by all columns
  group_by(value, diff) |>
  count()
graph1_naive<- graph1_naive |>
  mutate(perc = n/nrow(graph1))

write.csv(graph1_naive,"naive_2_1.csv")

graph1 <- grid1$bias_change_1 |>
  as_tibble()

graph1_change1 <- graph1 |>  # Group by all columns
  group_by(value, diff) |>
  count()
graph1_change1 <- graph1_change1|>
  mutate(perc = n/nrow(graph1))

write.csv(graph1_change1,"change1_2_1.csv")

grid1small1 <- grid1 |> filter(p_y1_u1 + p_m1_u1*p_y1_m1== p_y2_u2*p_u2_u1+p_y1_u1*p_y2_y1+
                                 p_m1_u1*p_y2_m1+p_m1_u1*p_y1_m1*p_y2_y1&
                                 p_y1_m1 + p_y1_u1*p_m1_u1 ==
                                 p_y1_m1*p_y2_y1 + p_y2_m1 + p_m1_u1*p_y1_u1*p_y2_y1+
                                 p_m1_u1*p_u2_u1*p_y2_u2 &
                                 p_m2_y1 == 0)


#overlaid densities of both biases



graph1 <- grid1small1$bias_ancova |>
  as_tibble()

graph1_ancova <- graph1 |>  # Group by all columns
  group_by(value) |>
  count()
graph1_ancova <- graph1_ancova |>
  mutate(perc = n/nrow(graph1))
write.csv(graph1_ancova,"ancova_2_1_c.csv")

graph1 <- grid1small1$bias_change |>
  as_tibble()

graph1_change <- graph1 |>  # Group by all columns
  group_by(value) |>
  count()
graph1_change <- graph1_change|>
  mutate(perc = n/nrow(graph1))

write.csv(graph1_change,"change_2_1_c.csv")

graph1 <- grid1small1$bias_change_2 |>
  as_tibble()

graph1_change2 <- graph1 |>  # Group by all columns
  group_by(value) |>
  count()
graph1_change2 <- graph1_change2|>
  mutate(perc = n/nrow(graph1))

write.csv(graph1_change2,"change2_2_1_c.csv")

graph1 <- grid1small1$bias_naive |>
  as_tibble()
graph1_naive <- graph1 |>  # Group by all columns
  group_by(value) |>
  count()
graph1_naive<- graph1_naive |>
  mutate(perc = n/nrow(graph1))

write.csv(graph1_naive,"naive_2_1_c.csv")

graph1 <- grid1$bias_change_1 |>
  as_tibble()

graph1_change1 <- graph1 |>  # Group by all columns
  group_by(value) |>
  count()
graph1_change1 <- graph1_change1|>
  mutate(perc = n/nrow(graph1))
write.csv(graph1_change1,"change1_2_1_c.csv")

grid1small1 <- grid1 |> filter(grid1$p_m2_u2* grid1$p_y2_u2==0)
graph1 <- grid1small1$bias_ancova |>
  as_tibble()

graph1_ancova <- graph1 |>  # Group by all columns
  group_by(value) |>
  count()
graph1_ancova <- graph1_ancova |>
  mutate(perc = n/nrow(graph1))
write.csv(graph1_ancova,"ancova_2_1_i.csv")

graph1 <- grid1small1$bias_change |>
  as_tibble()

graph1_change <- graph1 |>  # Group by all columns
  group_by(value) |>
  count()
graph1_change <- graph1_change|>
  mutate(perc = n/nrow(graph1))

write.csv(graph1_change,"change_2_1_i.csv")

graph1 <- grid1small1$bias_change_2 |>
  as_tibble()

graph1_change2 <- graph1 |>  # Group by all columns
  group_by(value) |>
  count()
graph1_change2 <- graph1_change2|>
  mutate(perc = n/nrow(graph1))

write.csv(graph1_change2,"change2_2_1_i.csv")

graph1 <- grid1small1$bias_naive |>
  as_tibble()
graph1_naive <- graph1 |>  # Group by all columns
  group_by(value) |>
  count()
graph1_naive<- graph1_naive |>
  mutate(perc = n/nrow(graph1))

write.csv(graph1_naive,"naive_2_1_i.csv")

graph1 <- grid1$bias_change_1 |>
  as_tibble()

graph1_change1 <- graph1 |>  # Group by all columns
  group_by(value) |>
  count()
graph1_change1 <- graph1_change1|>
  mutate(perc = n/nrow(graph1))
write.csv(graph1_change1,"change1_2_1_i.csv")