# the longitudinal mediation model

# Load required packages.
library(profvis)
library(tictoc)
library(plyr)
library(doParallel)
library(dagitty)
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

# Setting the number of cores to speed up the analysis
cl<-49
registerDoParallel(cl)

setDTthreads(percent = 90)
#### analytic bias from u1, u2
#### 9 grid points will only run on supercomputer
#### lower number of gridpoints yield similar results

profvis({

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
p_y2_m2 <- ng 
p_m2_m1 <- ng 
p_y2_y1 <- ng 
p_y2_u2 <- ng 
p_m2_u2 <- ng 
p_u2_u1 <- ng
p_m1_u1 <- ng
p_y1_u1 <- ng 
p_y1_m1 <- ng 
p_gm_m2 <- 1
p_gm_m1 <- -1
p_gy_y2 <- 1
p_gy_y1 <- -1


grid1 <- data_frame()
grid <- data.table()
for (i in 1:length(ng)){
  p_y2_x <- ng[[i]]
  grid <- CJ(p_m2_x, p_y2_x, p_y2_m2, p_m2_m1, p_y2_y1, p_y2_u2, p_m2_u2, p_y1_u1,
              p_m1_u1, p_u2_u1, p_y1_m1, p_gm_m2, p_gm_m1, p_gy_y2, p_gy_y1)
  grid[,`:=`(eps_m1 = 1-p_m1_u1^2,
              eps_y1 = 1-p_y1_u1^2-
                              p_y1_m1^2-
                              2*p_y1_u1*p_m1_u1*p_y1_m1,
              eps_m2 = 1-p_m2_m1^2-
                              p_m2_u2^2-
                              p_m2_x^2-
                              2*p_m2_m1*p_m2_u2*p_m1_u1*p_u2_u1,
              eps_y2 = 1-p_y2_y1^2-
                              p_y2_m2^2-
                              p_y2_u2^2-
                              p_y2_x^2-
                              2*p_y2_x*p_m2_x*p_y2_x-
                              2*p_y2_m2*p_m2_m1*(p_y1_m1*p_y2_y1+p_m1_u1*p_y2_u2*p_u2_u1+p_m1_u1*p_y1_u1*p_y2_y1)-
                              2*p_y2_u2*p_y1_u1*p_u2_u1*p_y2_y1)]
  grid <- grid[eps_m1 >= 0 & eps_m2 >= 0 & eps_y1 >= 0 & eps_y2 >= 0]
  grid1 <- rbind(grid1,grid)
}
  

            

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

grid1 <- grid1[sd_m2.m1 >= 0 & sd_y2.m1 >= 0 & sd_y1.m1 >= 0]

grid1 <- na.omit(grid1)
summarize(grid1)
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
sum(is.na(grid1))
grid1[,16:20] <- round_df(grid1[,16:20], 10)
grid1$Min_1 <- colnames(grid1[,16:19])[minCol(abs(grid1[,16:19]), ties.method = "first")]
grid1$Min_2 <- colnames(grid1[,16:19])[minCol(abs(grid1[,16:19]), ties.method = "last")]
grid1<- grid1 |> mutate(Min = ifelse(Min_1 != Min_2, "Tie", Min_1))


# RESULTS---------------------------------------------------------------------------------------------
# Final number of |Bias_naive| >= |Bias_fd|
# no constraint
(f.tab <- with(grid1, table(Min)))
prop.table(f.tab)

# p_m1_u == p_m2_u
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1,], 
               table(Min)))
prop.table(f.tab)

# p_y1_u == p_y2_u
(f.tab <- with(grid1[grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1,], 
               table(Min)))
prop.table(f.tab)

# p_y1_m1 == 0
(f.tab <- with(grid1[grid1$p_y1_m1 == 0,], 
               table(NearMin)))
prop.table(f.tab)

# p_m2_m1 == 0
(f.tab <- with(grid1[grid1$p_m2_m1 == 0,], 
               table(NearMin)))
prop.table(f.tab)

# p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_y2_y1 == 0,], 
               table(NearMin)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_u == p_y2_u
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1&
                       grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1 ,], 
               table(NearMin)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_m1==0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1 &
                       grid1$p_y1_m1 == 0,], 
               table(NearMin)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_m2_m1==0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1 &
                       grid1$p_m2_m1 == 0,], 
               table(NearMin)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1&
                       grid1$p_y2_y1 == 0,], 
               table(NearMin)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_y1_m1 == 0
(f.tab <- with(grid1[grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1&
                       grid1$p_y1_m1 == 0,], 
               table(NearMin)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_m2_m1 == 0
(f.tab <- with(grid1[grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1 &
                       grid1$p_m2_m1 == 0,], 
               table(NearMin)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1&
                       grid1$p_y2_y1 == 0,], 
               table(NearMin)))
prop.table(f.tab)

# p_y1_m1 == 0 & p_m2_m1 == 0
(f.tab <- with(grid1[grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0,], 
               table(NearMin)))
prop.table(f.tab)

# p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(NearMin)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_u == p_y2_u & p_y1_m1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1 &
                       grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1 &
                       grid1$p_y1_m1 == 0,], 
               table(NearMin)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_u == p_y2_u & p_m2_m1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1 &
                       grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1 &
                       grid1$p_m2_m1 == 0,], 
               table(NearMin)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_u == p_y2_u & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1 &
                       grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1 &
                       grid1$p_y2_y1 == 0,], 
               table(NearMin)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_m1 == 0 & p_m2_m1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1&
                       grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0,], 
               table(NearMin)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1&
                       grid1$p_y1_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(NearMin)))
prop.table(f.tab)


# p_m1_u == p_m2_u & p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1 &
                       grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(NearMin)))
prop.table(f.tab)

# p_y1_m1 == 0 & p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(NearMin)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_y1_m1 == 0 & p_m2_m1 == 0
(f.tab <- with(grid1[grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1 &
                       grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0,], 
               table(NearMin)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_y1_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1 &
                       grid1$p_y1_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(NearMin)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1 &
                       grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(NearMin)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_m1 == 0 & p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1&
                       grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(NearMin)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_y1_m1 == 0 & p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1 &
                       grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(NearMin)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_m1_u == p_m2_u & p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1 &
                       grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1 &
                       grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(NearMin)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_m1_u == p_m2_u & p_y1_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1 &
                       grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1 &
                       grid1$p_y1_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(NearMin)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_m1_u == p_m2_u & p_y1_m1 == 0 & p_m2_m1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1 &
                       grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1 &
                       grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0,], 
               table(NearMin)))
prop.table(f.tab)

# common trend assumption
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1 &
                       grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1 &
                       grid1$p_y1_m1 == 0,], 
               table(NearMin)))
prop.table(f.tab)

# common trend assumption
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1 &
                       grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1&
                       grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(NearMin)))
prop.table(f.tab)


####### graphs
#overlaid densities of both biases
graph1 <- pivot_longer(grid1, cols = c(bias_change, bias_naive, bias_ancova),
                       names_to = "bias_type", values_to = "bias_result")
p11 <- ggplot(graph1,aes(x=(bias_result),y=after_stat(count)/max(count), line_type = bias_type)) + 
  geom_density(bw=.13,trim=FALSE) +
  labs(x="Raw bias", y="Percent") + 
  xlim(-1.7,1.7)+
  theme_apa() + 
  theme(legend.position = c(.2,.8))

#tikz("biasdensity.tex",width = 4,height = 3, standAlone = TRUE)
pdf("figure_p11_2.pdf", height=6, width=6)
p11
dev.off()

graph2 <- pivot_longer(grid1, cols = c(bias_change_2, bias_naive, bias_ancova),
                       names_to = "bias_type", values_to = "bias_result")
p12 <- ggplot(graph2,aes(x=bias_result, y=..count../max(count), linetype = bias_type)) + 
  geom_density(bw=.13,trim=FALSE) + 
  labs(x="Raw bias", y="Percent") + 
  xlim(-1.7,1.7)+
  theme_apa() + theme(legend.position = c(.2,.8))
#tikz("biasdensity_2.tex",width = 4,height = 3, standAlone = TRUE)
pdf("figure_p12_2.pdf", height=6, width=6)
p12
dev.off()

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
graph3 <- pivot_longer(grid1small1, cols = c(bias_change_2, bias_naive, bias_ancova),
                       names_to = "bias_type", values_to = "bias_result")
p2 <- ggplot(graph3,aes(x=bias_result,y=..count../max(count), linetype = bias_type)) + 
  geom_histogram() + facet_grid(~bias_type)
labs(x="Raw bias", y="Percent") + 
  xlim(-1.7,1.7) + 
  theme_apa() + theme(legend.position = c(.2,.8))
#tikz("biasdensitynogamma.tex",width = 4,height = 3, standAlone = TRUE)
pdf("figure_p2_2.pdf", height=6, width=6)
p2
dev.off()
table(grid1small1$bias_change_1)
ptable2 <- prop.table(table(grid1small1$bias_ind_2))
ptable2b <- prop.table(table(grid1small1$bias_ind_eq_2))

})


stopCluster(cl)
