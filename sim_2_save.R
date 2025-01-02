library(profvis)
library(tictoc)
library(plyr)
library(doParallel)
#cl<-49
#registerDoParallel(cl)
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

for (i in length(ng)){
  p_y2_x <- ng[[i]]
  grid <- CJ(p_m2_x, p_y2_x, p_y2_m2, p_m2_m1, p_y2_y1, p_y2_u2, p_m2_u2, p_y1_u1,
             p_m1_u1, p_u2_u1, p_y1_m1, p_gm_m2, p_gm_m1, p_gy_y2, p_gy_y1)
  grid[,`:=`(eps_m1 = sqrt(1-p_m1_u1^2),
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
  grid <- grid[eps_m1>=0 & eps_y1>=0 & eps_m2>=0 & eps_y2>=0]
  grid1 <- rbind(grid1,grid)
}
grid1[1:10]
rm(grid)
grid1 <- na.omit(grid1)
# computing estimators
grid1[,`:=`(rm2y2 = p_y2_m2 +
              p_m2_m1 * p_y1_m1 * p_y2_y1 +
              p_m2_m1 * p_m1_u1 * p_u2_u1 * p_y2_u2 +
              p_m2_m1 * p_m1_u1 * p_y1_u1 * p_y2_y1 +
              p_m2_u2 * p_y2_u2 +
              p_m2_u2 * p_u2_u1 * p_m1_u1 * p_y1_m1 * p_y2_u2 +
              p_m2_u2 * p_u2_u1 * p_y1_u1 * p_y2_y1, 
            rm1m2 = p_m2_m1 + 
              p_m1_u1 * p_u2_u1 * p_m2_u2,
            rm1y2 = p_m2_m1*p_y2_m2 + 
              p_m1_u1*p_u2_u1*p_y2_u2 + 
              p_y1_m1*p_y2_y1 + 
              p_m1_u1*p_y1_u1*p_y2_y1 +
              p_m1_u1*p_u2_u1*p_m2_u2*p_y2_m2,
            ry1m2 = p_m2_m1*p_y1_m1 +
              p_m2_m1*p_m1_u1*p_y1_u1+
              p_y1_m1*p_m1_u1*p_u2_u1*p_m2_u2 +
              p_y1_u1*p_u2_u1*p_m2_u2,
            ry1y2 = p_y2_y1 +
              p_y1_m1*p_m2_m1*p_y2_m2 + 
              p_y1_m1*p_m1_u1*p_u2_u1*p_m2_u2*p_y2_m2 +
              p_y1_m1*p_m1_u1*p_u2_u1*p_y2_u2 +
              p_y1_u1*p_m1_u1*p_m2_m1*p_y2_m2 + 
              p_y1_u1*p_u2_u1*p_y2_u2 + 
              p_y1_u1*p_u2_u1*p_m2_u2*p_y2_m2,
            rm1y1 = p_y1_m1 + 
              p_m1_u1*p_y1_u1)]

grid1[,`:=`(ry1m2.m1 = (ry1m2-rm1y1*rm1m2)/sqrt((1-rm1y1^2)*(1-rm1m2^2)),
            ry1y2.m1 = (ry1y2-rm1y1*rm1y2)/sqrt((1-rm1y1^2)*(1-rm1y2^2)),
            rm2y2.m1 = (rm2y2-rm1m2*rm1y2)/sqrt((1-rm1m2^2)*(1-rm1y2^2)),
            sd_m2.m1 = sqrt(1-rm1m2^2),
            sd_y2.m1 = sqrt(1-rm1y2^2),
            sd_y1.m1 = sqrt(1-rm1y1^2))]

grid1 <- na.omit(grid1)

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
            ancova = p_m2_x*(((rm2y2.m1*sd_y2.m1/sd_m2.m1)-(ry1m2.m1*sd_y1.m1/sd_m2.m1)*(ry1y2.m1*sd_y2.m1/sd_y1.m1))/
                               (1-((ry1m2.m1*sd_y1.m1/sd_m2.m1)^2*sd_m2.m1^2/sd_y1.m1^2))),
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

grid1 <- grid1[,-(16:37)]

grid1[,16:20] <- round_df(grid1[,16:20], 4)
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
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1&
                       grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1 ,], 
               table(Min)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_m1==0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1 &
                       grid1$p_y1_m1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_m2_m1==0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1 &
                       grid1$p_m2_m1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1&
                       grid1$p_y2_y1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_y1_m1 == 0
(f.tab <- with(grid1[grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1&
                       grid1$p_y1_m1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_m2_m1 == 0
(f.tab <- with(grid1[grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1 &
                       grid1$p_m2_m1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1&
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
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1 &
                       grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1 &
                       grid1$p_y1_m1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_u == p_y2_u & p_m2_m1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1 &
                       grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1 &
                       grid1$p_m2_m1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_u == p_y2_u & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1 &
                       grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1 &
                       grid1$p_y2_y1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_m1 == 0 & p_m2_m1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1&
                       grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1&
                       grid1$p_y1_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(Min)))
prop.table(f.tab)


# p_m1_u == p_m2_u & p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1 &
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
(f.tab <- with(grid1[grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1 &
                       grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_y1_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1 &
                       grid1$p_y1_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1 &
                       grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_m1_u == p_m2_u & p_y1_m1 == 0 & p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1&
                       grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_y1_m1 == 0 & p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1 &
                       grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_m1_u == p_m2_u & p_m2_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1 &
                       grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1 &
                       grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_m1_u == p_m2_u & p_y1_m1 == 0 & p_y2_y1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1 &
                       grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1 &
                       grid1$p_y1_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(Min)))
prop.table(f.tab)

# p_y1_u == p_y2_u & p_m1_u == p_m2_u & p_y1_m1 == 0 & p_m2_m1 == 0
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1 &
                       grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1 &
                       grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0,], 
               table(Min)))
prop.table(f.tab)

# common trend assumption
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1 &
                       grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1 &
                       grid1$p_y1_m1 == 0,], 
               table(Min)))
prop.table(f.tab)

# common trend assumption
(f.tab <- with(grid1[grid1$p_m1_u1 == grid1$p_u2_u1*grid1$p_m2_u2 + grid1$p_m1_u1*grid1$p_m2_m1 &
                       grid1$p_y1_u1 == grid1$p_u2_u1*grid1$p_y2_u2 + grid1$p_y1_u1*grid1$p_y2_y1&
                       grid1$p_y1_m1 == 0 &
                       grid1$p_m2_m1 == 0 &
                       grid1$p_y2_y1 == 0,], 
               table(Min)))
prop.table(f.tab)

#####graph

####### graphs
#overlaid densities of both biases
graph1 <- grid1$bias_ancova |>
  as_tibble()

graph1_ancova <- graph1 |>  # Group by all columns
  group_by(value) |>
  count()
graph1_ancova <- graph1_ancova |>
  mutate(perc = n/nrow(graph1))
write.csv(graph1_ancova,"ancova_1.csv")

graph1 <- grid1$bias_change |>
  as_tibble()

graph1_change <- graph1 |>  # Group by all columns
  group_by(value) |>
  count()
graph1_change <- graph1_change|>
  mutate(perc = n/nrow(graph1))

write.csv(graph1_change,"change_1.csv")

graph1 <- grid1$bias_change_2 |>
  as_tibble()

graph1_change2 <- graph1 |>  # Group by all columns
  group_by(value) |>
  count()
graph1_change2 <- graph1_change2|>
  mutate(perc = n/nrow(graph1))

write.csv(graph1_change2,"change2_1.csv")

graph1 <- grid1$bias_naive |>
  as_tibble()
graph1_naive <- graph1 |>  # Group by all columns
  group_by(value) |>
  count()
graph1_naive<- graph1_naive |>
  mutate(perc = n/nrow(graph1))

write.csv(graph1_naive,"naive_1.csv")

grid1small1 <- grid1 |> filter(p_m1_u1 == p_u2_u1*p_m2_u2 + p_m1_u1*p_m2_m1 &
                                 p_y1_u1 == p_u2_u1*p_y2_u2 + p_y1_u1*p_y2_y1&
                                 p_y1_m1 == 0 &
                                 p_m2_m1 == 0 &
                                 p_y2_y1 == 0)


#overlaid densities of both biases
graph1 <- grid1small1$bias_ancova |>
  as_tibble()

graph1_ancova <- graph1 |>  # Group by all columns
  group_by(value) |>
  count()
graph1_ancova <- graph1_ancova |>
  mutate(perc = n/nrow(graph1))
write.csv(graph1_ancova,"ancova_2_1.csv")

graph1 <- grid1small1$bias_change |>
  as_tibble()

graph1_change <- graph1 |>  # Group by all columns
  group_by(value) |>
  count()
graph1_change <- graph1_change|>
  mutate(perc = n/nrow(graph1))

write.csv(graph1_change,"change_2_1.csv")

graph1 <- grid1small1$bias_change_2 |>
  as_tibble()

graph1_change2 <- graph1 |>  # Group by all columns
  group_by(value) |>
  count()
graph1_change2 <- graph1_change2|>
  mutate(perc = n/nrow(graph1))

write.csv(graph1_change2,"change2_2_1.csv")

graph1 <- grid1small1$bias_naive |>
  as_tibble()
graph1_naive <- graph1 |>  # Group by all columns
  group_by(value) |>
  count()
graph1_naive<- graph1_naive |>
  mutate(perc = n/nrow(graph1))

write.csv(graph1_naive,"naive_2_1.csv")
