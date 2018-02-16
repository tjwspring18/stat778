library(survival)

setwd("~/code/stat778/hw1")

r <- read.csv("HW1.dat", sep = " ", header = F)

fit <- survfit(Surv(V1, V2) ~ 1, data = r, conf.type = "plain")

plot(fit)

tjw <- read.csv("blargh", sep = ",", header = T)

points(x=tjw$time, y=tjw$survival, col = "red")
points(x=tjw$time, y=tjw$lower95, col = "red")
points(x=tjw$time, y=tjw$upper95, col = "red")
