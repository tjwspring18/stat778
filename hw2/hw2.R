setwd("~/code/stat778/hw2/")
df <- read.csv("hw2_2.csv")

a <- c()
for(i in 1:1000){
	d <- rnorm(50, -0.5, sqrt(2))
	a <- c(a, mean(d))
}
print(1-((length(a[a > -0.097]) + length(a[a < -0.877]))/1000))

a <- c()
for(i in 1:1000){
	d <- rnorm(50, -0.5, sqrt(2))
	a <- c(a, var(d))
}
print(1-((length(a[a > 2.79]) + length(a[a < 1.207]))/1000))

a <- c()
for(i in 1:1000){
	d <- rnorm(100, -0.5, sqrt(2))
	a <- c(a, mean(d))
}
print(1-((length(a[a > -0.220]) + length(a[a < -0.771]))/1000))

a <- c()
for(i in 1:1000){
	d <- rnorm(100, -0.5, sqrt(2))
	a <- c(a, var(d))
}
print(1-((length(a[a < 1.453]) + length(a[a > 2.544]))/1000))

a <- c()
for(i in 1:1000){
	d <- rnorm(200, -0.5, sqrt(2))
	a <- c(a, mean(d))
}
print(1-((length(a[a > -0.302]) + length(a[a < -0.693]))/1000))

a <- c()
for(i in 1:1000){
	d <- rnorm(100, -0.5, sqrt(2))
	a <- c(a, var(d))
}
print(1-((length(a[a < 1.600]) + length(a[a > 2.382]))/1000))
