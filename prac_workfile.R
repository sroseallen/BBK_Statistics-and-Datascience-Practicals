library(tidyr)
library(readr)
library(dplyr)
library(MASS)

if(!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install(c("nortest","moments","hopach","ape"))
library(nortest)
library(moments)
library(hopach)
library(ape)

# SESSION 1 ----
write(1:10, file="test.txt")
getwd()
help.start()
??read

3+5
123/135.6
2^9
pi
4*pi*2^3/3

v <- c(2,7,3,10,-1)
w <- 1:5
u <- v*w
u
exp(v)
log(v)
exp(log(v))
x <- seq(0,10,0.1) #0 to 10 with 0.1 as an interval.
x[2]
v[v>3]
sum(v>3) #how many instances of digits > 3
which(v>3) #what position in the vector are digits > 3

v <- 1:9
w <- 101:109
u <- c(v,w)
u
v <- NULL
rm (x)
ls()
rm(list=ls())

## General plotting ----
x <- seq(1,10,0.4)
y <- (x^3 - 4*x^2 - 35*x)
plot(x,y)
plot(x,y,type="o", col="blue", ylim=c(-200,300))
lines(c(0,10),c(0,0),col="black")
lines(x,(x^2 + 4*x + 9), col="red", lty=2)
lines(x, (-2*x^2 + 4*x + 14), col="green", lty=3)
title(main="polynomial plot", col.main="black", font.main=14)

interface_data <- read.table("../Resources_practicals/protein_interface_data.txt")
interface_data
names(interface_data) <- c("PDBcode", "protein_type", "area")
areas <- interface_data$area

## Boxplots ----
mean(areas)
median(areas)
quantile(areas)
box <- boxplot(quantile(areas))
original_par <- par()
par(pin=c(3,1), font=2, ps=10, family="sans")
box <- boxplot(areas, horizontal=TRUE, xlab="crystal interface areas")
par(original_part)

areas_dimer <- interface_data$area[interface_data$protein_type == "d"]
areas_monomer <- interface_data$area[interface_data$protein_type == "m"]
par(mfrow=c(1,2))
boxplot(areas_dimer, 
        ylab="monomer crystal interface areas", 
        ylim = c(0, 7000))
boxplot(areas_monomer, 
        ylab="dimer crystal interface areas", 
        ylim = c(0, 7000))

summary(areas_dimer)
summary(areas_monomer)

## Histograms ----
par(mfrow=c(1,1))
hist(areas_dimer)
hist(areas_dimer, breaks=10)

# force a specific number of bins
my.breaks <- seq(0, 8000, by=300)
hist(areas_dimer, breaks=my.breaks)

# pretty histograms :)
par(mfrow=c(1,2))
hist(areas_monomer, 
     breaks=seq(0, 2000, 400),
     xlab="Interface area (Angstroms)",
     ylab="Frequency",
     main="Monomers",
     xlim=c(0,8000))
hist(areas_dimer, 
     breaks=seq(0,8000,400),
     xlab="Interface area (Angstroms)",
     ylab="Frequency",
     main="Dimers",
     xlim=c(0,8000))

# exporting a plot
jpeg('../Resources_practicals/rplot.jpg')
par(mfrow=c(1,2))
boxplot(areas_dimer)
boxplot(areas_monomer)
dev.off()



# SESSION 2 ----

## Bionomial Distribution (Discrete) ----
choose(12,10) # computes the binomial coefficient - prob of binary outcome occuring exactly a k number of times out of n times [coefficient given as: n! / k! (n-k)!)]
sum(dbinom(0:0,12,0.8))  # k = 0
sum(dbinom(0:1,12,0.8))  # k <= 1
sum(dbinom(0:2,12,0.8))  # k <= 2
sum(dbinom(0:10,12,0.8)) # k <= 10
sum(dbinom(0:12,12,0.8)) # k <= 12

barplot(dbinom(0:12,size=12,prob=0.8)) # probability density binomial distribution
barplot(pbinom(0:12,size=12,prob=0.8),names.arg=c(0:12),cex.names=0.8,ylab="probability", xlab="Recovered") # cumulative binomial distribution - sum of all outcomes up to k number of times
qbinom() #quantiles; The value below which a specified fraction of the distribution lies
rbinom() #random

## Poisson Distribution (Discrete) ----
## The mean and variance of a random variable with a Poisson distribution are the same and equal to λ (i.e. s2=l and m=l)
barplot(dpois(1:20,10),names.arg=cat_name,cex.names=0.5,xlab="no. of events", ylab="probability")
p <- 2.5/1000000
n <- 2*1000000 # base pairs
lambda <- p*n # equals 5

dpois(10,lambda) # probability density poisson distribution
ppois(10,lambda) # cumulative poisson distribution: prob of occurring 10 times or fewer
ppois(10,lambda,lower=FALSE) # as above: prob of more than 10 times
qpois() #quantiles; The value below which a specified fraction of the distribution lies
rpois() #random

## Uniform Distribution ----
dunif() #probability density of getting value close to x (continuous distribution) or probability of getting exactly the value of x (discrete distribution)
punif() #probability of a value of x or less - P(X ≤ x)
qunif() #quantiles; The value below which a specified fraction of the distribution lies
runif() #random

## Practical 2 ----
barplot(dbinom(0:12,size=12,prob=0.8),names.arg=c(0:12),cex.names=0.8,ylab="probability", xlab="Recovered")

### 2.1 ----
# a
dbinom(14,22,0.7)
# b
barplot(dbinom(0:22, size=22, prob=0.7),names.arg=c(0:22), cex.names=0.7, ylab="probability",xlab="purine")
# c
barplot(pbinom(0:22, size=22, prob=0.7),names.arg=c(0:22), cex.names=0.7, ylab="probability",xlab="purine")
# d
pbinom(13,22,0.7)
# e
1 - pbinom(10, 22, 0.7)
#alt e
pbinom(10,size=22,prob=0.7,lower.tail=FALSE)
# f - NOTE THAT YOU NEED TO USE qbinom() TO ACHIEVE THE CORRECT QUANTILE!!
binom_dist <- pbinom(0:22, size=22, prob=0.7)
purines <- c(0:22)
plot(purines, binom_dist, xlab="purines", ylab="probability", type="l")

abline(v=qbinom(p=0.5, size=22, prob=0.7), col="red")

### 2.2 ----
# a
dbinom(1,4,0.25)
# b: 3 or 4 children
sum(dbinom(3:4,4,0.25))
1-pbinom(2,4,0.25)
#c
which.max(dbinom(0:100,100,0.25))
barplot(dbinom(0:100,100,0.25))

### 2.3 ----
# a
p <- 1/200000 #frequency
n <- 45000000
lambda <- p*n # mean number of errors in 45 000 000 n
# b (for poisson, only need to plot around 225, as variance = mean = 225 in this case)
plot(dpois(0:550,lambda))
# calculating median
abline(v=qpois(0.5, lambda))
# could also do a weighted sum to estimate the mean:
weighted_sum<-0
for(i in 0:999) {weighted_sum = weighted_sum + i*dpois(i,lambda=225)}
abline(v=weighted_sum)
# c
plot(ppois(0:550,lambda))
quantile_95 <- qpois(0.95,lambda)
abline(v=quantile_95)
# d
par(mfrow=c(1,2))
p <- 0.6/200000 #frequency
n <- 45000000
lambda <- p*n # mean number of errors in 45 000 000 n
plot(dpois(0:550,lambda))
plot(ppois(0:550,lambda))
quantile_95 <- qpois(p=0.95,lambda)
abline(h=quantile_95)

ppois(250,135) #output of 1 = certainty that there are less than 250 errors with a lambda distribution of 135.

### 2.4 ----
mybreaks <- seq(0,1,0.05)
uniform_data <- runif(10,0,1)
hist(uniform_data, breaks=mybreaks)

for (i in 1:10) {
  uniform_data <- runif(10*(2^i),0,1)
  hist(uniform_data, breaks=mybreaks)
  Sys.sleep(1)
}

samples <- numeric()
mean_values <- numeric()
for (i in 1:10) {
  uniform_data <- runif(10*(2*i),0,1)
  samples[i] <- 10*(2^i)
  mean_values[i] <- mean(uniform_data)
  hist(uniform_data,breaks=mybreaks)
  Sys.sleep(1)
}

plot(samples, mean_values, type="o")

runif(100,-1,1)
sources <- 1
total_error <- rep(0,100)
for (i in 1:sources) {
  total_error = total_error + runif(100,-1,1)
}
hist(total_error,breaks=20)



# SESSION 3 ----
## Normal Distribution ----
x <- seq(-2,12,0.1)
dnorm(6,mean=5,sd=2) # pdf: probability density (pdf) of a given value x (6)
plot(x, dnorm(x,mean=5,sd=2), ty="l", xlab = "x",ylab = "N(x,5,2)") #pdf plot

pnorm(6,mean=5,sd=2) # cdf: cumulative probability (cdf) that a random value x is less than or equal to a value (6)
plot(x, pnorm(x,mean=5,sd=2), ty="l", xlab = "x",ylab = "cdf") #cdf plot

pop_mean=5
pop_stdev=2
sample_obs=8
sample_mean=10.3
1-pnorm(sample_mean,pop_mean,pop_stdev/sqrt(sample_obs)) # how to calculate likelihood of your sample observations coming from the population distribution

qnorm() #quantiles; The value below which a specified fraction of the distribution lies
rnorm() #random

## t-distributions ----
## Essentially: Normal distributions (hypothetical) made when you don't know the population st.dev. and have to estimate st.dev. using the sample st.dev.
dt(x,df) #df = degrees of freedom, not dataframe
pt(x,df)
rt(x,df)
qt(x,df)

## One sample t-test ----
## defaults to two-sided test
t.test(c(1,2,3,4,5),mu=10)

#output: gives df (which is n-1), t-score, and p-value (to reject null if <0.05)
#gives the mean and 95% confidence intervals within which likely contains the true mean for your given sample

## Practical ----
### 3.1 ----
single_sample <- runif(100000,-1,1)
sd_single <- sd(single_sample)

sources <- 30 
nrepeats <- 100000 
total_error<-rep(0,nrepeats) # vector of 100k zeros 
for (i in 1:sources) { 
  # add 100k random values from a uniform distribution 
  total_error = total_error + runif(nrepeats,-1,1) 
} 
hist(total_error,breaks=40) 

# Laplace's formula: Proof
sd_sum <- sd(total_error)

print(sd_sum)
print(sd_single * sqrt(30))

### 3.2 ----
norm_samples <- seq(10,10000,by=10)
norm_means <- NULL

for (i in norm_samples) {
  norm_means <- c(norm_means, mean(rnorm(i,mean=100,sd=5)))
}
plot(norm_samples, norm_means)

### 3.3 ----
sample_num <- 2000

vec_10 <- rep(0,sample_num)
vec_100 <- rep(0,sample_num)
vec_1000 <- rep(0, sample_num)

for (i in 1:sample_num) {
  vec_10[i] <- mean(rnorm(10,100,5))
  vec_100[i] <- mean(rnorm(100,100,5))
  vec_1000[i] <- mean(rnorm(1000,100,5))
}

hist(vec_10,breaks=20,prob=TRUE,xlab="mean values",main="2000 samples, sample size 10") 
curve(dnorm(x,100,5/sqrt(10)), add=TRUE) #note: for theoretical dist, you need to calculate the SAMPLE SD, not use the POPULATION SD.

#alt method for adding a line
x<-seq(90,110,0.1)
lines(x,dnorm(x,100,5/sqrt(10)),col="red")

hist(vec_100,breaks=20,prob=TRUE,xlab="mean values",main="2000 samples, sample size 100") 
curve(dnorm(x,100,5/sqrt(100)), add=TRUE)

hist(vec_1000,breaks=20,prob=TRUE,xlab="mean values",main="2000 samples, sample size 1000") 
curve(dnorm(x,100,5/sqrt(1000)), add=TRUE)

### 3.4: Z-testing ----
sample_10 <- rnorm(10,20,3)
sample_100 <- rnorm(100,20,3)
sample_1000 <- rnorm(1000,20,3)

# H0: These samples are not drawn from a normal distribution with mean=18, sd=3
# HA: These samples are drawn from a normal distribution with mean=18, sd=3

z_10 <- (mean(sample_10)-18) / (3/sqrt(10))
z_100 <- (mean(sample_100)-18) / (3/sqrt(100)) #above critical Z with dist(0,1), therefore HA
z_1000 <- (mean(sample_1000)-18) / (3/sqrt(1000)) #above critical Z with dist(0,1), therefore HA

critical_z_pos <- qnorm(1-0.05/2,0,1)
critical_z_neg <- qnorm(0.05/2,0,1)

### 3.5: t-testing ----
t_10 <- (mean(sample_10)-18) / (sd(sample_10)/sqrt(10))
t_100 <- (mean(sample_100)-18) / (sd(sample_100)/sqrt(10))
t_1000 <- (mean(sample_1000)-18) / (sd(sample_1000)/sqrt(10))

critical_t_10 <- qt(1-0.05/2,9)
critical_t_100 <- qt(1-0.05/2,99)
critical_t_1000 <- qt(1-0.05/2,999)

# dont need to do the t-test manually like above!
t.test(sample_10,mu=18) 
# if p-value is less than the sig level, reject the null. 
# The null is that samples ARE DRAWN FROM A NORMAL DISTRIBUTION.
t.test(sample_100,mu=18)
t.test(sample_1000,mu=18)

### 3.6: dataframes ----
age <- c(60,43,72,35,47)
BMI <- c(28, 32, 21, 27, 35) #body mass index
bp <- c(124, 145, 127, 133, 140) #blood pressure
gender <- c("male", "female", "female", "male", "female") 
medical = cbind(age, BMI, bp) 
dim(medical)
mode(medical)
class(medical)
medical_df <- data.frame(medical, gender)

#base
medical_df[medical_df$"gender"=="female","age"]
#dplyr
medical_df %>% select(age) %>% filter(gender == "female")

#base
medical_df[medical_df$"age">43,"gender"]
#dplyr
medical_df %>% filter(age > 43) %>% select(gender)

#base
medical_df[medical_df$"gender"=="female",]
#dplyr
medical_df %>% filter(gender=="female")

# mean BMI for each factor in the col 'gender'
tapply(medical_df$BMI, medical_df$gender, mean)

### 3.7 ----
#### a ----
# H0 is the Pima population is not more obese than the general female population - i.e. not more significantly different

#### b ----
x=seq(0,50,0.1)
plot(x,dnorm(x,25,5),ty="l")

#### c ----
library(MASS)
pima <- Pima.tr
sample_mean <- mean(pima$bmi) 

#### d ----
sample_sd <- sd(pima$bmi)
# plot sample distribution if null is true as well
plot(x,dnorm(x,25,5/sqrt(200)),type="l",col="red")
lines(x,dnorm(x,sample_mean,sample_sd),ty="l")

#### e ----
t_pima <- (sample_mean-25) / (sample_sd/sqrt(200))
p_pima <- pt(t_pima, 199,lower.tail=FALSE)
critical_t_pima <- qt(1-0.05,199) #one-sided to negative skew (more positive value), therefore do not divide by 2

#or, simply:
t.test(Pima.tr$bmi,mu=25,alternative="greater")

# SESSION 4 ----
## Shapiro-Wilk ----
x <- c(0.8,2,3.7,4,4.8,5.4,6,7.9,8.1,8.6,9,10.1,12,14.9,17,19.2,26.3)
shapiro.test(x)

## Kolmogorov-Smirnov ----
x <- c(0.8,2,3.7,4,4.8,5.4,6,7.9,8.1,8.6,9,10.1,12,14.9,17,19.2,26.3)
ks.test(x,"pnorm",9.4,6.71) #where 9.4 is the mean of the population normal distribution, and 6.71 the standard deviation.
## two-sample KS test
ks.test(x,y)

## Two sample t-test ----
t.test(x,y,alternative="greater")

## Mann-Witney-Wilcoxon test ----
wilcox.test(x,y,alternative="greater") #one-sided

## Practical ----
### 4.1 ----
x <- seq(0,2,0.05)
plot(x,dbeta(x,3,3),type="l")
sample <- rbeta(1000,3,3)
skewness(sample)
kurtosis(sample)

beta_10 <- rbeta(10,2,2)
x <- seq(0,1,0.05)
hist(beta_10,breaks=5,prob=TRUE,main="10 sample beta plot")
curve(dbeta(x,2,2),add=TRUE)
#alt method: lines(x,dbeta(x,2,2),type="l")

shapiro.test(beta_10) # if not normal: likely lots of extreme values (large skew)
ks.test(beta_10,"pnorm",mean(beta_10),sd(beta_10)) # if not normal: likely not many values around the mean (large kurtosis)
#PLEASE NOTE THAT YOU NEED TO SPECIFY THE SAMPLE MEAN AND SD OR KS WILL NOT WORK!!!
ad.test(beta_10) # similar conclusions as the S-W test.

### 4.2 ----
beta_50 <- rbeta(50,2,2)
hist(beta_50,breaks=15,prob=TRUE,main="Beta plot")
shapiro.test(beta_50) 
ks.test(beta_50,"pnorm",mean(beta_50),sd(beta_50))
ad.test(beta_50)

### shapiro-wilks power
beta_sw_truepos <- 0
for (i in 1:1000) {
  beta_40sample <- rbeta(40,3,2)
  p_value <- as.numeric(shapiro.test(beta_40sample)[2])
  if (p_value < 0.05) {
    beta_sw_truepos = beta_sw_truepos + 1
  }
}
sw_power <- beta_sw_truepos/1000

### Kolmogorov-Smirnov power
beta_ks_truepos <- 0
for (i in 1:1000) {
  beta_40sample <- rbeta(40,3,2)
  p_value <- as.numeric(ks.test(beta_40sample,"pnorm",mean(beta_40sample),sd(beta_40sample))[2])
  if (p_value < 0.05) {
    beta_ks_truepos = beta_ks_truepos + 1
  }
}
ks_power <- beta_ks_truepos/1000

### Anderson-Darling power
beta_ad_truepos <- 0
for (i in 1:1000) {
  beta_40sample <- rbeta(40,3,2)
  p_value <- as.numeric(ad.test(beta_40sample)[2])
  if (p_value < 0.05) {
    beta_ad_truepos = beta_ad_truepos + 1
  }
}
ad_power <- beta_ad_truepos/1000

### 4.3 ----
data(golub)
head(golub)
?golub

tumour_class <- factor(golub.cl, levels=0:1, labels=c("ALL", "AML"))
golub_ALL <- golub[,tumour_class=="ALL"]

### Null hypothesis: average expression of Cyclin D3 is not greater than the average of all gene expression level changes
### t test
cyclinD3 <- golub_ALL[1042,]
all_expt_cyclinD3 <- golub_ALL[-1042,]
hist(cyclinD3)

# can use t-test off the bat because >30 samples and therefore central data theorum holds true
t.test(cyclinD3, all_expt_cyclinD3, alternative="greater")

# result: reject the null, average expression IS GREATER in Cyclin3

### 4.4 ----
not_normal <- 0
for (i in 1:3051) {
  current_gene_set <- golub_ALL[i,]
  p_value <- as.numeric(ad.test(current_gene_set)[2])
  if (p_value < 0.05) {
    not_normal <- not_normal + 1
  }
}
not_normal_proportion <- not_normal/3051

### 4.5 ----
golub_AML <- golub[,tumour_class=="AML"]
golub_ALL <- golub[,tumour_class=="ALL"]

par(pin=c(3,3),mfrow=c(1,2))
hist(golub_AML[1042, ], 
     breaks=5,
     main="CyclinD3 AML Expression",
     xlab="Gene Expression",
     xlim=range(-1.0,3.0),
     ylab="Frequency")
hist(golub_ALL[1042, ], 
     breaks=5,
     main="CyclinD3 ALL Expression",
     xlab="Gene Expression",
     xlim=range(-1.0,3.0),
     ylab="Frequency")

wilcox.test(golub_AML[1042, ],golub_ALL[1042, ],alternative="two.sided")
t.test(golub_AML[1042, ],golub_ALL[1042, ],alternative="two.sided")
ks.test(golub_AML[1042, ],golub_ALL[1042, ],alternative="two.sided")

# Null hypothesis for the t-test is that the MEANS for the two datasets are the same (equal to 0)
# Null hypothesis for the MWU test is that MEDIANS for the two datasets are the same
# Null hypothesis for the KS test is that both samples are from the same distribution

# first check the ratio of standard deviations of the two samples
sd(golub[1042,tumour_class=="ALL"])/sd(golub[1042,tumour_class=="AML"])
mean(golub[1042,tumour_class=="ALL"])
mean(golub[1042,tumour_class=="AML"])
# t.test for difference in means with equal variances
t.test(golub[1042,tumour_class=="ALL"],golub[1042,tumour_class=="AML"],alternative="two.sided",var.equal=TRUE)

# test for difference in medians
wilcox.test(golub[1042,tumour_class=="ALL"],golub[1042,tumour_class=="AML"])
# test for differences in distribution
ks.test(golub[1042,tumour_class=="ALL"],golub[1042,tumour_class=="AML"])

# in this case all test conclude that their Null hypothesies can be rejected (and importantly you can see in the plot that
# the two distributions are different!)

# SESSION 5 ----
## Bootstrapping ----
library(boot)
## an example loop for bootstrapping means
for (i in 1:1000) {
  bootstrap_sample <- sample(original_sample, replace=TRUE)
  bootstrap_mean <- mean(bootstrap_sample) #can be other parameters as well
}

# boot library can automatically implement CI corrections

# for difference in means between two samples
x <- rnorm(10,10,2)
y <- rnorm(10,15,3)

for (i in 1:10000) {
  boot_x <- sample(x,replace=TRUE)
  boot_y <- sample(y,replace=TRUE)
  difference <- c(difference,mean(boot_x)-mean(boot_y))
}


# SESSION 6 ----
## Chi Squared Test ----
dchisq()
pchisq()
qchisq()
rchisq()

p <- c(0.5,0.5)
obs <- c(12,18)
a=sum(obs)
exp <- a*p
x.sq <- sum((obs-exp)^2/exp)
x.sq
qchisq(0.95,1)

## Practical 5 ----
### 5.1 Bootstrapping ----
sample_size<-50
x<-rnorm(sample_size,100,5)
sample_means<-NULL
for (i in 1:5000) {
  sample_means <- c(sample_means,mean(sample(x,replace=TRUE)))
}
hist(sample_means,
     breaks=seq(90,110,0.4),
     main="bootstrap means",
     col="grey",
     xlab="mean value")
sprintf("Bootstrap 95 confidence interval lies between %06.3fand %06.3f", 
        quantile(sample_means,0.025),
        quantile(sample_means,0.975))

# comparison to t-distribution
mean(x) + qt(0.025, sample_size-1) * (sd(x) / sqrt(sample_size))
mean(x) + qt(0.975, sample_size-1) * (sd(x) / sqrt(sample_size))

### 5.2 ----
chol <- c(21.0, -3.25, -10.75, -13.75, -32.50, -39.50, -41.75, -56.75, -80.0)

mean_chol <- mean(chol)
sample_means<-NULL
for (i in 1:5000) {
  sample_means <- c(sample_means,mean(sample(chol,replace=TRUE)))
}

# standard error of the mean: half the central 68.2% of the distribution
quantile(sample_means, 0.841)
quantile(sample_means, 0.159)
SEM <- (quantile(sample_means, 0.841) - quantile(sample_means, 0.159)) / 2
# approx. equal to sd(sample_means)

quantile(sample_means,0.025)
quantile(sample_means,0.975)
# 0 is not within the 95% CI for the bootstrapped samples, therefore reject the Null - mean change is significantly different from 0

median_chol <- median(chol)
sample_medians<-NULL
for (i in 1:5000) {
  sample_medians <- c(sample_medians,median(sample(chol,replace=TRUE)))
}
quantile(sample_medians, 0.025)
quantile(sample_medians, 0.975)
# 0 is not within the 95% CI for the bootstrapped samples, therefore reject the Null - median change is significantly different from 0

### 5.3 ----
#separately resample the groups and produce means from this to compare
chol_f <- c(23.0, 10.25, -3.75, -13.0, -23.50, -29.50, -31.0, -43.75)

diff <- NULL
for (i in 1:5000) {
  sample_m <- sample(chol,replace=TRUE)
  sample_f <- sample(chol_f,replace=TRUE)
  diff <- c(diff, mean(sample_m)-mean(sample_f))
}
quantile(diff, 0.025)
quantile(diff, 0.975)
# 0 is within the 95% CIs, therefore accept the null. Male and Female patient average results are the same.

### 5.4 ----
# Null hypothesis: inherited characteristics are blended, not using a dominant/recessive model where Round is dominant
# 1 degree of freedom here, critical value = 3.84
qchisq(0.95,1)

p <- c(0.75,0.25)
obs <- c(5474, 1850)
a <- sum(obs)
exp <- a*p
x.sq <- sum((obs-exp)^2/exp)
x.sq
# accept the null, X^2 is smaller than the critical value of 3.84

chisq.test(obs, p=p)

### 5.5 ----
# 3 degrees of freedom
p <- c(0.25, 0.25, 0.25, 0.25)
obs <- c(410, 789, 573, 394)
total <- sum(obs)

x.sq <- 0
for (i in 1:4) {
  x.sq <- x.sq + ((obs[i] - (p[i] * total))^2 / (p[i] * total))
}
x.sq

qchisq(0.95,3)
# X^2 is larger than critical value for 0.95 and 3 degrees of freedom therefore reject the null.

### 5.6 ----
one_sample_test <- function(x, y) {
  exp <- sum(x)*y
  x.sq <- sum((x-exp)^2/exp)
  return (x.sq)
}
one_sample_test(c(5474, 1850),c(0.5,0.5))

### 5.7 ----
data <- matrix(c(12,10,19,4),2,byrow=TRUE)

uncorrected_fisher <- fisher.test(data)[[1]]
prob_observed <- (choose(23,19) * choose(22,12) / choose(45,31))
p_correction <- prob_observed/2
corrected_fisher <- uncorrected_fisher - p_correction

fisher_corrector <- function(data) {
  uncorrected_fisher <- fisher.test(data)[[1]]
  marg_row1 <- data[1,1] + data[1,2]
  marg_row2 <- data[2,1] + data[2,2]
  marg_col1 <- data[1,1] + data[2,1]
  total <- marg_row1 + marg_row2
  
  prob_observed <- (choose(marg_row2, data[2,1]) * choose(marg_row1, data[1,1]) / choose(total, marg_col1))
  corrected_fisher <- uncorrected_fisher - (prob_observed/2)
  
  return(corrected_fisher)
}

fisher_corrector(data)

# SESSION 7 ----
## One-way ANOVA and F distributions ----
df(f,2,12)
pf(f,1,8)
qf(p,1,13)
rf(n,3,29)
oneway.test(values~index, data, var.equal=FALSE) #var.equal = Welch's correction, for testing in two groups. 

## Pearson's correlation ----
## product-moment correlation for points x and y
cor.test(x,y)

## Practical 6 ----
### 6.1 ----
sample1 <- rnorm(5,0,1)
sample2 <- rnorm(5,0,1)
wilcox.test(sample1, sample2, alternative="two.sided", conf.level=0.95)

p_acceptance <- NULL
for (i in 1:1000) {
  sample1 <- rnorm(5,0,1)
  sample2 <- rnorm(5,0,1)
  wilcox <- wilcox.test(sample1, sample2, alternative="two.sided", conf.level=0.95)[3]
  if (wilcox < 0.05) {
    p_acceptance <- c(p_acceptance, wilcox)
  }
}
length(p_acceptance) #19 for this sample

p_acceptance <- NULL
for (i in 1:1000) {
  sample1 <- rnorm(1000,0,1)
  sample2 <- rnorm(1000,0,1)
  wilcox <- wilcox.test(sample1, sample2, alternative="two.sided", conf.level=0.95)[3]
  if (wilcox < 0.05) {
    p_acceptance <- c(p_acceptance, wilcox)
  }
}
length(p_acceptance) # 48 for this sample - 0.048 of the time, you reject the null. Which nearly matches the sig level of 0.05.


### 6.2 ----
### Calculate the mean expression level for each gene, the overall mean expression level and then the test statistic
geneA<-c(2, 3, 1, 2)
geneB<-c(8, 7, 9, 8)
geneC<-c(11, 12, 13, 12)
genes <- data.frame(A=geneA,B=geneB,C=geneC)

totmean <- mean(c(geneA,geneB,geneC))
circular <- c("A","B","C")
top_f <- 0
bottom_f <- 0

for (i in 1:3) {
  meani <- mean(genes[,i])
  top_f <- top_f + ((4*(meani - totmean)^2) / 2) #have now corrected, needed to * individual mean-totmean by number of values in the dataset (4 in this case)
  for (j in genes[,i]) {
    bottom_f <- bottom_f + (((j - meani)^2) / 9)
  }
}
Fstat <- top_f / bottom_f #calculated to be 38? Not correct, but gets the right result. Check answers when released.

#### Correct answer 6.2 ----
yA <- mean(geneA)
yB<- mean(geneB)
yC <- mean(geneC)
y <- mean(c(geneA,geneB,geneC))
F <- ((4*(yA-y)^2 + 4*(yB-y)^2 + 4*(yC-y)^2)/2)/((sum((geneA-yA)^2)+ sum((geneB-yB)^2)+sum((geneC-yC)^2))/9)

top_correct <- (4*(mean(geneA)-totmean)^2 + 4*(mean(geneB)-totmean)^2 + 4*(mean(geneC)-totmean)^2) / 2
bottom_correct <- (sum((geneA-mean(geneA))^2)+ sum((geneB-mean(geneB))^2)+sum((geneC-mean(geneC))^2)) / 9
Fcorrect <- top_correct / bottom_correct

### Plot the cumulative distribution of the appropriate F-distribution for testing the null hypothesis that the means of all three samples are equal.
x <- seq(0,10,0.1)
plot(x,pf(x,2,9),type="l")
segments(0,0.95,x1=qf(0.95,2,9),col="red")
segments(qf(0.95,2,9),0.0,qf(0.95,2,9),0.95,col="black",lty=2)

### What is the critical value for rejecting the null if the chosen significance level is 5%?
qf(0.95,2,9) # 4.256495

### Are these means of these three samples equal?
### if you do it right, F=152, which is much higher than the critical value and therefore you reject the null
### means are not the same
oneway.test(values~ind,stack(genes),var.equal=TRUE)

### multiple testing
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("mutoss"))
library(mutoss)

result <- regwq(values~ind,stack(genes),alpha=0.05)
pairwise_rejection <- cbind(result$confIntervals[,0],result$rejected)
# all pairs have a pairwise_rejection of 1, therefore reject the null for all pairs of groups (no equal means for any pair of groups)


### 6.3 ----
BiocManager::install(c("ALL")) 
library(ALL) 
data(ALL) 
ALL_stagesB1B2B3<- ALL[,ALL$BT %in% c("B1","B2","B3")] 
expression_SKI <- exprs(ALL_stagesB1B2B3)["1866_g_at",] 
expression_SKI_stages<-data.frame(cbind(expression_SKI, ALL_stagesB1B2B3$BT))
oneway.test(expression_SKI~V2,expression_SKI_stages,var.equal=TRUE)[[3]] 

### How many patients are in each of the B, B1, B2, B3 and B4 groups
table(ALL$BT)

### Construct a subset of the ALL data frame containing the B-cell ALL patients in stage B, B1, B2, B3, B4 from the ALL data. 
ALL_Bstages<- ALL[,ALL$BT %in% c("B","B1","B2","B3","B4")] 
table(ALL_Bstages$BT)

### Using a for loop over the row numbers (or the names), collect (in a vector) the p-values from a one-way ANOVA across these five groups for each and every gene. 
qf(0.95,4,90) # F stat critical value: 2.472927

p_values <- NULL
for (i in row.names(ALL)) {
  expression_values <- exprs(ALL_Bstages)[i,]
  expression_gene_i <- data.frame(cbind(expression_values, ALL_Bstages$BT))
  
  bottom_f <- 0
  top_f <- 0
  for (i in 1:5) {
    values_table <- pivot_wider(expression_gene_i, names_from=V2, values_from=expression_values)
    values_totmean <- mean(c(unlist(values_table[[1]]),
                             unlist(values_table[[2]]),
                             unlist(values_table[[3]]),
                             unlist(values_table[[4]]),
                             unlist(values_table[[5]])))
    meani <- mean(unlist(values_table[[i]]))
    top_f <- top_f + ((4*(meani - values_totmean)^2) / 4)
    for (j in unlist(values_table[[i]])) {
      bottom_f <- bottom_f + (((j - meani)^2) / 90)
    }
  }
  Fstat <- top_f / bottom_f
  if (Fstat > 2.472927) {
    p_values <- c(p_values, oneway.test(expression_values~V2, expression_gene_i, var.equal=FALSE)[[3]])
  } else {
    p_values <- c(p_values, oneway.test(expression_values~V2, expression_gene_i, var.equal=TRUE)[[3]])
  }
}

### a version that assumes equal variances for each expressed gene
p_values <- NULL
for (i in row.names(ALL)) {
  expression_values <- exprs(ALL_Bstages)[i,]
  expression_gene_i <- data.frame(cbind(expression_values, ALL_Bstages$BT))
  p_values <- c(p_values, oneway.test(expression_values~V2, expression_gene_i, var.equal=TRUE)[[3]])
}

### Which genes have smaller p-values than 0.000004 (the Bonferroni value for a family wise error of 0.05)?
count <- 0
for (i in p_values) {
  if (i < 0.000004) {
    count <- count + 1
  }
}
count

### Collect the p-values from the Kruskal-Wallis test (kruskal.test) in a vector. 
p_kw_values <- NULL
for (i in row.names(ALL)) {
  expression_values <- exprs(ALL_Bstages)[i,]
  expression_gene_i <- data.frame(cbind(expression_values, ALL_Bstages$BT))
  p_kw_values <- c(p_kw_values, kruskal.test(expression_gene_i$expression_values, expression_gene_i$V2)[[3]])
}

### How many genes have p-values smaller than 0.001 from both ANOVA and Krusal-Wallis? Which genes are they?
sig_genes <- NULL
sig_gene_index <- NULL
for (i in 1:12625) {
  if (p_values[[i]] < 0.001 & p_kw_values[[i]] < 0.001) {
    sig_gene_index <- c(sig_gene_index, i)
    sig_genes <- c(sig_genes, rownames(ALL_Bstages[i]))
  }
}


# SESSION 8 ----
## Regression ----
## sample data from true distribution y=3-2x plus a simulated normally distributed measurement error
x <- seq(1,10,0.5)
y <- 3.0 - 2*x + rnorm(length(x),0,1.2)
plot(x,y,xlim=c(0,10), ylim=c(-20,5))

## fit data to a linear model
model <- lm(y~x)
model_quadratic <- lm(y~x+xsq) #how to do a quadratic model
summary(model)

# Residual standard error = estimated sd of residual distribution
# F-statistic = use to test null hypothesis

model<-lm(y~x) # create a model object where y depends linearly on x
abline(model,col="red") # least squares regression line
coefficients(model) # recover estimate of true parameters
fitted(model) # estimated values of y for the model
residuals(model) # return the residuals (y observed – y estimated)
summary(model) # returns lengthy summary of model (inc. t tests)
anova(model) # analysis of variance and F-statistic for parameter(s)
anova(base,full) # partial-F test of two models ‘base’ and ‘full’
cor(x,y) # compute correlation coefficient (default: Pearson)

## plot of the residuals
plot (x,residuals(model), xlim=c(0,10))


## Practical 7 ----
### 7.1 ----
# Create a scatterplot of these data points
x <- runif(5,0,10)
y <- 1.2 + 0.5*x + rnorm(5,0,1.3)
plot(x, y)

# Compute the Pearson correlation coefficient. Is the correlation coefficient significantly different from zero ? Hint : use the appropriate t statistic and t test from Lecture 7
r <- cor(x,y) # correlation not sig. different from 0
t <- r * sqrt((5-2)/(1-r^2)) ;t # this is the t statistic for the correlation
qt(0.025, 5-2) # two-tailed test against critical values of t for N-2 degrees of freedom
qt(0.975, 5-2)

cor.test(x,y) #alt version

# Fit these data points to a straight line model and add the best fit line to the scatterplot
model<- lm(y~x)
abline(model,col="red")

# Find the fitted parameters for the straight line. 
coefficients(model)
# Is the slope significantly different from zero? 
anova(model)
# Are the parameters of the model close to those of the equation from which the data were derived ?

residuals <- residuals(model)
plot(x,residuals)

### 7.2 simple linear regression ----
# For each datasets (A,B,C,D) find the coefficients of the best fit straight line. 
# For each dataset is the slope significantly different from zero? 
# Assess the ‘goodness’ of fit of each the straight line to each dataset. 

anscombeA_x<-c(10,8,13,9,11,14,6,4,12,7,5)
anscombeA_y<-c(8.04,6.95,7.58,8.81,8.33,9.96,7.24,4.26,10.84,4.82,5.68)
anscombeB_x<-c(10,8,13,9,11,14,6,4,12,7,5)
anscombeB_y<-c(9.14,8.14,8.74,8.77,9.26,8.1,6.13,3.10,9.13,7.26,4.74)
anscombeC_x<-c(10,8,13,9,11,14,6,4,12,7,5)
anscombeC_y<-c(7.46,6.77,12.74,7.11,7.81,8.84,6.08,5.39,8.15,6.42,5.73)
anscombeD_x<-c(8,8,8,8,8,8,8,19,8,8,8)
anscombeD_y<-c(6.58,5.76,7.71,8.84,8.47,7.04,5.25,12.5,5.56,7.91,6.89)

model_A <- lm(anscombeA_y~anscombeA_x)
model_B <- lm(anscombeB_y~anscombeB_x)
model_C <- lm(anscombeC_y~anscombeC_x)
model_D <- lm(anscombeD_y~anscombeD_x)

# Dataset A
coefficients(model_A)
anova(model_A)
summary(model_A) # slope sig = 0.00217, significantly different from 0
# F-stat p-value: 17.99, p-value low, reject the null. slope sig. different from 0. 

# plot of the residuals
plot(anscombeA_x, residuals(model_A)) # good plot, appears normal about 0
plot(anscombeB_x, residuals(model_B)) # functional form looks quadratic
plot(anscombeC_x, residuals(model_C)) # there is an outlier, all other data would fit a straight line
plot(anscombeD_x, residuals(model_D)) # The variance is dependant on x


### 7.3 polynomial regression ----
### Attempt to improve the model for anscombeB by adding a term in x2 to the model 
anscombeB_xsq<-anscombeB_x*anscombeB_x 
quad_modelB<-lm(anscombeB_y~anscombeB_x + anscombeB_xsq) 
coefficients(quad_modelB) 

### Add the best fit quadratic curve to a scatterplot of the anscombeB data + plot the residuals
plot(anscombeB_x,anscombeB_y) 
sortedx<-sort(anscombeB_x) 
lines(sortedx, -5.9957343 + (2.7808392*sortedx)-(0.1267133 * sortedx*sortedx),col="red") 
plot(anscombeB_x,residuals(quad_modelB))

### Is the goodness of fit R2 improved by adding the term in x2? Use a partial-f test to test if the added term significantly improves the fit. 
summary(model_B) 
summary(quad_modelB) 
anova(model_B,quad_modelB) # looking at Pr(>F), fit is partially improved, can reject the null that the additional parameter is 0. 


### 7.4 multiple linear regression ----
data.home <- "http://people.cryst.bbk.ac.uk/~ubcg66a/statistics/" 
expdata.url <- paste(data.home,"ecoli_expression.dat",sep="") 
expdata <- read.table(url(expdata.url),header=T) 
expdata$llen <- log(expdata$len) 
names(expdata)

# method below shows how to step up from the simplest model, or you could start with all parameters and remove gradually.
# Show that len and llen are (unsurprisingly) the most correlated explanatory variables by testing for correlation between all pairs of the explanatory variables. 
pairs(expdata)
cor(expdata, method="spearman") # only using spearman because of the shape of the pairs of scatter plots (non-normal)

# First regress exp against len and exp against llen. 
# Which of len and llen gives a better fit? What is the evidence for this? 
# Which of len and llen should be omitted from the model? 
reg.len <- lm(exp ~ len,data=expdata) 
summary(reg.len) # r squared = 0.14

reg.llen <- lm(exp ~ llen,data=expdata) 
summary(reg.llen) # r squared = 0.22, fit is better with llen, remove len from model

# Now add gc, cai, gravy, arom to a linear model. 
# Use an ANOVA table(s) to decide which variables contribute to the explanation and which not. Can one leave out any of the variables without losing too much in F value? 
reg.lm <- lm(exp ~ llen + gc + cai + gravy + arom,data=expdata) 
anova(reg.lm) 
# all variables look significant but some are correlated from the correlation test.

# Could use an F test for specific subsets to provide evidence for/against their selection.
exp_llen_gc_cai_gravy_arom_llencai = 
  lm(exp~llen+gc+cai+gravy+arom+llen:cai,data=expdata) 
exp_llen_gc_cai_arom_llencai = lm(exp~llen+gc+cai+arom+llen:cai,data=expdata) 
anova( exp_llen_gc_cai_arom_llencai, 
       exp_llen_gc_cai_gravy_arom_llencai)$"Pr(>F)"[2] 

# Now add pairwise interactions between the remaining terms, which ones are significant? 
reg.lm <- lm(exp ~ llen + gc + cai + llen:gc + llen:cai + gc:cai,data=expdata) 
anova(reg.lm) # llen and cai are most significant interactions

# Final proposed model 1:
reg.lm_1 <- lm(exp ~ llen + gc + cai + llen:cai, data=expdata) 
reg.lm_1
# Final proposed model 2:
reg.lm_2 <- lm(exp ~ llen + gc + arom + cai + llen:cai, data=expdata) 
reg.lm_2

anova(reg.lm_1, reg.lm_2)
# This shows that the arom parameter is significantly different from 0, even with a small effect size, and so is beneficial to include

# Plot your predicted values for expression levels against the experimental values
plot(9.319-0.999*expdata$llen-9.792*expdata$gc-2.886*expdata$arom-1.577*expdata$cai+0.895*expdata$llen*expdata$cai,expdata$exp) 


# SESSION 9 ----
prcomp(datamatrix) # create a PCA matrix or frame of data
kmeans(datamatrix, centers=k) # perform k-means clustering
cmdscale(datamatrix, method =”name_of_method” k = value) # metric MDS
isoMDS(datamatrix,k = value) # isoMAP (non-metric MDS)

## PCA ----
## 2-dim correlated data points from bi-normal distribution
expX<-rnorm(100,6,3) 
expY<- 2+1.1*expX+rnorm(100,0,1)

## plot data
par(pin=c(3,3),font=2,ps=10,family="sans") 
plot(expX,expY,pch=21, xlab="Expression level gene X",ylab="Expression level gene Y") 
expmatrix<-cbind(expX,expY) 

## PCA analysis
exppca <- prcomp(expmatrix, center = TRUE) 
plot(exppca,main = "PCA eigenvalues")
par(pin=c(5,2),font=2,ps=10,family="sans") 
plot(exppca$x[, 1], exppca$x[, 2], pch=21, main = "PCA transformed data", xlab = "PC1", ylab = "PC2") 

## PCA with scaling of data (equalising variances)
exppca <- prcomp(expmatrix, center = TRUE, scale = TRUE) 
plot(exppca,main = "PCA eigenvalues")
par(pin=c(5,2),font=2,ps=10,family="sans") 
plot(exppca$x[, 1], exppca$x[, 2], pch=21, main = "PCA transformed data", xlab = "PC1", ylab = "PC2") 

## example with 3 groups
geneX_groupA <- rnorm(200,3.2,0.8) 
geneY_groupA <- rnorm(200,1.6,0.8) 
geneZ_groupA <- rnorm(200,1.6,0.8) 
geneX_groupB <- rnorm(200,5.0,0.8) 
geneY_groupB <- rnorm(200,3.0,0.8) 
geneZ_groupB <- rnorm(200,3,0.8) 
geneX_groupC <- rnorm(60,4.0,0.6) 
geneY_groupC <- rnorm(60,3.0,0.6) 
geneZ_groupC <- rnorm(60,0.0,0.6) 

collected_Xvalues <- c(geneX_groupA,geneX_groupB,geneX_groupC) 
collected_Yvalues <- c(geneY_groupA,geneY_groupB,geneY_groupC) 
collected_Zvalues <- c(geneZ_groupA,geneZ_groupB,geneZ_groupC) 
grouped_data <- data.frame(x=collected_Xvalues,y=collected_Yvalues,z=collected_Zvalues, group=factor(c(rep("A",200),rep("B",200),rep("C",60)))) 

plot(grouped_data[,1],grouped_data[,2], pch=21, xlab="Expression level gene X",ylab="Expression level gene Y",bg=c("blue","cyan","green")[grouped_data[,4]], main="") 
plot(grouped_data[,1],grouped_data[,3], pch=21, xlab="Expression level gene X",ylab="Expression level gene Z",bg=c("blue","cyan","green")[grouped_data[,4]], main="") 
plot(grouped_data[,2],grouped_data[,3], pch=21, xlab="Expression level gene Y",ylab="Expression level gene Z",bg=c("blue","cyan","green")[grouped_data[,4]], main="") 

### PCA analysis on the first 2 eigenvectors
pca <- prcomp(grouped_data[1:3], center = TRUE, scale = TRUE) 
plot(pca$x[, 1], pca$x[, 2], pch=21, bg=c("blue","cyan","green")[grouped_data[,4]], main = "PCA", xlab = "PC1", ylab = "PC2") 


## k-means ----
## dataframe with PCA output coordinates (from above example)
mydata <- cbind(collected_Xvalues, collected_Yvalues, collected_Zvalues)

## k-means repeated for each cluster to find the optimum number of clusters - you want one which minimises within group sum of squares and cluster number
withingroups_stress <- (nrow(mydata)-1)*sum(apply(mydata,2,var)) 
for (i in 2:8) withingroups_stress[i] <- 
  sum(kmeans(mydata,centers=i)$withinss) 
par(pin=c(2,2),font=2,ps=10,family="sans") 
plot(1:8, withingroups_stress[1:8], type="b", xlab="Number of Clusters", ylab="Within groups sum of squares") 

## perform k-means with chosen cluster number
myclusters <- kmeans(mydata,3)

## plot clusters
library(cluster) 
clusplot(mydata, myclusters$cluster, color=TRUE, shade=TRUE, labels=1, lines=0) 

## investigate the agreement or not between the original group id and the results of clustering
mynewdata <- data.frame(mydata, myclusters$cluster,group=factor(c(rep("A",200),rep("B",200),rep("C",60)))) 
mynewdata


## MDS (Multidimensional Clustering) ----
## Calculate euclidean distances between the rows
distances <- dist(mydata, method = "euclidian")

## scaling function and plotting
mdsdata <- cmdscale(distances,eig=TRUE, k=2) 
plot(mdsdata$points[,1], mdsdata$points[,2], type = "p", pch=21, 
     bg=c("blue","cyan","green")[grouped_data[,4]], 
     main = "MDS of gene expression levels", xlab = "Coordinate 1", 
     ylab = "Coordinate 2") 

## Non-metric MDS (Isomap)
library(MASS) 
isomdsdata <- isoMDS(distances, k=2) 
plot(isomdsdata$points[,1], isomdsdata$points[,2], type = "p", pch=21, 
     bg=c("blue","cyan","green")[grouped_data[,4]], 
     main = "isoMDS of gene expression levels", xlab = "Coordinate 1", 
     ylab = "Coordinate 2") 

