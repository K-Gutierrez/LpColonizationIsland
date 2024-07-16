install.packages("lawstat")
library(lawstat)

group1 <- c(1.68, 4.84, 3.28, 2.22, 0.65, 1.71)
group2 <- c(4.51, 3.16, 3.81, 1.36, 1.93, 0.80, 2.32, 5.47, 0.10, 6.94, 1.60, 5.29, 3.12, 4.84, 3.26, 5.05, 5.98, 3.85, 2.13, 4.10)


result <- brunner.munzel.test(group1, group2)
print(result)
