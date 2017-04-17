# Perform statistical analysis and t tests on simulated datasets.
w_data <-read.csv("fb_weighting_analysis.csv")

tapply(w_data$Infected, w_data$Weighted, FUN=mean)
tapply(w_data$Infected, w_data$Weighted, FUN=sd)

tapply(w_data$Dead, w_data$Weighted, FUN=mean)
tapply(w_data$Dead, w_data$Weighted, FUN=sd)

t.test(w_data$Infected[w_data$Weighted=="False"],
       w_data$Infected[w_data$Weighted=="True"])

t.test(w_data$Dead[w_data$Weighted=="False"],
       w_data$Dead[w_data$Weighted=="True"])

u_data <-read.csv("random_weighting_analysis.csv")

tapply(u_data$Infected, u_data$Weighted, FUN=mean)
tapply(u_data$Infected, u_data$Weighted, FUN=sd)

tapply(u_data$Dead, u_data$Weighted, FUN=mean)
tapply(u_data$Dead, u_data$Weighted, FUN=sd)

t.test(u_data$Infected[u_data$Weighted=="False"],
       u_data$Infected[u_data$Weighted=="True"])

t.test(u_data$Dead[u_data$Weighted=="False"],
       u_data$Dead[u_data$Weighted=="True"])
