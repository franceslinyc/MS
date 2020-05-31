library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)
?registerDoParallel
foreach(i = 1:3) %dopar% sqrt(i)

x <- iris[which(iris[,5] != "setosa"), c(1,5)]
trials <- 10000

ptime <- system.time({
    r <- foreach(icount(trials), .combine=cbind) %dopar% {
        ind <- sample(100, 10, replace = TRUE)
        result1 <- glm(x[ind,2] ~ x[ind,1], family = binomial(link = "logit"))
        coefficients(result1)
    }
})

stime <- system.time({
    r <- foreach(icount(trials), .combine=cbind) %do% {
        ind <- sample(100, 10, replace = TRUE)
        result1 <- glm(x[ind,2] ~ x[ind,1], family = binomial(link = "logit"))
        coefficients(result1)
    }
})

ptime; stime
