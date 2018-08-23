# simple lm vs logistic regression
y = rbinom(1000,1,0.8)
x = matrix(rnorm(length(y)*100),nrow=length(y))
for(j in 1:10){
  x[,j] = x[,j] + y*2
}

ps1 = c();ps2 = c()
for(j in 1:ncol(x)){
  currx = x[,j]
  l1 = summary(lm(y~currx))
  ps1[j] = l1$coefficients[2,4]
  l2 = summary(glm(y~currx,family=binomial(link='logit')))
  ps2[j] = l2$coefficients[2,4]
}
plot(ps1,ps2)
hist(ps1)
