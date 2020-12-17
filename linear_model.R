

set.seed(50)

data<-read.csv("trainingData.csv")


trainingRows <- sample(1:nrow(data), 0.975 * nrow(data))
trainingData <- data[trainingRows, ]
testData<-data[-trainingRows,]
fit <- lm(keystoneness ~ rel.abun*cstr*cc.rel*md*eigen.centr*cb, data=trainingData)
summary(fit)
new.fit<-predict(fit,testData[,-1])
plot(testData$keystoneness,new.fit)




