# Bad-SVM
bad_svm <- function() {

# Generates bad_svm dataset in 2-D space
innerCircle <- generateCircle(500,1,2)
outerCircle <- generateCircle(700,4,5)
concentricCircles <- rbind(innerCircle,outerCircle)
concentricCircles1 <- cbind(concentricCircles, cluster)
trainData <- rbind(concentricCircles1[1:300,],concentricCircles1[501:800,])
testData <- rbind(concentricCircles1[301:500,],concentricCircles1[801:1200,])
cluster <- as.matrix(c(rep(1,500),rep(2,700)))
plot(concentricCircles, col = cluster, main="BAD SVM Dataset", xlab = "X-Axis" ,ylab = "Y-Axis")

# Performs un-kernelized Support Vector Machines on bad_svm data
sv <- ksvm(concentricCircles,cluster, type = "C-svc",kernel="vanilladot")
#sv <- ksvm(trainData[,1:2], trainData[,3], type = "C-svc",kernel="vanilladot")
svmCluster <- round(predict(sv,concentricCircles))
plot(sv, data = concentricCircles, main="BAD SVM linear classification Plot", xlab = "X-Axis" ,ylab = "Y-Axis")

# Calculating performance metrics of the un-kernelized svm
confusionMatrixSVM <- getConfusionMatrix(cluster,svmCluster,unique(cluster))
confusionMatrixSVM
contingencyTableMetrics(confusionMatrixSVM)
twoCrossConfusionMatrixMetrics(confusionMatrixSVM)

# Using kernel SVM methods on bad_svm data and calculating performance metrics
ksv <- ksvm(concentricCircles, cluster, type = "C-svc",kernel="laplacedot")
plot(ksv, data = concentricCircles)
ksvmCluster <- predict(ksv,concentricCircles)
plot(concentricCircles, col = ksvmCluster )

confusionMatrixKSVM <- getConfusionMatrix(cluster,ksvmCluster,unique(cluster))
confusionMatrixKSVM
contingencyTableMetrics(confusionMatrixKSVM)
twoCrossConfusionMatrixMetrics(confusionMatrixKSVM)


ksv <- ksvm(concentricCircles, cluster, type = "C-svc",kernel="rbfdot")
plot(ksv, data = concentricCircles)
ksvmCluster <- predict(ksv,concentricCircles)
plot(concentricCircles, col = ksvmCluster )

confusionMatrixKSVM1 <- getConfusionMatrix(cluster,ksvmCluster,unique(cluster))
confusionMatrixKSVM1
contingencyTableMetrics(confusionMatrixKSVM1)
twoCrossConfusionMatrixMetrics(confusionMatrixKSVM1)
}