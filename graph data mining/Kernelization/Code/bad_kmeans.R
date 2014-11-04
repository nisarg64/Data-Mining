# Bad Kmeans
bad_kmeans <- function() {

# Generates bad_kmeans dataset in 2-D
Sphere1 <- generate2DSphere(500, 5,c(4,4))
Sphere2 <- generate2DSphere(700,2,c(9,9))
Spheres <- rbind(Sphere1,Sphere2)
plot(Spheres, col = cluster, xlab = "X-Axis" ,ylab = "Y-Axis", main = "BAD-kmeans Dataset")

# Performs k-means clustering in bad_kmeans dataset
kc1 <- kmeans(Spheres, 2)
plot(Spheres,col=kc1$cluster, xlab = "X-Axis" ,ylab = "Y-Axis", main = "Bad Kmeans Clustering")

# Calculating performance metrics of the un-kernelized kmeans
confusionMatrix1 <- getConfusionMatrix(cluster,kc1$cluster,unique(cluster))
confusionMatrix1
contingencyTableMetrics(confusionMatrix1)
twoCrossConfusionMatrixMetrics(confusionMatrix1)

# Using kernel kmeans on bad_kmeans data and calculating performance metrics
kc2 <- kkmeans(Spheres,2,kernel = "rbfdot")
plot(Spheres,col=kc2,main="Kernel Kmeans Clustering - rbfdot kernel", xlab = "X-Axis" ,ylab = "Y-Axis")
confusionMatrix2 <- getConfusionMatrix(cluster,kc2,unique(cluster))
confusionMatrix2
contingencyTableMetrics(confusionMatrix2)
twoCrossConfusionMatrixMetrics(confusionMatrix2)

kc3 <- kkmeans(Spheres,2,kernel = "anovadot")
plot(Spheres,col=kc3,main="Kernel Kmeans Clustering - anovadot kernel", xlab = "X-Axis" ,ylab = "Y-Axis")
confusionMatrix3 <- getConfusionMatrix(cluster,kc3,unique(cluster))
confusionMatrix3
contingencyTableMetrics(confusionMatrix3)
twoCrossConfusionMatrixMetrics(confusionMatrix3)

kc4 <- kkmeans(Spheres,2,kernel = "besseldot")
plot(Spheres,col=kc4,main="Kernel Kmeans Clustering - besseldot kernel", xlab = "X-Axis" ,ylab = "Y-Axis")
confusionMatrix4 <- getConfusionMatrix(cluster,kc4,unique(cluster))
confusionMatrix4
contingencyTableMetrics(confusionMatrix4)
twoCrossConfusionMatrixMetrics(confusionMatrix4)

}