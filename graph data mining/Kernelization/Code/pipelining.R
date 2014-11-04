# Bad Kmeans 7 dimensional Pipelining
pipelining <- function(){

# This function generates a circle with given inner and outer radius, 
# center and no of points in a 7-D space
generateCirc <- function(n,innerRadius,outerRadius,centre = c(0,0,0,0,0,0,0)){
  x0 <- centre[1] ; y0 <- centre[2] ; z0 <- centre[3] ; a0 <- centre[4] ; b0 <- centre[5] ; c0 <- centre[6] ; d0 <- centre[7]
  theta <- 2*pi*runif(n) 
  r <- sqrt(runif(n,min =innerRadius,max = outerRadius )) 
  cbind(x=(outerRadius - innerRadius)*r*cos(theta)+x0, y=(outerRadius - innerRadius)*r*sin(theta)+y0,z=(outerRadius - innerRadius)*0*sin(theta)+z0,a=(outerRadius - innerRadius)*0*sin(theta)+a0,b=(outerRadius - innerRadius)*0*sin(theta)+b0,c=(outerRadius - innerRadius)*0*sin(theta)+c0,d=(outerRadius - innerRadius)*0*sin(theta)+d0) 
  
}

# Generates bad_kmeans dataset in 7-D
NDSphere1 <- generate7DSphere(500,2)
NDSphere2 <- generateCirc(700,5,8)
NDSpheres <- rbind(NDSphere1,NDSphere2)
plot(NDSpheres,col=cluster, xlab = "X-Axis", ylab = "Y-Axis", main = "Bad_Kmeans Dataset in 7D")

# Performs k-means clustering in bad_kmeans dataset
Nkc <- kmeans(NDSpheres,2)
plot(NDSpheres, col = Nkc$cluster, xlab = "X-Axis", ylab = "Y-Axis", main = "Unkernelized Kmeans Clustering in 7D on Bad_Kmeans Data")

# Calculating performance metrics of the un-kernelized kmeans
ConfusionMatrix <- getConfusionMatrix(cluster,Nkc$cluster, unique(cluster))
ConfusionMatrix
contingencyTableMetrics(ConfusionMatrix)
twoCrossConfusionMatrixMetrics(ConfusionMatrix)

# Using kernel PCA on bad_kmenas data set as a pre processing
Nkpca <- kpca(NDSpheres, kernel = "rbfdot", kpar=list(sigma = 0.03))
eig(Nkpca)
eigenRatio1 <- eig(Nkpca)/sum(eig(Nkpca))
plot(pcv(Nkpca)[,1],pcv(Nkpca)[,2],col = cluster, xlab = "Principal Component 1" , ylab = "Principal Component 2",  main = "Kernel PCA - rbfdot kernel on Bad_Kmeans Data")

# Performs k-means clustering after applying kernel PCA, on M-dimensions
Knkc <- kmeans(pcv(Nkpca)[,1:2],2)
plot(pcv(Nkpca),col = Knkc$cluster,xlab = "X-Axis", ylab = "Y-Axis", main = "Kmeans Clustering after applying Kernelized PCA")

# Calculating performance metrics of the kernelized pca and kmeans
ConfusionMatrix1 <- getConfusionMatrix(cluster,Knkc$cluster,unique(cluster))
ConfusionMatrix1
contingencyTableMetrics(ConfusionMatrix1)
twoCrossConfusionMatrixMetrics(ConfusionMatrix1)

# Using kernel kmeans on bad_kmeans data and calculating performance metrics
kkm <- kkmeans(NDSpheres,2)
plot(NDSpheres,col=kkm,xlab = "X-Axis", ylab = "Y-Axis", main = "Kernelized Kmeans Clustering on Bad_kmeans data")

ConfusionMatrix2 <- getConfusionMatrix(cluster,kkm,unique(cluster))
ConfusionMatrix2
contingencyTableMetrics(ConfusionMatrix2)
twoCrossConfusionMatrixMetrics(ConfusionMatrix2)

}