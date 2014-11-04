#BAD PCA
bad_pca <- function(){

# Generates bad_pca dataset in 2-D space
outerCircle1 <- generateCircle(500,6,7)
Sphere21 <- generate2DSphere(700,2,c(0,0))
pcaCircle <- rbind(outerCircle1,Sphere21)
plot(pcaCircle,col=cluster, xlab = "X-Axis" ,ylab = "Y-Axis", main = "BAD-pca Dataset")

# Performs un-kernelized principal component analysis on bad_pca data
pca <- princomp(pcaCircle,centers = TRUE)
summary(pca)
plot(pca, main = "BAD-pca variance")
plot(pca$scores[,1],col = cluster,ylab = "Principal Comp. 1", main = "Dataset after applying PCA 1-D")

# Using kernel PCA on bad_pca data set and calculating performance metrics 
kpc <- kpca(~.,as.data.frame(pcaCircle),kernel ="besseldot",kpar=list(sigma = 0.18))
eig(kpc)
plot(pcv(kpc)[,3],pcv(kpc)[,2], col = cluster,main="Kernel PCA - besseldot kernel", xlab = "Principal Component 1" ,ylab = "Principal Component 2")
eigenRatio1 <- eig(kpc)/sum(eig(kpc))

kpc <- kpca(~.,as.data.frame(pcaCircle),kernel ="laplacedot",kpar=list(sigma = 0.0008))
plot(pcv(kpc)[,3],pcv(kpc)[,2], col = cluster,main="Kernel PCA - laplacedot kernel", xlab = "Principal Component 1" ,ylab = "Principal Component 2")
eigenRatio2 <- eig(kpc)/sum(eig(kpc))

}