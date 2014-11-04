contingencyTableMetrics <-
function(mat)
{
if((nrow(mat)!=2)&(ncol(mat)!=2)){
print("Not a 2X2 Matrix")
}
else
{
f11=mat[1,1]
f10=mat[1,2]
f01=mat[2,1]
f00=mat[2,2]
T=sum(mat)
Rand = (f11+f00)/T
JC = (f11)/(f11+f10+f01)
return(list("Rand"=Rand,"Jaccard Coefficient"=JC))
}

}

