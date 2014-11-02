# Checks if the required packages are installed. I f not, install the packages listed

packages <- c("igraph","cluster","fpc","Matrix","MASS")

for(pkg in packages){
  if(!is.element(pkg, installed.packages()[,1]))
  {install.packages(pkg, repos="http://cran.fhcrc.org")
  } else {print(paste(pkg, " library already installed"))}
}