## load package for qRF modelling
library(ranger)

rm(list=ls())

#### ---------------------------------------
#### SET SEED
#### ---------------------------------------
set.seed(10)

## Number of iterations
IT = 1

SEEDr = sample(1:10000,IT,replace=FALSE)

#######################################################
##### USEFUL FUNCTIONS
#######################################################
source("./utils/sim_utils.R")

#### ---------------------------------------
#### DATA
#### ---------------------------------------
## Define here the directory to find the data from
## https://zenodo.org/records/6513429
infolder  <- "../../clusteredData/data"

## CHOOSE OCS or AGB
CASE0 = "OCS"

## directory where to store training and test data
outfolder <- paste0("./data_EUR/",CASE0,"/")

## directory where to store results
resfolder <- paste0("./resu_EUR/",CASE0,"/")

##number of observations
NN = 500

##number of Monte Carlo simulations
nmc = 10

##list of random clusters
CASE.list = read.table("CASE.list_Nc1.txt")

## Choose here the number of observations within the cluster
N.list = c(400)

for (N in N.list){## loop on the number of observations

for (it in 1:IT){## loop on number of iterations

## CASE = geostratum is considered a cluster
CASE = CASE.list[it,1]

## set the name for storage
nom = paste0(resfolder,"/FIXE_CLUST_N",NN,"_nc",N,"NumC_",1,"_CASE",CASE)

## generate
SEED = SEEDr[it]
source("./utils/run_extract_oneCluster_Europe.R")

## define the training dataset
load(paste0(outfolder,"/OCSdata_clust_train.RData"))
x_train = OCSdata.train[,-c(ncol(OCSdata.train))]
y_train = OCSdata.train[,c(ncol(OCSdata.train))]
df.tr = data.frame(
  X = x_train,
  OCS = y_train
)
names(df.tr) = c(names(x_train),"OCS")
df.tr = na.omit(df.tr)

## define the index of the observatiosn within the cluster
id_clust1 = 1:nnc
id_ext = (nnc+1):(nrow(df.tr))
id = c(
	rep(1,length(id_clust1)),
	rep(2,length(id_ext))
	)

## define the test dataset
load(paste0(outfolder,"/OCSdata_regular.RData"))
x_test = OCSdata.regular[,-ncol(OCSdata.regular)]
y_test = OCSdata.regular[,ncol(OCSdata.regular)]

df.te = data.frame(
  X = x_test,
  OCS = y_test
)
names(df.te) = c(names(x_train),"OCS")
df.te = na.omit(df.te)

## define the geographical coordinates
XY.te = df.te[,c(1,2)]

#### ---------------------------------------
#### correlation analysis
#### ---------------------------------------
C = NULL
if (CASE0 == "AGB") quel = 13
if (CASE0 == "OCS") quel = 10
ou = (1:(ncol(df.tr)-1))#[-c(1,2,10)]
for (i in ou){
	if (i != quel) C[i] = cor(df.tr[,i],df.tr$OCS, method="kendall")
	if (i == quel){
		xo = as.factor(df.tr[,i])
		yo = df.tr$OCS
		df0 = data.frame(xo,yo)
		C[i] = sqrt(summary(lm(yo ~ xo,df0))$r.squared)
	}
}
oo = order(abs(C),  decreasing = T)

#### ---------------------------------------
#### Spatial proxies
#### ---------------------------------------
add_dist <- function(df.tr,XY,mini,maxi){
	df.tr2 = df.tr
	corner = matrix(0,5,2)
	corner[1, ] = c(mini[1],mini[2])
	corner[2, ] = c(mini[1],maxi[2])
	corner[3, ] = c(maxi[1],mini[2])
	corner[4, ] = c(maxi[1],maxi[2])
	corner[5, ] = c((mini[1]+mini[2])/2,(maxi[1]+maxi[2])/2)

	mydist<-function(u,v=rep(0,2)){
		sqrt((u[1]-v[1])^2+(u[2]-v[2])^2)
	}
	DD = matrix(0,nrow(df.tr),5)
	for (ii in 1:5) DD[,ii] = apply(XY,1,mydist,v=corner[ii,])
	df.tr2$D1 = DD[,1]
	df.tr2$D2 = DD[,2]
	df.tr2$D3 = DD[,3]
	df.tr2$D4 = DD[,4]
	df.tr2$D5 = DD[,5]
	return(df.tr2)
}

mini = c(min(x_test[,1]),min(x_test[,2]))
maxi = c(max(x_test[,1]),max(x_test[,2]))

## geograhical coordinates
XY = df.tr[,c(1,2)]

#### ---------------------------------------
#### LOOP
#### ---------------------------------------
ncv = 10
## Define here the levels / values of the random variables
NCV = matrix(sample(1:nrow(df.tr),nrow(df.tr),replace=F),ncol = ncv)
MTRY = c("sqrt","1/3","max")
NVAR = seq(2,round((ncol(df.tr)-1)/2))
NODE = seq(1,10)
SAMPL = c("variance", "extratrees", "maxstat")
ERR = c("Y","N")
SPA = c("Y","N")

## Sample here
rtr = sample(1:ncv,nmc,replace=T)
rmtry = sample(MTRY,nmc,replace=T)
rnode = sample(NODE,nmc,replace=T)
rerr = sample(ERR,nmc,replace=T)
rvar = sample(NVAR,nmc,replace=T)
rsampl = sample(SAMPL,nmc,replace=T)
rspa = sample(SPA,nmc,replace=T)

M = Q1 = Q3 = matrix(0,nmc,nrow(df.te))
qq = seq(0.05,0.95,by=0.05)
qqlevel = abs(qq[1:(length(qq)/2)] - qq[length(qq):(length(qq)/2+1)]) 
QQ = array(0,c(nmc,nrow(df.te),length(qq)))

for (ii in 1:nmc){

	print(ii)

	dfB = df.tr[-NCV[,rtr[ii]],]
	dfB = dfB[,c(oo[1:rvar[ii]],ncol(dfB))]
	XYB = XY[-NCV[,rtr[ii]],]

	id2 = id[-NCV[,rtr[ii]]]

	if (rerr[ii] == "N"){
		w = rep(1,nrow(dfB))
	}else{
		w = rep(1,nrow(dfB))
		w[(id2==2)] = 1/4
	}

	if (rspa[ii] == "Y"){
	
		dfB2 = add_dist(dfB,XYB,mini,maxi)

		if (rmtry[ii] == "sqrt"){
			rmtry2 = sqrt(ncol(dfB2)-1)
		}else if (rmtry[ii] == "max"){
			rmtry2 = ncol(dfB2)-1
		}else if (rmtry[ii] == "1/3"){
			rmtry2 = round((ncol(dfB2)-1)/3)
		}

		modelB <- ranger(
			formula = OCS~., 
 			data = dfB2,
  			mtry=rmtry2,
  			min.node.size = rnode[ii],
 			num.trees = 1000,
  			respect.unordered.factors = 'order',
  			quantreg = TRUE,
			splitrule = rsampl[ii],
			case.weights = w 
			)

		df.te2 = add_dist(df.te[,oo[1:rvar[ii]]],XY.te,mini,maxi)
		names(df.te2) = names(dfB2)[-which(names(dfB2) %in% "OCS")]

		PREB = predict(modelB,df.te2,type = "quantiles", quantiles = qq)
		M[ii,] = PREB$predictions[,10]##median
		Q1[ii,] = PREB$predictions[,1]##percentile at 0.05
		Q3[ii,] = PREB$predictions[,19]##percentile at 0.95
		QQ[ii,,] = PREB$predictions##all percentiles

	}else{
	
		dfB2 = dfB

		if (rmtry[ii] == "sqrt"){
			rmtry2 = sqrt(ncol(dfB2)-1)
		}else if (rmtry[ii] == "max"){
			rmtry2 = ncol(dfB2)-1
		}else if (rmtry[ii] == "1/3"){
			rmtry2 = round((ncol(dfB2)-1)/3)
		}

		modelB <- ranger(
			formula = OCS~., 
 			data = dfB2,
  			mtry = rmtry2,
  			min.node.size = rnode[ii],
 			num.trees = 1000,
  			respect.unordered.factors = 'order',
  			quantreg = TRUE,
			splitrule = rsampl[ii],
			case.weights = w 
			)

		df.te2 = df.te[,oo[1:rvar[ii]]]
		names(df.te2) = names(dfB2)[-ncol(dfB2)]

		PREB = predict(modelB,df.te2,type = "quantiles", quantiles = qq)
		M[ii,] = PREB$predictions[,10]##median
		Q1[ii,] = PREB$predictions[,1]##percentile at 0.05
		Q3[ii,] = PREB$predictions[,19]##percentile at 0.95
		QQ[ii,,] = PREB$predictions##all percentiles

	}

}

XR = data.frame(
	tr = rtr,
	nvar = rvar,
	sampl = rsampl,
	mtry = rmtry,
	node = rnode,
	err = rerr,
	spa =rspa
	)

save(XR,M,Q1,Q3,QQ,file=paste0(nom,"_",1,".RData"))##save

}##it

}### cases
############################################################
