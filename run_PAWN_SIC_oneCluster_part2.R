library(ranger)

rm(list=ls())

#### ---------------------------------------
#### CONFIGURATION
#### ---------------------------------------
IT = 20
CASE0 = "OCS"

qq = seq(0.05,0.95,by=0.05)
qqlevel = abs(qq[1:(length(qq)/2)] - qq[length(qq):(length(qq)/2+1)]) 

NN = 500 ##total number of observations
Nc = 1 ## number of clusters
nmc = 1000 ## number of Monte-Carlo iterations

outfolder = paste0("./resu_EUR/",CASE0,"/")
resfolder = paste0("./resu_EUR/",CASE0,"/CLUST1/")

## load pre stored CASES
CASE.list = read.table("CASE.list_Nc1.txt")

N.list = c(250,400) ##number of observaitons within the cluster

#### ---------------------------------------
#### USEFUL FUNCTIONS
#### ---------------------------------------
CAfct <- function(y,lo,up){
	id <- NULL
	for (i in 1:length(y)){
		id[i] <- ifelse(y[i] >= lo[i] & y[i] <= up[i], 1, 0)
	}
	sum(id)/length(y)
}

#### ---------------------------------------
#### LOAD TRAIING / TEST DATASET FOR SETTINGS
#### ---------------------------------------
load(paste0("./data_EUR/",CASE0,"/OCSdata_regular.RData"))

x_test = OCSdata.regular[,-ncol(OCSdata.regular)]
y_test = OCSdata.regular[,ncol(OCSdata.regular)]

df.te = data.frame(
	X = x_test,
	OCS = y_test
)
load(paste0("./data_EUR/",CASE0,"/OCSdata_clust_train.RData"))
x_train = OCSdata.train[,-c(ncol(OCSdata.train))]

#### ---------------------------------------
#### LOOP
#### ---------------------------------------

for (N in N.list){

##################################################
### LOOP
for (it in 1:IT){

	CASE = CASE.list[it,1]
	nom = paste0(outfloder,"FIXE_CLUST_N",NN,"_nc",N,"NumC_",Nc,"_CASE",CASE) ## name where load data
	load(paste0(nom,"_",1,".RData"))
	nomS = paste0(resfolder,"FIXE_CLUST_N",NN,"_nc",N,"NumC_",Nc,"_CASE",CASE)#"_RSAMPL",it) ## name where store data

SS.crps = SS.is50 = SS.is90 = matrix(0,nrow(df.te),ncol(XR))

for (ii in 1:nrow(df.te)){

	print(ii)	
	AE = abs(M[,ii] - df.te$OCS[ii])

	## Approximated crps
	### https://search.r-project.org/CRAN/refmans/fabletools/html/distribution_accuracy_measures.html
	cov.pi = matrix(0,length(qq),nmc)
	for (kk in 1:(length(qq))){
		quant = QQ[,ii,kk]
		for (imc in 1:nmc){
			if (df.te$OCS[ii] < quant[imc]){
				cov.pi[kk,imc] = (-df.te$OCS[ii]+quant[imc])*(1-qq[kk])
			}
			if (quant[imc] <= df.te$OCS[ii]){
				cov.pi[kk,imc] = (df.te$OCS[ii]-quant[imc])*(qq[kk])
			}
		}
	}
	CRPS = 2*apply(cov.pi,2,sum)*0.05
	
	## Accuracy plot
	cov.pi = matrix(0,(length(qq)/2),nmc)
	for (kk in 1:(length(qq)/2)){
		lo = QQ[,ii,kk]
		up = QQ[,ii,length(qq) - kk]
		for (imc in 1:nmc) cov.pi[kk,imc] = CAfct(df.te$OCS[ii],lo[imc],up[imc])
	}
	CA = apply(cov.pi,2,mean)

## PI 1,19 ou 5,15
	Q1 = QQ[,ii,1]
	Q3 = QQ[,ii,19]
	IS90 = scoringutils:::interval_score(true_values=df.te$OCS[ii], lower=Q1,  upper=Q3,  interval_range=90)

	Q1 = QQ[,ii,5]
	Q3 = QQ[,ii,15]
	IS50 = scoringutils:::interval_score(true_values=df.te$OCS[ii], lower=Q1,  upper=Q3,  interval_range=50)

	### CRPS
	df = data.frame(XR,y=CRPS)
	for (j in 1:ncol(XR)){
		ll = unique(df[,j])
		SS0 = PVAL0 = NULL
		#plot(ecdf(df$y))			
		for (i in 1:length(ll)){
			ff = which(df[,j]==ll[i])
			#lines(ecdf(df$y[ff]),col=2)
			ks = ks.test(x=df$y[ff], y=df$y,alternative="g")
			SS0[i] = ks$statistic
			PVAL0[i] = ks$p.value
		}
		SS.crps[ii,j] = max(SS0)
	}

	### IS50
	df = data.frame(XR,y=IS50)
	for (j in 1:ncol(XR)){
		ll = unique(df[,j])
		SS0 = PVAL0 = NULL
		#plot(ecdf(df$y))			
		for (i in 1:length(ll)){
			ff = which(df[,j]==ll[i])
			#lines(ecdf(df$y[ff]),col=2)
			ks = ks.test(x=df$y[ff], y=df$y,alternative="g")
			SS0[i] = ks$statistic
			PVAL0[i] = ks$p.value
		}
		SS.is50[ii,j] = max(SS0)
	}

	### IS90
	df = data.frame(XR,y=IS90)
	for (j in 1:ncol(XR)){
		ll = unique(df[,j])
		SS0 = PVAL0 = NULL
		#plot(ecdf(df$y))			
		for (i in 1:length(ll)){
			ff = which(df[,j]==ll[i])
			#lines(ecdf(df$y[ff]),col=2)
			ks = ks.test(x=df$y[ff], y=df$y,alternative="g")
			SS0[i] = ks$statistic
			PVAL0[i] = ks$p.value
		}
		SS.is90[ii,j] = max(SS0)
	}

}

save(SS.crps,SS.is90,SS.is50,file=paste0(nomS,"PAWN_",1,"CRPS.RData"))

}##it

}#}}### cases
