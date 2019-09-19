library(nlme)

hrvData	<- read.table('ChamberStudyData.txt', header = TRUE)
head(hrvData)

## primary analysis ##
outs	<- c('hr_nn', 'RMSSD', 'SDNN', 'HF', 'LF')

table2	<- matrix(0, nrow = length(outs), ncol = 6)
rownames(table2)	<- c('HR', 'RMSSD', 'SDNN', 'HF', 'LF')
colnames(table2)	<- c('Control', 'Flight', 'Rel. Diff.', 'Lower', 'Upper', 'P-value')

for(i in 1:length(outs)){
	
	hrvData$out	<- hrvData[,outs[i]]
	
	## fit gls ##
	model		<-	gls(log(out) ~ type*cruising + STARTH + as.factor(day) + as.factor(ID), data = hrvData, correlation = corAR1(form = ~ 1 | perExp), na.action = na.omit)

	## extract results ##
	cruise		<-	which(rownames(summary(model)$tTable) == "cruising")
	last		<-	dim(summary(model)$tTable)[1]
	intPval		<-	summary(model)$tTable[dim(summary(model)$tTable)[1],4]

	# Control
	contj		<-	summary(model)$tTable[cruise,1]
	stContj		<-	summary(model)$tTable[cruise,2]
	lContj		<-	contj - qnorm(0.975)*stContj
	uContj		<-	contj + qnorm(0.975)*stContj

	# Flight
	fltj		<-	summary(model)$tTable[cruise,1] + summary(model)$tTable[last,1]
	cruTypCov	<-	summary(model)$varBeta[last,cruise]
	vFltj		<-	summary(model)$tTable[last,2]^2 + stContj^2 + 2*cruTypCov
	stFltj		<-	sqrt(vFltj)
	lFltj		<-	fltj - qnorm(0.975)*stFltj
	uFltj		<-	fltj + qnorm(0.975)*stFltj

	inputj			<- c(contj, fltj)
	table2[i,1:2]	<- 100*(exp(inputj)-1)
	
	ssB12		<- vFltj + cruTypCov
	SigB		<- matrix(ssB12, nrow = 2, ncol = 2)
	SigB[1,1]	<- vFltj
	SigB[2,2]	<- stContj^2
	
	gB			<- matrix(c(exp(fltj), -exp(contj)), nrow = 2)
	varRD		<- t(gB)%*%SigB%*%gB
	
	table2[i,3]	<- 100*(sum(gB))
	table2[i,4]	<- 100*((sum(gB) - qnorm(0.975)*sqrt(varRD)))
	table2[i,5]	<- 100*((sum(gB) + qnorm(0.975)*sqrt(varRD)))
	zi			<- (sum(gB)^2)/varRD
	table2[i,6]	<- round(1-pchisq(zi, df = 1), 4)
}

# Table 2 #
table2

## time-varying analysis ##
outs	<- c('RMSSD', 'SDNN', 'hr_nn', 'LF', 'HF')

grids		<- vector(mode='list', length = length(outs))
timeVar		<- vector(mode='list', length = length(outs))

for(i in 1:length(outs)){
	
	hrvData$out			<- hrvData[,outs[i]]
	hrvData$typeCruis	<- hrvData$type*hrvData$cruising

	## fit gls ##
	model		<-	gls(log(out) ~ type*Texp + cruising + typeCruis + STARTH + as.factor(day) + as.factor(ID), data = hrvData, correlation = corAR1(form = ~ 1 | perExp), na.action = na.omit)
	modTable	<-	summary(model)$tTable
	varCov		<-	summary(model)$varBeta
	varNames	<-	rownames(modTable)
	b2Loc		<-	which(varNames == "Texp")
	b3Loc		<-	which(varNames == "cruising")
	b4Loc		<-	which(varNames == "typeCruis")
	b42Loc		<-	which(varNames == "type:Texp")
	
	# set up newdata for prediction
	gridTime	<- seq(min(hrvData$Texp), max(hrvData$Texp), length = 6*60)

	contEstT	<-	modTable[b3Loc,1] + gridTime*modTable[b2Loc,1]
	contVarT	<-	varCov[b3Loc,b3Loc] + (gridTime^2)*varCov[b2Loc,b2Loc] + 2*gridTime*varCov[b3Loc,b2Loc]
	contLowT	<-	contEstT - 1.96*sqrt(contVarT)
	contUppT	<-	contEstT + 1.96*sqrt(contVarT)
	contTimT	<- cbind(contEstT, contLowT, contUppT)

	fltEstT		<-	modTable[b3Loc,1] + modTable[b4Loc,1] + gridTime*(modTable[b2Loc,1] + modTable[b42Loc,1])
	fltVarT		<-	varCov[b3Loc,b3Loc] + varCov[b4Loc,b4Loc] + 2*varCov[b3Loc,b4Loc] + (gridTime^2)*(varCov[b2Loc,b2Loc] + varCov[b42Loc,b42Loc] + 2*varCov[b2Loc,b42Loc]) + 2*gridTime*(varCov[b3Loc,b2Loc] + varCov[b3Loc,b42Loc] + varCov[b4Loc,b2Loc] + varCov[b4Loc,b42Loc])
	fltLowT		<-	fltEstT - 1.96*sqrt(fltVarT)
	fltUppT		<-	fltEstT + 1.96*sqrt(fltVarT)
	fltTimT		<-	cbind(fltEstT, fltLowT, fltUppT)
	
	full			<-	cbind(contTimT, fltTimT)
	grids[[i]]		<- gridTime
	timeVar[[i]]	<- full
}

# Figure 2 #
mainNames	<- c("RMSSD","SDNN","HR","LF","HF")
ylims		<- vector('list', length = length(mainNames))
ylims[[1]]	<- c(-0.1, 0.4)
ylims[[2]]	<- c(-0.05, 0.45)
ylims[[3]]	<- c(-0.085,0.07)
ylims[[4]]	<- c(-0.05, 1)
ylims[[5]]	<- c(-0.16, 0.8)

layout(mat = matrix(c(1,1,2,2,3,3,
                      0,4,4,5,5,0), nrow = 2, byrow = TRUE))
for(i in 1:5){
	
	mini	<- min(timeVar[[i]], na.rm = TRUE)
	maxi	<- max(timeVar[[i]], na.rm = TRUE)
	xtloc	<- seq(mini, maxi, length = 5)
	labs	<- round(xtloc, 3)
	
	datt	<- timeVar[[i]]
	meansi	<- datt[, c('contEstT', 'fltEstT')]
	cont	<- datt[, c('contLowT', 'contUppT')]
	flt		<- datt[, c('fltLowT', 'fltUppT')]
	
	par(mai = rep(0.4,4))
	matplot(grids[[i]], meansi, type = 'l', lty = 1, pch = 16, ylab = 'Post - Pre', main = mainNames[i], xaxt = "n", cex.lab = 1.5, mgp = c(2,1,0), col=c(rgb(0, 0.75, 0), col=rgb(0, 0, 0.75)), axes=FALSE, ylim = ylims[[i]], xlab = 'Time (hours)',  cex.lab = 1)
	abline(h = 0, col = 'darkgray', lty = 2)
	axis(2, at=xtloc, labels = labs)
	axis(side = 1, at = 0:6)
	box()
	polygon(c(grids[[i]],rev(grids[[i]])), c(cont[,1],rev(cont[,2])), col=rgb(0, 0.75, 0, alpha = 0.25), border=NA)
	polygon(c(grids[[i]],rev(grids[[i]])), c(flt[,1],rev(flt[,2])), col=rgb(0, 0, 0.75, alpha = 0.25), border=NA)
}


