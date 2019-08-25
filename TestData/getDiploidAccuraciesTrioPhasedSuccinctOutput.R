getGenotype = function(row){
	row = row[order(row)]
	return(paste(row, collapse=""))
}

solution = read.table(paste("admSolutionChr1.bgl.phased"), stringsAsFactors=F, header=T)
solution = solution[, -c(1:2)]
numAdm = ncol(solution)/2

output = as.matrix(read.table("../outputTrioPhased.0.Viterbi.txt", stringsAsFactors=F))

numSnps = nrow(solution)
totalLocs = numAdm*numSnps
newTotal = 0

snpsPerWindow = read.table("../outputTrioPhased.0.SNPsPerWindow.txt", stringsAsFactors=F)[,1]
numWindows = length(snpsPerWindow)
outputLong = matrix(NA, nrow=numSnps, ncol=numAdm*2)
currentSNP = 1
for(windex in 1:numWindows){
	for(snpIndex in 1:snpsPerWindow[windex]){
		outputLong[currentSNP,] = output[windex,]
		currentSNP = currentSNP + 1
	}
}

output = outputLong

for(adm in 1:numAdm){
	#cat(adm, " ")
	called = output[,c(1+2*(adm-1), 2*(adm-1)+2)]
	true = solution[,c(1+2*(adm-1), 2*(adm-1)+2)]
	
	calledGen = apply(called,1, getGenotype)
	trueGen = apply(true,1, getGenotype)
	
	newAccuracy = length(which(calledGen == trueGen))
	newTotal = newTotal + newAccuracy
	
	cat(newAccuracy/numSnps, "\t")
}

cat("\nMean Accuracy", newTotal/totalLocs, "\n", sep="\t")


