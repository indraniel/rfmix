getGenotype = function(row){
	row = row[order(row)]
	return(paste(row, collapse=""))
}

solution = read.table(paste("admSolutionChr1.bgl.phased"), stringsAsFactors=F, header=T)
solution = solution[, -c(1:2)]
numAdm = ncol(solution)/2

output = read.table("../outputPopPhased.0.Viterbi.txt", stringsAsFactors=F)

numSnps = nrow(solution)
totalLocs = numAdm*numSnps
newTotal = 0

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


