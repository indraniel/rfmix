for(chr in 1:1){
adm = read.table(paste("admChr", chr, ".bgl.phased",sep=""), header=T, stringsAsFactors=F)
adm = adm[, -c(1:2)]
admT = t(adm)
numAdm = nrow(admT)

ceu = read.table(paste("ceuRefChr", chr, "bgl.phased", sep=""), header=T, stringsAsFactors=F)
ceu = ceu[, -c(1:2)]
ceuT = t(ceu)
numCeu = nrow(ceuT)

yri = read.table(paste("yriRefChr", chr, "bgl.phased", sep=""), header=T, stringsAsFactors=F)
yri = yri[, -c(1:2)]
yriT = t(yri)
numYri = nrow(yriT)

jpt = read.table(paste("jpt+chbRefChr", chr, "bgl.phased", sep=""), header=T, stringsAsFactors=F)
jpt = jpt[, -c(1:2)]
jptT = t(jpt)
numJpt = nrow(jptT)

all = rbind(admT, ceuT, yriT, jptT)

binarize = function(snp){
	asFactor = factor(snp)
	allele1 = levels(asFactor)[1]
	allele2 = levels(asFactor)[2]
	result = rep(NA, length(snp))
	result[snp == allele1] = '0'
	result[snp == allele2] = '1'
	return(result)
}

all01 = apply(all, 2, binarize)

# Transpose back
alleles = t(all01)

write.table(alleles, paste("../alleles", chr, ".txt", sep=""), quote=F, sep="", col.names=F, row.names=F)
}

# Classes
classes = c(rep(0,numAdm), rep(1, numCeu), rep(2, numYri), rep(3,numJpt))
classesM = matrix(classes, nrow=1)
write.table(classesM, "../classes.txt", quote=F, row.names=F, col.names=F)

