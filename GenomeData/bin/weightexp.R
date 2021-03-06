suppressMessages(library(hash))
arg = commandArgs(T)



annotation_in <- read.table(arg[2], header = FALSE, sep = '\t', stringsAsFactors = FALSE) 


con <- file(arg[1], "r") #expressmatrix
line1 = readLines(con, n = 1)
line1 = unlist(strsplit(line1,split='\t'))
line1 = (line1[-(1:2)])
annohash = hash()

annotation_in = annotation_in[-1,]
.set(annohash, keys = annotation_in[,1], values = annotation_in[,2])
colallnum = as.numeric(unname(which(has.key(line1, annohash))))


write.table(t(c('number', 'name', 'weight')), file = arg[3], row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE, append = FALSE)

lineCnt = 0
output = matrix(ncol =3)
output = as.data.frame(output)
output[1,] = 1:3
output[,1] = as.character(output[,1])
output[,2] = as.character(output[,2])
output[,3] = as.character(output[,3])

while(1){
	oneline = readLines(con, n = 1)
	if(length(oneline) == 0){
		break
	}
	dataline = unlist(strsplit(oneline,split='\t'))
	names = dataline[(1:2)]
	expdata = dataline[-(1:2)]
	expdata = as.numeric(expdata)
	expdataselect = expdata[colallnum]

	names[1] = gsub(pattern = '\\.[^\\.]+$', replacement = '', names[1])
	#write.table(t(c(names, value)), file = arg[3], row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE, append = TRUE)
	output = rbind(output, c(names, mean(expdataselect)))
	lineCnt = lineCnt+1
}
close(con)
output = output[-1,]
output[,3] = as.numeric(output[,3])
#print(output[,3])
output = output[order(-output[,3]),]
outputnew = output[,3]
if (length(which(outputnew > 0.005)) >= 1) {
	if (length(which(outputnew <= 0.005)) >= 1) {
		zerorows = which(outputnew <= 0.005)
		points = quantile(outputnew[-zerorows], c(0.2, 0.4, 0.6, 0.8))
	}else{
		zerorows = c()
		points = quantile(outputnew[-zerorows], c(0.2, 0.4, 0.6, 0.8))
	}
	otherrows1 = which(outputnew <= points[1] & outputnew > 0.005)
	otherrows2 = which(outputnew <= points[2] & outputnew > points[1])
	otherrows3 = which(outputnew <= points[3] & outputnew > points[2])
	otherrows4 = which(outputnew <= points[4] & outputnew > points[3])
	otherrows5 = which(outputnew > points[4])
	outputnew[zerorows] <- 0
	outputnew[otherrows1] <- 1
	outputnew[otherrows2] <- 2
	outputnew[otherrows3] <- 3
	outputnew[otherrows4] <- 4
	outputnew[otherrows5] <- 5
}
count = c(zerorows, otherrows1, otherrows2, otherrows3, otherrows4, otherrows5)
length(count)
output[,3] = outputnew

write.table(output, file = arg[3], row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE, append = TRUE)

warnings()

