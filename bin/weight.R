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
	expdatabackground = expdata[-colallnum]
	ratio = (median(expdataselect) + 0.0001) / (median(expdatabackground)+0.0001)
	if (ratio >= 10000) {
		ratio = 10000
	}else if (ratio <= 0.0001) {
		ratio = 0.0001
	}
	
	if (max(expdataselect) >= 0.005){
	value = 1
	if (wilcox.test(expdataselect, expdatabackground, alternative = 'greater')$p.value <= 0.0005 ) {
		if (ratio <= 1) {
			value = 1
		}else{
			value = ratio
		}
	}else{
	if (wilcox.test(expdataselect, expdatabackground, alternative = 'less')$p.value <= 0.0005 ) {
		if (ratio >= 1) {
			value = 1
		}else{
			value = ratio
		}
	}else{
		value  = 1
	}

	}
	}else{
		value = 0
	}
	
	names[1] = gsub(pattern = '\\.[^\\.]+$', replacement = '', names[1])
	write.table(t(c(names, value)), file = arg[3], row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE, append = TRUE)

	lineCnt = lineCnt+1
}
close(con)
warnings()

