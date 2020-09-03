args <- commandArgs(TRUE)
logratiofile <- args[1]
outfile<- args[2]

lr.df = read.table(logratiofile,header=TRUE,sep='\t')
lr.h <- colnames(lr.df)
sampleName = strsplit(logratiofile,'\\.')[[1]][1]
print(sampleName)

# RATIO col
lr.col <- NULL
flag <- 0
for (i in lr.h){
	flag <- flag + 1
	if (i == "RATIO"){
		lr.col <- flag
	}
}

lr.df[,lr.col][is.infinite(lr.df[,lr.col])] <- 1
lr.df[,lr.col][is.na(lr.df[,lr.col])] <- 1
cn.val <- 2 * (lr.df[,lr.col])

# add CN col
lr.df$CN <- cn.val

# CN range
CN.range <- range(lr.df$CN)

y.adj <- CN.range[2] * 0.02

xAxisMaxValue <- nrow(lr.df)
smp <- strsplit(outfile,'\\/')[[1]]
smp.len <- length(smp)
print(smp.len)
smp.pfx <- strsplit(smp[smp.len],'\\.')[[1]][1]
print(smp.pfx)

pdf(file = outfile, width = 15, height = 8)
par(mar = c(5,5,5,10), las=1)
plot(NULL, xlim=c(0, xAxisMaxValue), ylim=c(0,4), xlab="", ylab="Copy Number", xaxt = "n", yaxt = "n", xaxs = "i", col="dimgrey", cex=2,main=smp.pfx,cex.lab=2,cex.main=3)

gene <- unique(lr.df$GENE)

for (i in gene){
    if (i == 'BRCA1'){
        index <- which(lr.df$GENE==i)
        yVal <- lr.df$CN[index] # get CN
        points(index, yVal, pch=20, cex=1.5, col="blue")
    }
    else if (i == 'BRCA2'){
        index <- which(lr.df$GENE==i)
        yVal <- lr.df$CN[index]
        points(index, yVal, pch=20, cex=1.5, col="red")
    }
    else {
        index <- which(lr.df$GENE==i)
        yVal <- lr.df$CN[index]
        points(index, yVal, pch=20, cex=1, col="dimgrey")
    }
}

chrs <- 1:22
for (i in chrs){
    index <- which(lr.df$CHR==i)
    index.last <- index[length(index)]
    abline(v=index.last+0.5, lty=1, lwd=1, col="darkgray")

}

yCN <- c(0,1,2,3,4)
axis(2, at = yCN, label = yCN, cex.axis = 2)

for (i in yCN){
	abline(h=i, lty=4, lwd=1, col="darkgray")
}

abline(h=2, lty=1, lwd=4, col="blue")

for (i in chrs){
	index <- which(lr.df$CHR==i)
	rect(index[1],4,index[length(index)],5,col="grey") # xleft/ybottom/xright/ytop
	text(median(index), 4.09, labels=i, cex=1.5) # x,y,lab,cex
}


legend.x <- nrow(lr.df) + 20
legend.y <- 2.8

legend(legend.x, legend.y, legend=c("BRCA1","BRCA2","others"), title="", pt.cex=2, cex=2, col=c('blue',"red","dimgrey"), pch=20, ncol=1, bty="n",xpd=TRUE)


dev.off()
print(warnings())

