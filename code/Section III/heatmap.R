library(reshape2)
library(grid)

# source("../settings.R")
dev.off()
mat_chi <- read.table("/mat_chi_all.txt")
mat_cramer <- read.table("/mat_cramer_all.txt")

ms_n <-read.table("/pos_n_all.txt")
hd_n <-read.table("/neg_n_all.txt")

M <- mat_cramer
n.pos <- ncol(M)
l.from <- which.min(apply(M,2,function(x)all(is.na(x))))
l.to <- nrow(M)
Msk <- mat_chi>0.05
Msk[is.na(Msk)] <- FALSE
M[Msk] <- 0

cols <- colorRampPalette(c("lightgrey", "skyblue","yellow","orange", "tomato"))(100)

image(x=1:n.pos,y=l.from:l.to,as.matrix(M[,l.from:l.to]),bty="n", xlab="position",
      xaxt="n", yaxt="n", ylab="CDR3 length (aa)", xaxs="i", yaxs="i", 
      xlim=c(1,20), zlim=c(0,1),
      col=cols)
axis( 1, at=1:n.pos )
axis( 2, at=l.from:l.to, las=2 )
text(l.from:l.to, l.from:l.to, as.character(ms_n$x+hd_n$x),pos=4)
text(21,22,"n",xpd=TRUE, pos=4)
title(main = "HLA_DRB1*15", cex.main = 1.8, font.main= 4, col.main='black'
)

for (x in  0.5:22.5) {
  abline(v = x, col = "white", lty = 2) # Vertical lines at x positions
}

# Add horizontal lines (at each y)
for (y in 5.5:22.5) {
  abline(h = y, col = "white", lty = 2) # Horizontal lines at y positions
}


par( fig=c(0.85,0.87,0.3,0.6), mar=rep(0,4), new=TRUE )
image( x=1,z=matrix(seq(0,1,length.out=23),nrow=1,ncol=23), xlim=c(0,1),
       zlim=c(0,1), bty="o", xaxt="n", yaxt="n", col=cols,
       xlab="", ylab="")
mtext( "Cramer's V", 2 )
axis(4 ,at=0:1, labels=c("0/ns","1"), tick=FALSE, las=2, line=-.5)

#
