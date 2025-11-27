library(ggplot2)
library(reshape2)
library(grid)
library(ggseqlogo)

n=14

test_residuals <-read.table('/twins_residuals.txt')
residuals_n <- test_residuals[test_residuals$length == n,]
# residuals_n <- twins_residuals
residuals_n$residual[residuals_n$residual < 2] <- 0
residuals_n$residual[residuals_n$chi.p.value > 0.01] <- 0

x2 <- melt(residuals_n, c("aa", "celltype", "position"), "residual")
x2 <- x2[order(x2[,'aa'], x2[,'position'], -x2[,'value']),]
x2 <- x2[!duplicated(cbind(x2$aa, x2$position)),]

index <- x2$celltype == "pos"
x2$value[index] = -abs(x2$value[index])

x3 <- dcast(x2, aa ~ position)
x3[is.na(x3)] <- 0
x4 <- x3[,-1]
rownames(x4) <- x3[,1]

annotation <- data.frame(
  x = c(-1,-1),
  y = c(500,-500),
  label = c("HLA-positive enriched", "HLA-negative enriched")
)

#quartz( type="pdf", file="plots/sequencemotif.pdf", width=3, height=2.2, pointsize=8)
par( mar=c(3,4,0.2,0.5), family="sans" )

cs1 = make_col_scheme(chars=c('D', "E", 'R', 'H', 'K','M', 'W', 'L', 'V', 'F', 'A', 'I', 'X', 'P',
                              'Q', 'N', 'S', 'G', 'T', 'Y', 'C'), 
                      groups=c('acidic', 'acidic', 'basic', 'basic', 'basic', 'hydrophobic', 
                               'hydrophobic', 'hydrophobic', 'hydrophobic', 'hydrophobic', 'hydrophobic', 
                               'hydrophobic', 'hydrophobic', 'hydrophobic','neutral', 
                               'neutral', "polar", "polar", "polar", "polar", "polar"), 
                      cols=c('#d11141', '#d11141', '#ffc425', '#ffc425', '#ffc425', '#00aedb',
                                      '#00aedb', '#00aedb','#00aedb','#00aedb','#00aedb','#00aedb','#00aedb','#00aedb',
                                      '#f37735', '#f37735', '#00b159', '#00b159', '#00b159', '#00b159', 
                                      '#00b159'))
                                      
theme_replace(axis.line = element_line(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank())

limits <- max(abs(x4)) * 3.5
plt <- ggseqlogo(as.matrix(x4), method='custom', seq_type='aa', col_scheme=cs1) +
  theme(text = element_text(size = 12, face = "bold", family = "Helvetica Neue Light", color="black"),
        plot.title = element_text(size = 8, family = "Helvetica Neue Light", color="black"),
        plot.subtitle = element_text(size = 8, family = "Helvetica", color="black"),
        plot.caption = element_text(size = 8, family = "Helvetica", color="black"),
        axis.title=element_text(size=8,family = "Helvetica Neue Light", color="black"),
        axis.text=element_text(size=8,family = "Helvetica Neue Light", color="black")) +
  
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  xlab("position\n") +
  # ggtitle("MS01 vs. HD02") +
  ggtitle("TWINs") +
  annotate("text", x = c(-0.4, -0.8, -0.8), y = c(0, limits*0.75,limits*-0.75), 
           label = c(expression("amino acid\nenrichment"),
                     expression("HLA-negative"), expression("HLA-positive")), 
           size=6, family="Helvetica Neue Light", 
           angle = 90, lineheight = 1) +
  geom_segment(aes(x=1,xend=n,y=0,yend=0), size=0.3, show.legend = F) +
  geom_segment(aes(x = 0, y = limits, xend = 0, yend = limits * -1), size=0.3,
               arrow = arrow(length = unit(0.02, "npc"), ends = "both")) +
  geom_segment(aes(x = n, y = 0, xend = n, yend = 5), size=0.3,
               arrow = arrow(length = unit(0.02, "npc"))) +
  annotate("text", x = n+0.1, y = 2.5, 
           label = 5, 
           size=8, family="Helvetica Neue Light") + 
  theme(legend.title = element_text(size = 8), 
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.5, "lines"),
        legend.direction='horizontal',
        legend.box.margin = margin(t=10),
        legend.position = c(0.52, -0.21))

plt
#dev.off()
