library(bio3d)
library(ggplot2)
library(ggpubr)

memory.limit(800000)
#--------
#Loading files
#-------
pdb_WT=read.pdb("WT_NPRS1_popcwi.pdb")
dcd_WT= read.dcd("WT_nprs1_popcwieq-04.dcd")

pdb_Asn107Ile= read.pdb("Asn107Ile_NPRS1_popcwi.pdb")
dcd_Asn107Ile <- read.dcd("Asn107Ile_nprs1_popcwieq-04.dcd")


print(pdb_WT)
print(pdb_Asn107Ile)

attributes(pdb_WT)
attributes(pdb_Asn107Ile)

print(dcd_WT)
print(dcd_Asn107Ile)


#----------
# All aa
#-----------
ca.inds_WT <- atom.select(pdb_WT, elety="CA")
ca.inds_Asn107Ile <- atom.select(pdb_Asn107Ile, elety="CA")

xyz_WT <- fit.xyz(fixed=pdb_WT$xyz, mobile=dcd_WT,
                  fixed.inds=ca.inds_WT$xyz,
                  mobile.inds=ca.inds_WT$xyz)

xyz_Asn107Ile<- fit.xyz(fixed=pdb_Asn107Ile$xyz, mobile=dcd_Asn107Ile,
               fixed.inds=ca.inds_Asn107Ile$xyz,
               mobile.inds=ca.inds_Asn107Ile$xyz)

dim(xyz_WT) == dim(dcd_WT)
dim(xyz_Asn107Ile) == dim(dcd_Asn107Ile)
#--------
#RMSD
#--------
rd_WT <- rmsd(xyz_WT[1,ca.inds_WT$xyz], xyz_WT[,ca.inds_WT$xyz])
rd_Asn107Ile <- rmsd(xyz_Asn107Ile[1,ca.inds_Asn107Ile$xyz], xyz_Asn107Ile[,ca.inds_Asn107Ile$xyz])

plot(rd_WT, typ="l", ylab="RMSD", xlab="Time ns", main="RMSD", col="black", ylim = c(0,3.5))
points(rd_Asn107Ile, typ="l", col="red", lty=2, lwd=1)
points(lowess(rd_WT), typ="l", col="blue", lty=2, lwd=2)
points(lowess(rd_Asn107Ile), typ="l", col="blue", lty=2, lwd=2)

legend(800,1, legend=c("WT","Asn107Ile"), col=c("black","red"), lty=1:1, cex=0.8)
#------

hist(rd_WT, breaks=40, freq=FALSE, main="RMSD Histogram - WT", xlab="RMSD")
lines(density(rd_WT), col="gray", lwd=3)

hist(rd_Asn107Ile, breaks=40, freq=FALSE, main="RMSD Histogram - Asn107Ile", xlab="RMSD")
lines(density(rd_Asn107Ile), col="gray", lwd=3)


summary(rd_WT)
summary(rd_Asn107Ile)
#----------
plot_multi_histogram <- function(df, feature, label_column) {
        plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
                geom_histogram(alpha=0.7, position="identity", aes(y = ..density..), color="black") +
                geom_density(alpha=0.7) +
                geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", linetype="dashed", size=1) +
                labs(x=feature, y = "Density")
        plt + guides(fill=guide_legend(title=label_column))
}


a = data.frame(RMSD=rd_WT, Category=rep('WT'))
b = data.frame(RMSD=rd_Asn107Ile, Category=rep("Asn107Ile"))
many_distros <- do.call('rbind', list(a,b))


options(repr.plot.width = 20, repr.plot.height = 2)
plot_multi_histogram(many_distros, 'RMSD', 'Category')
#----------
rf_WT <- rmsf(xyz_WT[,ca.inds_WT$xyz])
rf_Asn107Ile <- rmsf(xyz_Asn107Ile[,ca.inds_Asn107Ile$xyz])

plot(rf_WT, ylab="RMSF", xlab="Residue Position", typ="l", main="RSMF", col="black", ylim = c(0,5))
points(rf_Asn107Ile, typ="l", col="red", lty=1, lwd=2)

legend(300, 4, legend=c("WT","Asn107Ile"), col=c("black", "red"), lty=1:1, cex=0.8)

#-------
pc_WT <- pca.xyz(xyz_WT[,ca.inds_WT$xyz])
plot(pc_WT, col=bwr.colors(nrow(xyz_WT)))

pc_Asn107Ile <- pca.xyz(xyz_Asn107Ile[,ca.inds_Asn107Ile$xyz])
plot(pc_Asn107Ile, col=bwr.colors(nrow(xyz_Asn107Ile)))




a = data.frame(PCA1=pc_WT$z[,1],PCA2=pc_WT$z[,2], Category=rep('WT'))
b = data.frame(PCA1=pc_Asn107Ile$z[,1],PCA2=pc_Asn107Ile$z[,2], Category=rep('Asn107Ile'))
many_distros <- do.call('rbind', list(a,b))

# Scatter plot colored by groups 
sp <- ggscatter(many_distros, x = "PCA1", y = "PCA2",
                color = "Category", palette = "hallmarks_light",
                size = 1, alpha = 0.5)+
        border()                                         
# Marginal density plot of x (top panel) and y (right panel)
xplot <- ggdensity(many_distros, "PCA1", fill = "Category",
                   palette = "hallmarks_light")
yplot <- ggdensity(many_distros, "PCA2", fill = "Category", 
                   palette = "hallmarks_light")+
        rotate()
# Cleaning the plots
yplot <- yplot + clean_theme() 
xplot <- xplot + clean_theme()
# Arranging the plot
ggarrange(xplot, NULL, sp, yplot, 
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(2, 1), heights = c(1, 2),
          common.legend = TRUE)

#-----
#Rg
#-----
rg_WT <- rgyr(xyz_WT)
rg_Asn107Ile <- rgyr(xyz_Asn107Ile)

plot(rg_WT, typ="l", ylab="Rg", xlab="Frame No.", main="Radius of Gyration", col="black")
points(rg_Asn107Ile, typ="l", col="red", lty=2, lwd=1)
points(lowess(rg_WT), typ="l", col="blue", lty=2, lwd=2)
points(lowess(rg_Asn107Ile), typ="l", col="blue", lty=2, lwd=2)

legend(900,63.2, legend=c("WT","Asn107Ile"), col=c("black", "red"), lty=1:1, cex=0.8)

#----------
#----------

#----------
# 100-330
#----------
ca.inds_WT <- atom.select(pdb_WT, elety="CA", resno=100:330)
ca.inds_Asn107Ile <- atom.select(pdb_Asn107Ile, elety="CA", resno=100:330)

xyz_WT <- fit.xyz(fixed=pdb_WT$xyz, mobile=dcd_WT,
                  fixed.inds=ca.inds_WT$xyz,
                  mobile.inds=ca.inds_WT$xyz)

xyz_Asn107Ile<- fit.xyz(fixed=pdb_Asn107Ile$xyz, mobile=dcd_Asn107Ile,
                        fixed.inds=ca.inds_Asn107Ile$xyz,
                        mobile.inds=ca.inds_Asn107Ile$xyz)

dim(xyz_WT) == dim(dcd_WT)
dim(xyz_Asn107Ile) == dim(dcd_Asn107Ile)

#--------
#RMSD
#--------
rd_WT <- rmsd(xyz_WT[1,ca.inds_WT$xyz], xyz_WT[,ca.inds_WT$xyz])
rd_Asn107Ile <- rmsd(xyz_Asn107Ile[1,ca.inds_Asn107Ile$xyz], xyz_Asn107Ile[,ca.inds_Asn107Ile$xyz])

plot(rd_WT, typ="l", ylab="RMSD", xlab="Time ns", main="RMSD (100-330)", col="black", ylim = c(0,3.5))
points(rd_Asn107Ile, typ="l", col="red", lty=2, lwd=1)
points(lowess(rd_WT), typ="l", col="blue", lty=2, lwd=2)
points(lowess(rd_Asn107Ile), typ="l", col="blue", lty=2, lwd=2)

legend(700,1, legend=c("WT","Asn107Ile"), col=c("black","red"), lty=1:1, cex=0.8)

#----------
plot_multi_histogram <- function(df, feature, label_column) {
  plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
    geom_histogram(alpha=0.7, position="identity", aes(y = ..density..), color="black") +
    geom_density(alpha=0.7) +
    geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", linetype="dashed", size=1) +
    labs(x=feature, y = "Density")
  plt + guides(fill=guide_legend(title=label_column))
}


a = data.frame(RMSD=rd_WT, Category=rep('WT'))
b = data.frame(RMSD=rd_Asn107Ile, Category=rep("Asn107Ile"))
many_distros <- do.call('rbind', list(a,b))


options(repr.plot.width = 20, repr.plot.height = 2)
plot_multi_histogram(many_distros, 'RMSD', 'Category')
#----------
rf_WT <- rmsf(xyz_WT[,ca.inds_WT$xyz])
rf_Asn107Ile <- rmsf(xyz_Asn107Ile[,ca.inds_Asn107Ile$xyz])

plot(rf_WT, ylab="RMSF", xlab="Residue Position", typ="l", main="RSMF", col="black", ylim = c(0,5))
points(rf_Asn107Ile, typ="l", col="red", lty=1, lwd=2)

legend(50, 4, legend=c("WT","Asn107Ile"), col=c("black", "red"), lty=1:1, cex=0.8)

#-------
pc_WT <- pca.xyz(xyz_WT[,ca.inds_WT$xyz])
plot(pc_WT, col=bwr.colors(nrow(xyz_WT)))

pc_Asn107Ile <- pca.xyz(xyz_Asn107Ile[,ca.inds_Asn107Ile$xyz])
plot(pc_Asn107Ile, col=bwr.colors(nrow(xyz_Asn107Ile)))




a = data.frame(PCA1=pc_WT$z[,1],PCA2=pc_WT$z[,2], Category=rep('WT'))
b = data.frame(PCA1=pc_Asn107Ile$z[,1],PCA2=pc_Asn107Ile$z[,2], Category=rep('Asn107Ile'))
many_distros <- do.call('rbind', list(a,b))

# Scatter plot colored by groups 
sp <- ggscatter(many_distros, x = "PCA1", y = "PCA2",
                color = "Category", palette = "hallmarks_light",
                size = 1, alpha = 0.5)+
  border()                                         
# Marginal density plot of x (top panel) and y (right panel)
xplot <- ggdensity(many_distros, "PCA1", fill = "Category",
                   palette = "hallmarks_light")
yplot <- ggdensity(many_distros, "PCA2", fill = "Category", 
                   palette = "hallmarks_light")+
  rotate()
# Cleaning the plots
yplot <- yplot + clean_theme() 
xplot <- xplot + clean_theme()
# Arranging the plot
ggarrange(xplot, NULL, sp, yplot, 
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(2, 1), heights = c(1, 2),
          common.legend = TRUE)


#-----
#Rg
#-----
rg_WT <- rgyr(xyz_WT)
rg_Asn107Ile <- rgyr(xyz_Asn107Ile)

plot(rg_WT, typ="l", ylab="Rg", xlab="Frame No.", main="Radius of Gyration", col="black")
points(rg_Asn107Ile, typ="l", col="red", lty=2, lwd=1)
points(lowess(rg_WT), typ="l", col="blue", lty=2, lwd=2)
points(lowess(rg_Asn107Ile), typ="l", col="blue", lty=2, lwd=2)

legend(900,63.2, legend=c("WT","Asn107Ile"), col=c("black", "red"), lty=1:1, cex=0.8)