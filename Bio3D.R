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

#----------
#----------

#----------
# Specific selection (107, 108,109,124,128,129,132,198,214,290,291,294,313,316,317,320)
#----------
ca.inds_WT_res107 <- atom.select(pdb_WT, elety="CA", resno=107)
ca.inds_WT_res108 <- atom.select(pdb_WT, elety="CA", resno=108)
ca.inds_WT_res109 <- atom.select(pdb_WT, elety="CA", resno=109)
ca.inds_WT_res124 <- atom.select(pdb_WT, elety="CA", resno=124)
ca.inds_WT_res128 <- atom.select(pdb_WT, elety="CA", resno=128)
ca.inds_WT_res129 <- atom.select(pdb_WT, elety="CA", resno=129)
ca.inds_WT_res132 <- atom.select(pdb_WT, elety="CA", resno=132)
ca.inds_WT_res198 <- atom.select(pdb_WT, elety="CA", resno=198)
ca.inds_WT_res214 <- atom.select(pdb_WT, elety="CA", resno=214)
ca.inds_WT_res290 <- atom.select(pdb_WT, elety="CA", resno=290)
ca.inds_WT_res291 <- atom.select(pdb_WT, elety="CA", resno=291)
ca.inds_WT_res294 <- atom.select(pdb_WT, elety="CA", resno=294)
ca.inds_WT_res313 <- atom.select(pdb_WT, elety="CA", resno=313)
ca.inds_WT_res316 <- atom.select(pdb_WT, elety="CA", resno=316)
ca.inds_WT_res317 <- atom.select(pdb_WT, elety="CA", resno=317)
ca.inds_WT_res320 <- atom.select(pdb_WT, elety="CA", resno=320)

ca.inds_WT <- combine.select(ca.inds_WT_res107, ca.inds_WT_res108, ca.inds_WT_res109, ca.inds_WT_res124, ca.inds_WT_res128, 
                             ca.inds_WT_res129, ca.inds_WT_res132, ca.inds_WT_res198, ca.inds_WT_res214, ca.inds_WT_res290,
                             ca.inds_WT_res291, ca.inds_WT_res294, ca.inds_WT_res313, ca.inds_WT_res316, ca.inds_WT_res317,
                             ca.inds_WT_res320, operator="+")



ca.inds_Asn107Ile_res107 <- atom.select(pdb_Asn107Ile, elety="CA", resno=107)
ca.inds_Asn107Ile_res108 <- atom.select(pdb_Asn107Ile, elety="CA", resno=108)
ca.inds_Asn107Ile_res109 <- atom.select(pdb_Asn107Ile, elety="CA", resno=109)
ca.inds_Asn107Ile_res124 <- atom.select(pdb_Asn107Ile, elety="CA", resno=124)
ca.inds_Asn107Ile_res128 <- atom.select(pdb_Asn107Ile, elety="CA", resno=128)
ca.inds_Asn107Ile_res129 <- atom.select(pdb_Asn107Ile, elety="CA", resno=129)
ca.inds_Asn107Ile_res132 <- atom.select(pdb_Asn107Ile, elety="CA", resno=132)
ca.inds_Asn107Ile_res198 <- atom.select(pdb_Asn107Ile, elety="CA", resno=198)
ca.inds_Asn107Ile_res214 <- atom.select(pdb_Asn107Ile, elety="CA", resno=214)
ca.inds_Asn107Ile_res290 <- atom.select(pdb_Asn107Ile, elety="CA", resno=290)
ca.inds_Asn107Ile_res291 <- atom.select(pdb_Asn107Ile, elety="CA", resno=291)
ca.inds_Asn107Ile_res294 <- atom.select(pdb_Asn107Ile, elety="CA", resno=294)
ca.inds_Asn107Ile_res313 <- atom.select(pdb_Asn107Ile, elety="CA", resno=313)
ca.inds_Asn107Ile_res316 <- atom.select(pdb_Asn107Ile, elety="CA", resno=316)
ca.inds_Asn107Ile_res317 <- atom.select(pdb_Asn107Ile, elety="CA", resno=317)
ca.inds_Asn107Ile_res320 <- atom.select(pdb_Asn107Ile, elety="CA", resno=320)

ca.inds_Asn107Ile <- combine.select(ca.inds_Asn107Ile_res107, ca.inds_Asn107Ile_res108, ca.inds_Asn107Ile_res109, ca.inds_Asn107Ile_res124, 
                                    ca.inds_Asn107Ile_res128, ca.inds_Asn107Ile_res129, ca.inds_Asn107Ile_res132, ca.inds_Asn107Ile_res198,
                                    ca.inds_Asn107Ile_res214, ca.inds_Asn107Ile_res290, ca.inds_Asn107Ile_res291, ca.inds_Asn107Ile_res294,
                                    ca.inds_Asn107Ile_res313, ca.inds_Asn107Ile_res316, ca.inds_Asn107Ile_res317,ca.inds_Asn107Ile_res320,
                                    operator="+")


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

joint_rmsf = data.frame(WT = rf_WT,Asn107Ile = rf_Asn107Ile)

a = data.frame(RMSF = rf_WT, Category=rep('WT'), 
               ResID= c("107","108","109","124","128","129","132","198","214","290","291","294","313","316","317","320"))

b = data.frame(RMSF = rf_Asn107Ile, Category=rep('Asn107Ile'), 
               ResID= c("107","108","109","124","128","129","132","198","214","290","291","294","313","316","317","320"))

many_distros <- do.call('rbind', list(a,b))

p<-ggplot(data=many_distros, aes(x=ResID, y=RMSF, fill=Category)) +
  geom_bar(stat="identity", position=position_dodge())


p

#p+scale_fill_brewer(palette="Dark2")

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


