##############
#Bray curtis
##############
library(ape)
library(ggplot2)
library(vegan)

#read in abundance table, transpose
taxa <- read.table("biom.txt", header=T, sep="\t", row.names=1)
taxa.t <- t(taxa)

#read in metadata file
taxa.md <- read.table("map.txt", header=T, sep="\t", row.names=1)

#merge metadata with biom file
merge <- transform(merge(taxa.t, taxa.md, by=0), row.names=Row.names)

#get braycurtis distance matrix on your count table (first need to determine which rows you want to pull from your merged dataframe)
n <- ncol(taxa.t) - 1

bcdist <- vegdist(merge[,2:n], method="bray")

#get your coordinates from your dissimilarity matrix
pcoa.taxa <- pcoa(bcdist, correction="none")

#what proportion of variance is explained by each PC?
modSource <- capscale(bcdist ~1)
eigSource <- eigenvals(modSource)
eigSource/sum(eigSource)

#plot your data
pdf("pcoa_bcdist.pdf")
qplot(pcoa.taxa$vectors[,1], pcoa.taxa$vectors[,2], color=merge$Phyl_Group, shape=merge$Genus, ylab="PC2 (5%)", xlab="PC1 (6%)", ylim=c(-0.6, 0.5), xlim=c(-0.6, 0.5)) + scale_shape_manual(values=1:13) + theme_classic() + coord_fixed()
dev.off()

###############
#Unifrac
###############
library("phyloseq");packageVersion("phyloseq")

#convert biom to matrix format, read in tree
taxa <- as.matrix(taxa)
tree <- read.tree("tree.tre")

#generate your phyloseq object
otutable <- otu_table(taxa, taxa_are_rows=T)
map <- sample_data(data.frame(taxa.md))
physeq <- phyloseq(otutable)
physeq <- merge_phyloseq(physeq, tree, map)

#plot data using unifrac distance
#dev.new()
pdf("pcoa_unifrac.pdf")
plot_ordination(physeq, ordinate(physeq, "PCoA", "unifrac", weighted=T), color="Phyl_Group", shape="Genus") + theme_classic() + coord_fixed() + scale_shape_manual(values=1:13) + scale_x_continuous(limits = c(-0.3, 0.55)) + scale_y_continuous(limits = c(-0.3, 0.55))
dev.off()

#if you want the matrix itself
#UniFrac(physeq, weighted=T, normalized=F)








