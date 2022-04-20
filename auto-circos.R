library(circlize)
library(StructuralVariantAnnotation)

#list of file names

##############################
### ADD TXT FILE NAME HERE ###
##############################
## TXT FILE AND VCFS NEED TO BE IN THE SAME DIRECTORY AND FOLDER AS THIS SCRIPT ##
file_names <- read.table("[txtfile with sample file names here]")


for (row in 1:nrow(file_names)){

file_name <- file_names[row,1]
pdf_name <- paste(file_name, ".pdf", sep="")

#####################
### ADD PATH HERE ### 
#####################
pdf_location <- paste("[add path here]", pdf_name, sep="")


#converts vcf to GRanges
gRange_raw_data <- pairs2breakpointgr(rtracklayer::import(file_name))


#del -> pairs
gRange_del_data <- subset(gRange_raw_data, gRange_raw_data$NA. == "DEL")
head(gRange_del_data)
pairs_del <- breakpointgr2pairs(gRange_del_data)
head(pairs_del)


#ins -> pairs
gRange_ins_data <- subset(gRange_raw_data, gRange_raw_data$NA. == "INS")
pairs_ins <- breakpointgr2pairs(gRange_ins_data)

#dup -> pairs
gRange_dup_data <- subset(gRange_raw_data, gRange_raw_data$NA. == "DUP")
pairs_dup <- breakpointgr2pairs(gRange_dup_data)

#bnd -> pairs
gRange_bnd_data <- subset(gRange_raw_data, gRange_raw_data$NA. == "BND")
pairs_bnd <- breakpointgr2pairs(gRange_bnd_data)
head(pairs_bnd)


#build the circos plot

pdf(file = pdf_location,
    width = 8.50,
    height = 9)

layout(matrix(c(1,0,2,3,4,5), 3, 2, byrow = TRUE))

#whole circos
circos.initializeWithIdeogram(species = "hg38")
circos.genomicLink(as.data.frame(S4Vectors::first(pairs_del)), as.data.frame(S4Vectors::second(pairs_del)), col = "darkmagenta")
circos.genomicLink(as.data.frame(S4Vectors::first(pairs_ins)), as.data.frame(S4Vectors::second(pairs_ins)), col = "black")
circos.genomicLink(as.data.frame(S4Vectors::first(pairs_dup)), as.data.frame(S4Vectors::second(pairs_dup)), col = "darkorange")
circos.genomicLink(as.data.frame(S4Vectors::first(pairs_bnd)), as.data.frame(S4Vectors::second(pairs_bnd)), col = "forestgreen")

#del
circos.initializeWithIdeogram(species = "hg38")
circos.genomicLink(as.data.frame(S4Vectors::first(pairs_del)), as.data.frame(S4Vectors::second(pairs_del)), col = "darkmagenta")

#ins
circos.initializeWithIdeogram(species = "hg38")
circos.genomicLink(as.data.frame(S4Vectors::first(pairs_ins)), as.data.frame(S4Vectors::second(pairs_ins)), col = "black")

#dup
circos.initializeWithIdeogram(species = "hg38")
circos.genomicLink(as.data.frame(S4Vectors::first(pairs_dup)), as.data.frame(S4Vectors::second(pairs_dup)), col = "darkorange")

#bnd
circos.initializeWithIdeogram(species = "hg38")
circos.genomicLink(as.data.frame(S4Vectors::first(pairs_bnd)), as.data.frame(S4Vectors::second(pairs_bnd)), col = "forestgreen")


dev.off()
}
