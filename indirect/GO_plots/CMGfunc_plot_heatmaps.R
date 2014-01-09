#!/bin/Rscript --slave
#------------------------------------------------- GET OPTIONS
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

#library(ggplot2)
library(plyr)
library(reshape)
library(pheatmap)
library(gplots)
library(grid)     ## Need to attach (and not just load) grid package
library(RColorBrewer)

#-------FUNCTION: annotate ---------
annotate <- function(freq) {
	annotationInfo <- colnames(freq)
	annotationInfo <- as.data.frame(annotationInfo)
	rownames(annotationInfo) <- colnames(freq)
	annotationInfo$annotationInfo <- factor(annotationInfo$annotationInfo, levels=unique(colnames(freq)))

	generateColours <- brewer.pal(n = 7, name = "RdYlBu")

	annotationColours <- generateColours[1:length(unique(colnames(freq)))]
	names(annotationColours) <- unique(colnames(freq))
	annotationColours <- list(annotationInfo=annotationColours)
	output <- list(info = annotationInfo, colors = annotationColours) 
	return (output)
}

#-------FUNCTION: reformatFreq ---------
reformatFreq <- function(freq_Long, column) {
	freq <- reshape(freq_Long, idvar=column,timevar='V1',v.names='V6', direction="wide")
	rownames(freq) <- freq[,1]
	freq <- as.matrix(freq)
	freq <- freq[,-1, drop=FALSE]
	mode(freq) <- 'numeric'
	freq[is.na(freq)] <- 0
	#sd <- apply(freq,1,sd)
	#freq <- freq[which(sd>=0.1),]
	return (freq)
}

#-------FUNCTION: reformatPer ---------
reformatPer <- function(per_Long, column) {
	per <- reshape(per_Long, idvar=column,timevar='V1',v.names='V7', direction="wide")
	rownames(per) <- per[,1]
	per <- as.matrix(per)
	per <- per[,-1, drop=FALSE]
	mode(per) <- 'numeric'
	per[is.na(per)] <- 0
	#sd <- apply(per,1,sd)
	#per <- per[which(sd>=0.1),]
	return (per)
}
#-------FUNCTION: reformatGO ---------
reformatGO <- function(go_Long) {
	go <- reshape(go_Long, idvar='V3',timevar='V1',v.names='V6', direction="wide")
	rownames(go) <- go[,1]
	go <- as.matrix(go)
	go <- go[,-1, drop=FALSE]
	mode(go) <- 'numeric'
	go[is.na(go)] <- 0
	#go <- apply(freq,1,sd)
	#go <- go[which(sd>=0.1),]
	return (go)
}

#############################################################################################
# START MAIN
#############################################################################################
#ffa500 (orange), #9e002b (red), #0000ff (purple), #0000ff (blue)

df <- read.table(args[1] , sep="\t")
freqcut <- as.numeric(args[2]) 
percut <- as.numeric(args[3]) 
widthsingle <- as.numeric(args[4]) 
widthmulti <- as.numeric(args[5]) 
sdcut <- as.numeric(args[6]) 

df$V6 <- as.integer(df$V6)
df$V7 <- as.numeric(as.character(df$V7))

# Frequencies
freq_mf <- subset(df, df$V6 >= freqcut & !grepl("GO", df$V2) & grepl("MF", df$V8), select=c(V1,V6,V8))
freq_cc <- subset(df, df$V6 >= freqcut & !grepl("GO", df$V2) & grepl("CC", df$V9), select=c(V1,V6,V9))
freq_bp <- subset(df, df$V6 >= freqcut & !grepl("GO", df$V2) & grepl("BP", df$V10), select=c(V1,V6,V10))
freq_mf_Long <- freq_mf[,1:3]
freq_bp_Long <- freq_bp[,1:3]
freq_cc_Long <- freq_cc[,1:3]

# Percentages
per_mf <- subset(df, df$V7 >= percut & !grepl("GO", df$V2) & grepl("MF", df$V8), select=c(V1,V7,V8))
per_cc <- subset(df, df$V7 >= percut & !grepl("GO", df$V2) & grepl("CC", df$V9), select=c(V1,V7,V9))
per_bp <- subset(df, df$V7 >= percut & !grepl("GO", df$V2) & grepl("BP", df$V10), select=c(V1,V7,V10))

per_mf_Long <- per_mf[,c(1:3)]
per_bp_Long <- per_bp[,c(1:3)]
per_cc_Long <- per_cc[,c(1:3)]

# GO term frequencies
go_mf <- subset(df, df$V6 >= freqcut & grepl("GO", df$V2) & df$V3!="NA" & df$V8=="molecular_function", select=c(V1,V6,V3))
go_cc <- subset(df, df$V6 >= freqcut & grepl("GO", df$V2) & df$V3!="NA" & df$V8=="cellular_component", select=c(V1,V6,V3))
go_bp <- subset(df, df$V6 >= freqcut & grepl("GO", df$V2) & df$V3!="NA" & df$V8=="biological_process", select=c(V1,V6,V3))
gomf_Long <- go_mf[,c(1:3)]
gocc_Long <- go_cc[,c(1:3)]
gobp_Long <- go_bp[,c(1:3)]

lapply((paste(c("GO_mf", length(unique( gomf_Long$V3))), collapse="::")), write, "_stats", append=F)
lapply((paste(c("GO_cc", length(unique( gocc_Long$V3))), collapse="::")), write, "_stats", append=T)
lapply((paste(c("GO_bp", length(unique( gobp_Long$V3))), collapse="::")), write, "_stats", append=T)

lapply((paste(c("Freq_mf", length(unique( freq_mf_Long$V8))), collapse="::")), write, "_stats", append=T)
lapply((paste(c("Freq_cc", length(unique( freq_cc_Long$V9))), collapse="::")), write, "_stats", append=T)
lapply((paste(c("Freq_bp", length(unique( freq_bp_Long$V10))), collapse="::")), write, "_stats", append=T)

lapply((paste(c("Per_mf", length(unique( per_mf_Long$V8))), collapse="::")), write, "_stats", append=T)
lapply((paste(c("Per_cc", length(unique( per_cc_Long$V9))), collapse="::")), write, "_stats", append=T)
lapply((paste(c("Per_bp", length(unique( per_bp_Long$V10))), collapse="::")), write, "_stats", append=T)

#############################################################################################
# SINGLE GENOME
#############################################################################################
if (length(levels(freq_cc_Long$V1)) == 1 & length(levels(freq_bp_Long$V1)) == 1 & length(levels(freq_cc_Long$V1)) == 1) {
	
	#############################################################################################
	# MF freq - SINGLE GENOME
	#############################################################################################
	if (dim(freq_mf_Long)[1] > 1 ) {
	
		# Reformat
		freq <- reformatFreq(freq_mf_Long, "V8")
		
		# Create annotations
		annotationInfo <- annotate(freq)[[1]]
		annotationColours <- annotate(freq)[[2]]

		# Set color range
		color_range <- length(seq(min(freq_mf_Long$V6), max(freq_mf_Long$V6), 5))

		# If all frequencies are the same, it is not possible to create a heatmap
		if (max(freq_mf_Long$V6) == min(freq_mf_Long$V6)) {
			pdf("plot_freqsingle_mf.pdf", paper="a4r",width = 0, height = 0)
			barplot(freq, col=c("#00ffa5"), main="Frequency of functions with MF GO, one genome", ylab="Frequency")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(freq, cluster_rows=T, cluster_cols=F, main = "Frequency of functions with MF GO, one genome",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white", "#ffa500","#00ffa5"))(color_range), 
			cellheight=5,filename ="plot_freqsingle_mf.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}
	} 
	#############################################################################################
	# BP freq - SINGLE GENOME
	#############################################################################################
	if (dim(freq_bp_Long)[1] > 1) {

		# Reformat
		freq <- reformatFreq(freq_bp_Long, "V10")

		# Create annotations
		annotationInfo <- annotate(freq)[[1]]
		annotationColours <- annotate(freq)[[2]]

		# Set color range
		color_range <- length(seq(min(freq_bp_Long$V6), max(freq_bp_Long$V6), 1))

		# If all frequencies are the same, it is not possible to create a heatmap
		if (max(freq_bp_Long$V6) == min(freq_bp_Long$V6)) {
			pdf("plot_freqsingle_bp.pdf", paper="a4r",width = 0, height = 0)
			barplot(freq, col=c("#ffa500"), main="Frequency of functions with BP GO, one genome", ylab="Frequency")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(freq, cluster_rows=T, cluster_cols=F, main = "Frequency of functions with BP GO, one genome",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white", "#ffa500","#005aff"))(color_range), 
			cellheight=5,filename="plot_freqsingle_bp.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}
	}
	#############################################################################################
	# CC freq - SINGLE GENOME
	#############################################################################################
	if (dim(freq_cc_Long)[1] > 1) {

		# Reformat
		freq <- reformatFreq(freq_cc_Long, "V9")

		# Create annotations
		annotationInfo <- annotate(freq)[[1]]
		annotationColours <- annotate(freq)[[2]]

		# Set color range
		color_range <- length(seq(min(freq_cc_Long$V6), max(freq_cc_Long$V6), 1))

		# If all frequencies are the same, it is not possible to create a heatmap
		if (max(freq_bp_Long$V6) == min(freq_bp_Long$V6)) {
			pdf("plot_freqsingle_cc.pdf", paper="a4r",width = 0, height = 0)
			barplot(freq, col=c("#ffa500"), main="Frequency of functions with CC GO, one genome", ylab="Frequency")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(freq, cluster_rows=T, cluster_cols=F, main = "Frequency of functions with CC GO, one genome",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white", "#ffa500","#a500ff"))(color_range), 
			cellheight=5,filename="plot_freqsingle_cc.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}
	}
}	

#############################################################################################
# Percentage - SINGLE GENOME
#############################################################################################
if (length(levels(per_mf_Long$V1)) == 1 & length(levels(per_bp_Long$V1)) == 1  & length(levels(per_cc_Long$V1)) == 1) {
	#############################################################################################
	# MF percentage - SINGLE GENOME
	#############################################################################################
	if (dim(per_mf_Long)[1] > 1 ) {
	
		# Reformat
		per <- reformatPer(per_mf_Long, "V8")

		# Create annotations
		annotationInfo <- annotate(per)[[1]]
		annotationColours <- annotate(per)[[2]]

		# Set color range
		color_range <- length(seq(min(per_mf_Long$V7), max(per_mf_Long$V7), 0.5))

		# If all peruencies are the same, it is not possible to create a heatmap
		if (max(per_mf_Long$V7) == min(per_mf_Long$V7)) {
			pdf("plot_persingle_mf.pdf", paper="a4r",width = 0, height = 0)
			barplot(per, col=c("#00ffa5"), main="Percentages of functions with MF GO, one genome", ylab="peruency")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(per, cluster_rows=T, cluster_cols=F, main = "Percentages of functions with MF GO, one genome",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white", "#ff005a","#00ffa5"))(color_range), 
			cellheight=5,filename ="plot_persingle_mf.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}
	} 
	#############################################################################################
	# BP percentage - SINGLE GENOME
	#############################################################################################
	if (dim(per_bp_Long)[1] > 1) {

		# Reformat
		per <- reformatPer(per_bp_Long, "V10")

		# Create annotations
		annotationInfo <- annotate(per)[[1]]
		annotationColours <- annotate(per)[[2]]

		# Set color range
		color_range <- length(seq(min(per_bp_Long$V7), max(per_bp_Long$V7), 0.5))

		# If all percentage are the same, it is not possible to create a heatmap
		if (max(per_bp_Long$V7) == min(per_bp_Long$V7)) {
			pdf("plot_persingle_bp.pdf", paper="a4r",width = 0, height = 0)
			barplot(per, col=c("#ae002b"), main="Percentages of functions with BP GO, one genome", ylab="percentage")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(per, cluster_rows=T, cluster_cols=F, main = "Percentages of functions with BP GO, one genome",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white", "#ff005a","#005aff"))(color_range), 
			cellheight=5,filename="plot_persingle_bp.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}
	}
	#############################################################################################
	# CC percentage - SINGLE GENOME
	#############################################################################################
	if (dim(per_cc_Long)[1] > 1) {

		# Reformat
		per <- reformatPer(per_cc_Long, "V9")

		# Create annotations
		annotationInfo <- annotate(per)[[1]]
		annotationColours <- annotate(per)[[2]]

		# Set color range
		color_range <- length(seq(min(per_cc_Long$V7), max(per_cc_Long$V7), 0.5))

		# If all percentage are the same, it is not possible to create a heatmap
		if (max(per_cc_Long$V7) == min(per_cc_Long$V7)) {
			pdf("plot_persingle_cc.pdf", paper="a4r",width = 0, height = 0)
			barplot(per, col=c("#ae002b"), main="Percentages of functions with CC GO, one genome", ylab="percentage")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(per, cluster_rows=T, cluster_cols=F, main = "Percentages of functions with CC GO, one genome",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white", "#ff005a","#a500ff"))(color_range), 
			cellheight=5,filename="plot_persingle_cc.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}
	}
}

#############################################################################################
# Go terms - single genome
#############################################################################################
if (length(levels(gomf_Long$V1)) == 1 & length(levels(gocc_Long$V1)) == 1 & length(levels(gobp_Long$V1)) == 1) {
	
	#############################################################################################
	# Molecular function frequency GO- SINGLE GENOME
	#############################################################################################
	if (dim(gomf_Long)[1] > 1 ) {

		# Reformat
		gomf <- reformatGO(gomf_Long)

		# Create annotations
		annotationInfo <- annotate(gomf)[[1]]
		annotationColours <- annotate(gomf)[[2]]

		# Set color range
		color_range <- length(seq(min(gomf_Long$V6), max(gomf_Long$V6), 1))

		# If all gomfuencies are the same, it is not possible to create a heatmap
		if (max(gomf_Long$V6) == min(gomf_Long$V6)) {
			pdf("plot_gomfsingle_mf.pdf", paper="a4r",width = 0, height = 0)
			barplot(gomf, col=c("#00a0ea"), main="Frequency of Molecular Function GO terms, one genome", ylab="frequency gomf")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(gomf, cluster_rows=T, cluster_cols=F, main = "Frequency of Molecular Function GO terms, one genome",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white", "#ffe500","#00ffa5"))(color_range), 
			cellheight=5,filename ="plot_gomfsingle_mf.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}

	} 
	#############################################################################################
	# Cellular component frequency GO- SINGLE GENOME
	#############################################################################################
	if (dim(gocc_Long)[1] > 1) {

		# Reformat
		gocc <- reformatGO(gocc_Long)

		# Create annotations
		annotationInfo <- annotate(gocc)[[1]]
		annotationColours <- annotate(gocc)[[2]]

		# Set color range
		color_range <- length(seq(min(gocc_Long$V6), max(gocc_Long$V6), 1))

		# If all goccuencies are the same, it is not possible to create a heatmap
		if (max(gocc_Long$V6) == min(gocc_Long$V6)) {
			pdf("plot_goccsingle_cc.pdf", paper="a4r",width = 0, height = 0)
			barplot(gocc, col=c("#00a0ea"), main="Frequency of Cellular Component GO terms, one genome", ylab="frequency gocc")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(gocc, cluster_rows=T, cluster_cols=F, main = "Frequency of Cellular Component GO terms, one genome",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white", "#ffe500","#a500ff"))(color_range), 
			cellheight=5,filename="plot_goccsingle_cc.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}
	}
	#############################################################################################
	# Biological Process frequency GO- SINGLE GENOME
	#############################################################################################
	if (dim(gobp_Long)[1] > 1) {

		# Reformat
		gobp <- reformatGO(gobp_Long)

		# Create annotations
		annotationInfo <- annotate(gobp)[[1]]
		annotationColours <- annotate(gobp)[[2]]

		# Set color range
		color_range <- length(seq(min(gobp_Long$V6), max(gobp_Long$V6), 1))

		# If all frequencies are the same, it is not possible to create a heatmap
		if (max(gobp_Long$V6) == min(gobp_Long$V6)) {
			pdf("plot_gobpsingle_bp.pdf", paper="a4r",width = 0, height = 0)
			barplot(gobp, col=c("#00a0ea"), main="Frequency of Biological Process GO terms, one genome", ylab="frequency gobp")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(gobp, cluster_rows=T, cluster_cols=F, main = "Frequency of Biological Process GO terms, one genome",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white", "#ffe500","#005aff"))(color_range), 
			cellheight=5,filename="plot_gobpsingle_bp.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}
	}
}	


#############################################################################################
# MULTIPLE GENOMES
#############################################################################################
if (length(levels(freq_cc_Long$V1)) > 1 & length(levels(freq_bp_Long$V1)) > 1 & length(levels(freq_cc_Long$V1)) > 1) {
	
	#############################################################################################
	# MF freq - MULTIPLE GENOMES
	#############################################################################################
	if (dim(freq_mf_Long)[1] > 1 ) {
	
		# Reformat
		freq <- reformatFreq(freq_mf_Long, "V8")
		sd <- apply(freq,1,sd)
		lapply((paste(c("Freq_mf_sdlow", length(sd[sd<sdcut])), collapse="::")), write, "_stats", append=T)
		freq <- freq[which(sd>=0.1),]

		# Create annotations
		annotationInfo <- annotate(freq)[[1]]
		annotationColours <- annotate(freq)[[2]]

		# Set color range
		color_range <- length(seq(min(freq_mf_Long$V6), max(freq_mf_Long$V6), 5))

		# If all frequencies are the same, it is not possible to create a heatmap
		if (max(freq_mf_Long$V6) == min(freq_mf_Long$V6)) {
			pdf("plot_freqmultiple_mf.pdf", paper="a4r",width = 0, height = 0)
			barplot(freq, col=c("#00ffa5"), main="Frequency of functions with MF GO, multiple genomes", ylab="Frequency")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(freq, cluster_rows=T, cluster_cols=T, main = "Frequency of functions with MF GO, multiple genomes",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white", "#ffa500","#00ffa5"))(color_range), 
			cellheight=5,filename ="plot_freqmultiple_mf.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}
	} 
	#############################################################################################
	# BP freq - MULTIPLE GENOMES
	#############################################################################################
	if (dim(freq_bp_Long)[1] > 1) {

		# Reformat
		freq <- reformatFreq(freq_bp_Long, "V10")
		sd <- apply(freq,1,sd)
		lapply((paste(c("Freq_bp_sdlow", length(sd[sd<sdcut])), collapse="::")), write, "_stats", append=T)
		freq <- freq[which(sd>=0.1),]

		# Create annotations
		annotationInfo <- annotate(freq)[[1]]
		annotationColours <- annotate(freq)[[2]]

		# Set color range
		color_range <- length(seq(min(freq_bp_Long$V6), max(freq_bp_Long$V6), 1))

		# If all frequencies are the same, it is not possible to create a heatmap
		if (max(freq_bp_Long$V6) == min(freq_bp_Long$V6)) {
			pdf("plot_freqmultiple_bp.pdf", paper="a4r",width = 0, height = 0)
			barplot(freq, col=c("#ffa500"), main="Frequency of functions with BP GO, multiple genomes", ylab="Frequency")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(freq, cluster_rows=T, cluster_cols=T, main = "Frequency of functions with BP GO, multiple genomes",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white", "#ffa500","#005aff"))(color_range), 
			cellheight=5,filename="plot_freqmultiple_bp.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}
	}
	#############################################################################################
	# CC freq - MULTIPLE GENOMES
	#############################################################################################
	if (dim(freq_cc_Long)[1] > 1) {

		# Reformat
		freq <- reformatFreq(freq_cc_Long, "V9")
		sd <- apply(freq,1,sd)
		lapply((paste(c("Freq_cc_sdlow", length(sd[sd<sdcut])), collapse="::")), write, "_stats", append=T)
		freq <- freq[which(sd>=0.1),]

		# Create annotations
		annotationInfo <- annotate(freq)[[1]]
		annotationColours <- annotate(freq)[[2]]

		# Set color range
		color_range <- length(seq(min(freq_cc_Long$V6), max(freq_cc_Long$V6), 1))

		# If all frequencies are the same, it is not possible to create a heatmap
		if (max(freq_bp_Long$V6) == min(freq_bp_Long$V6)) {
			pdf("plot_freqmultiple_cc.pdf", paper="a4r",width = 0, height = 0)
			barplot(freq, col=c("#ffa500"), main="Frequency of functions with CC GO, multiple genomes", ylab="Frequency")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(freq, cluster_rows=T, cluster_cols=T, main = "Frequency of functions with CC GO, multiple genomes",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white", "#ffa500","#a500ff"))(color_range), 
			cellheight=5,filename="plot_freqmultiple_cc.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}
	}
}	

#############################################################################################
# Percentage - MULTIPLE GENOMES
#############################################################################################
if (length(levels(per_mf_Long$V1)) > 1& length(levels(per_bp_Long$V1)) > 1 & length(levels(per_cc_Long$V1)) > 1) {
	#############################################################################################
	# MF percentage - MULTIPLE GENOMES
	#############################################################################################
	if (dim(per_mf_Long)[1] > 1 ) {
	
		# Reformat
		per <- reformatPer(per_mf_Long, "V8")
		sd <- apply(per,1,sd)
		lapply((paste(c("Per_mf_sdlow", length(sd[sd<sdcut])), collapse="::")), write, "_stats", append=T)
		per <- per[which(sd>=0.1),]

		# Create annotations
		annotationInfo <- annotate(per)[[1]]
		annotationColours <- annotate(per)[[2]]

		# Set color range
		color_range <- length(seq(min(per_mf_Long$V7), max(per_mf_Long$V7), 0.5))

		# If all peruencies are the same, it is not possible to create a heatmap
		if (max(per_mf_Long$V7) == min(per_mf_Long$V7)) {
			pdf("plot_permultiple_mf.pdf", paper="a4r",width = 0, height = 0)
			barplot(per, col=c("#00ffa5"), main="Percentages of functions with MF GO, multiple genomes", ylab="peruency")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(per, cluster_rows=T, cluster_cols=T, main = "Percentages of functions with MF GO, multiple genomes",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white", "#ff005a","#00ffa5"))(color_range), 
			cellheight=5,filename ="plot_permultiple_mf.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}
	} 
	#############################################################################################
	# BP percentage - MULTIPLE GENOMES
	#############################################################################################
	if (dim(per_bp_Long)[1] > 1) {

		# Reformat
		per <- reformatPer(per_bp_Long, "V10")
		sd <- apply(per,1,sd)
		lapply((paste(c("Per_bp_sdlow", length(sd[sd<sdcut])), collapse="::")), write, "_stats", append=T)
		per <- per[which(sd>=0.1),]

		# Create annotations
		annotationInfo <- annotate(per)[[1]]
		annotationColours <- annotate(per)[[2]]

		# Set color range
		color_range <- length(seq(min(per_bp_Long$V7), max(per_bp_Long$V7), 0.5))

		# If all percentage are the same, it is not possible to create a heatmap
		if (max(per_bp_Long$V7) == min(per_bp_Long$V7)) {
			pdf("plot_permultiple_bp.pdf", paper="a4r",width = 0, height = 0)
			barplot(per, col=c("#ae002b"), main="Percentages of functions with BP GO, multiple genomes", ylab="percentage")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(per, cluster_rows=T, cluster_cols=T, main = "Percentages of functions with BP GO, multiple genomes",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white", "#ff005a","#005aff"))(color_range), 
			cellheight=5,filename="plot_permultiple_bp.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}
	}
	#############################################################################################
	# CC percentage - MULTIPLE GENOMES
	#############################################################################################
	if (dim(per_cc_Long)[1] > 1) {

		# Reformat
		per <- reformatPer(per_cc_Long, "V9")
		sd <- apply(per,1,sd)
		lapply((paste(c("Per_cc_sdlow", length(sd[sd<sdcut])), collapse="::")), write, "_stats", append=T)
		per <- per[which(sd>=0.1),]


		# Create annotations
		annotationInfo <- annotate(per)[[1]]
		annotationColours <- annotate(per)[[2]]

		# Set color range
		color_range <- length(seq(min(per_cc_Long$V7), max(per_cc_Long$V7), 0.5))

		# If all percentage are the same, it is not possible to create a heatmap
		if (max(per_cc_Long$V7) == min(per_cc_Long$V7)) {
			pdf("plot_permultiple_cc.pdf", paper="a4r",width = 0, height = 0)
			barplot(per, col=c("#ae002b"), main="Percentages of functions with CC GO, multiple genomes", ylab="percentage")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(per, cluster_rows=T, cluster_cols=T, main = "Percentages of functions with CC GO, multiple genomes",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white", "#ff005a","#a500ff"))(color_range), 
			cellheight=5,filename="plot_permultiple_cc.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}
	}
}
	

#############################################################################################
# Go terms - MULTIPLE GENOMES
#############################################################################################
if (length(levels(gomf_Long$V1)) > 1 & length(levels(gocc_Long$V1)) > 1 & length(levels(gobp_Long$V1)) > 1) {
	
	#############################################################################################
	# Molecular function frequency GO- MULTIPLE GENOMES
	#############################################################################################
	if (dim(gomf_Long)[1] > 1 ) {

		# Reformat
		gomf <- reformatGO(gomf_Long)
		sd <- apply(gomf,1,sd)
		lapply((paste(c("GO_mf_sdlow", length(sd[sd<sdcut])), collapse="::")), write, "_stats", append=T)
		gomf <- gomf[which(sd>=0.1),]

		# Create annotations
		annotationInfo <- annotate(gomf)[[1]]
		annotationColours <- annotate(gomf)[[2]]

		# Set color range
		color_range <- length(seq(min(gomf_Long$V6), max(gomf_Long$V6), 1))

		# If all gomfuencies are the same, it is not possible to create a heatmap
		if (max(gomf_Long$V6) == min(gomf_Long$V6)) {
			pdf("plot_gomfmultiple_mf.pdf", paper="a4r",width = 0, height = 0)
			barplot(gomf, col=c("#00a0ea"), main="Frequency of Molecular Function GO terms, multiple genomes", ylab="frequency gomf")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(gomf, cluster_rows=T, cluster_cols=T, main = "Frequency of Molecular Function GO terms, multiple genomes",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white", "#ffe500","#00ffa5"))(color_range), 
			cellheight=5,filename ="plot_gomfmultiple_mf.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}

	} 
	#############################################################################################
	# Cellular component frequency GO- MULTIPLE GENOMES
	#############################################################################################
	if (dim(gocc_Long)[1] > 1) {

		# Reformat
		gocc <- reformatGO(gocc_Long)
		sd <- apply(gocc,1,sd)
		lapply((paste(c("GO_cc_sdlow", length(sd[sd<sdcut])), collapse="::")), write, "_stats", append=T)
		gocc <- gocc[which(sd>=0.1),]

		# Create annotations
		annotationInfo <- annotate(gocc)[[1]]
		annotationColours <- annotate(gocc)[[2]]

		# Set color range
		color_range <- length(seq(min(gocc_Long$V6), max(gocc_Long$V6), 1))

		# If all goccuencies are the same, it is not possible to create a heatmap
		if (max(gocc_Long$V6) == min(gocc_Long$V6)) {
			pdf("plot_goccmultiple_cc.pdf", paper="a4r",width = 0, height = 0)
			barplot(gocc, col=c("#00a0ea"), main="Frequency of Cellular Component GO terms, multiple genomes", ylab="frequency gocc")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(gocc, cluster_rows=T, cluster_cols=T, main = "Frequency of Cellular Component GO terms, multiple genomes",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white", "#ffe500","#a500ff"))(color_range), 
			cellheight=5,filename="plot_goccmultiple_cc.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}
	}
	#############################################################################################
	# Biological Process frequency GO- MULTIPLE GENOMES
	#############################################################################################
	if (dim(gobp_Long)[1] > 1) {

		# Reformat
		gobp <- reformatGO(gobp_Long)
		sd <- apply(gobp,1,sd)
		lapply((paste(c("GO_bp_sdlow", length(sd[sd<sdcut])), collapse="::")), write, "_stats", append=T)
		gobp <- gobp[which(sd>=0.1),]

		# Create annotations
		annotationInfo <- annotate(gobp)[[1]]
		annotationColours <- annotate(gobp)[[2]]

		# Set color range
		color_range <- length(seq(min(gobp_Long$V6), max(gobp_Long$V6), 1))

		# If all frequencies are the same, it is not possible to create a heatmap
		if (max(gobp_Long$V6) == min(gobp_Long$V6)) {
			pdf("plot_gobpmultiple_bp.pdf", paper="a4r",width = 0, height = 0)
			barplot(gobp, col=c("#00a0ea"), main="Frequency of Biological Process GO terms, multiple genomes", ylab="frequency gobp")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(gobp, cluster_rows=T, cluster_cols=T, main = "Frequency of Biological Process GO terms, multiple genomes",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white", "#ffe500","#005aff"))(color_range), 
			cellheight=5,filename="plot_gobpmultiple_bp.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}
	}

}	





	
