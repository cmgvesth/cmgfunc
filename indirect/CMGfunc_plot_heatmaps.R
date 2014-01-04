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

#------------------------------------------------- FUNTION END

#############################################################################################
# START MAIN
#############################################################################################
#ffa500 (orange), #9e002b (red), #0000ff (purple), #0000ff (blue)

df <- read.table(args[1] , sep="\t")
freqcut <- as.numeric(args[2]) 
percut <- as.numeric(args[3]) 
widthsingle <- as.numeric(args[4]) 
widthmulti <- as.numeric(args[5]) 

#df <- read.table("CMGfunc/_tmpfunc" , sep="\t") # args[1]
#freqcut <- as.numeric(10) # args[2]

df$V6 <- as.integer(df$V6)
df$V7 <- as.numeric(as.character(df$V7))

# Frequencies
freq_high <- subset(df, df$V6 >= freqcut & !grepl("GO", df$V2), select=c(V1,V6,V8))
freq_low <- subset(df, df$V6 < freqcut & !grepl("GO", df$V2), select=c(V1,V6,V8))

freq_high_Long <- freq_high[,1:3]
freq_low_Long <- freq_low[,1:3]

# Percentages
per_high <- subset(df, df$V7 >= percut & !grepl("GO", df$V2), select=c(V1,V7,V8))
per_low  <- subset(df, df$V7 < percut & !grepl("GO", df$V2), select=c(V1,V7,V8))

per_high_Long <- per_high[,c(1:3)]
per_low_Long  <- per_low[,c(1:3)]

# GO term frequencies
go_mf <- subset(df, df$V6 >= freqcut & grepl("GO", df$V2) & df$V3!="NA" & df$V8=="molecular_function", select=c(V1,V6,V3))
go_cc <- subset(df, df$V6 >= freqcut & grepl("GO", df$V2) & df$V3!="NA" & df$V8=="cellular_component", select=c(V1,V6,V3))
go_bp <- subset(df, df$V6 >= freqcut & grepl("GO", df$V2) & df$V3!="NA" & df$V8=="biological_process", select=c(V1,V6,V3))

gomf_Long <- go_mf[,c(1:3)]
gocc_Long <- go_cc[,c(1:3)]
gobp_Long <- go_bp[,c(1:3)]

#############################################################################################
# SINGLE GENOME
#############################################################################################
if (length(levels(freq_high_Long$V1)) == 1 & length(levels(freq_low_Long$V1)) == 1) {
	
	#############################################################################################
	# High freq - SINGLE GENOME
	#############################################################################################
	if (dim(freq_high_Long)[1] > 1 ) {
	
		freq <- reshape(freq_high_Long, idvar='V8',timevar='V1',v.names='V6', direction="wide")
		rownames(freq) <- freq[,1]
		freq <- as.matrix(freq)
		freq <- freq[,-1, drop=FALSE]
		mode(freq) <- 'numeric'
		freq[is.na(freq)] <- 0

		# Create annotations
		annotationInfo <- colnames(freq)
		annotationInfo <- as.data.frame(annotationInfo)
		rownames(annotationInfo) <- colnames(freq)
		annotationInfo$annotationInfo <- factor(annotationInfo$annotationInfo, levels=unique(colnames(freq)))
		
		generateColours <- brewer.pal(n = 7, name = "RdYlBu")

		annotationColours <- generateColours[1:length(unique(colnames(freq)))]
		names(annotationColours) <- unique(colnames(freq))
		annotationColours <- list(annotationInfo=annotationColours)

		# Set color range
		if (dim(freq_low_Long)[1] > 1) { from <- max(freq_low_Long$V6) 
		} else from <-  min(freq_high_Long$V6) 
		color_range <- length(seq(from, max(freq_high_Long$V6), 5))


		# If all frequencies are the same, it is not possible to create a heatmap
		if (max(freq_high_Long$V6) == min(freq_high_Long$V6)) {
			string <- paste("plot_freqsingle_high.pdf")
			pdf(string, paper="a4r",width = 0, height = 0)
			barplot(freq, col=c("#00ffa5"), main="High frequency clans, one genome", ylab="Frequency")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(freq, cluster_rows=T, cluster_cols=F, main = "High frequency functions, one genome",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("#ffa500","#00ffa5"))(color_range), 
			cellheight=5,filename ="plot_freqsingle_high.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}
	} 
	#############################################################################################
	# Low freq - SINGLE GENOME
	#############################################################################################
	if (dim(freq_low_Long)[1] > 1) {

		freq <- reshape(freq_low_Long, idvar='V8',timevar='V1',v.names='V6', direction="wide")
		rownames(freq) <- freq[,1]
		freq <- as.matrix(freq)
		freq <- freq[,-1, drop=FALSE]
		mode(freq) <- 'numeric'
		freq[is.na(freq)] <- 0

		# Create annotations
		annotationInfo <- colnames(freq)
		annotationInfo <- as.data.frame(annotationInfo)
		rownames(annotationInfo) <- colnames(freq)
		annotationInfo$annotationInfo <- factor(annotationInfo$annotationInfo, levels=unique(colnames(freq)))
		
		generateColours <- brewer.pal(n = 7, name = "RdYlBu")

		annotationColours <- generateColours[1:length(unique(colnames(freq)))]
		names(annotationColours) <- unique(colnames(freq))
		annotationColours <- list(annotationInfo=annotationColours)

		# Set color range
		if (dim(freq_high_Long)[1] > 1) { to <- min(freq_high_Long$V6) 
		} else to <-  max(freq_low_Long$V6) 
		color_range <- length(seq(min(freq_low_Long$V6), to, 1))

		# If all frequencies are the same, it is not possible to create a heatmap
		if (max(freq_low_Long$V6) == min(freq_low_Long$V6)) {
			string <- paste("plot_freqsingle_low.pdf")
			pdf(string, paper="a4r",width = 0, height = 0)
			barplot(freq, col=c("#ffa500"), main="Low frequency clans, one genome", ylab="Frequency")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(freq, cluster_rows=T, cluster_cols=F, main = "Low frequency functions, one genome",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white","#ffa500"))(color_range), 
			cellheight=5,filename="plot_freqsingle_low.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}
	}
}	

#############################################################################################
# Percentage - SINGLE GENOME
#############################################################################################
if (length(levels(per_high_Long$V1)) == 1 & length(levels(per_low_Long$V1)) == 1) {
	#############################################################################################
	# High percentage - SINGLE GENOME
	#############################################################################################
	if (dim(per_high_Long)[1] > 1 ) {
	
		per <- reshape(per_high_Long, idvar='V8',timevar='V1',v.names='V7', direction="wide")
		rownames(per) <- per[,1]
		per <- as.matrix(per)
		per <- per[,-1, drop=FALSE]
		mode(per) <- 'numeric'
		per[is.na(per)] <- 0

		# Create annotations
		annotationInfo <- colnames(per)
		annotationInfo <- as.data.frame(annotationInfo)
		rownames(annotationInfo) <- colnames(per)
		annotationInfo$annotationInfo <- factor(annotationInfo$annotationInfo, levels=unique(colnames(per)))
		
		generateColours <- brewer.pal(n = 7, name = "RdYlBu")

		annotationColours <- generateColours[1:length(unique(colnames(per)))]
		names(annotationColours) <- unique(colnames(per))
		annotationColours <- list(annotationInfo=annotationColours)

		# Set color range
		if (dim(per_low_Long)[1] > 1) { from <- max(per_low_Long$V7) 
		} else from <-  min(per_high_Long$V7) 
		color_range <- length(seq(from, max(per_high_Long$V7), 0.5))

		# If all peruencies are the same, it is not possible to create a heatmap
		if (max(per_high_Long$V7) == min(per_high_Long$V7)) {
			string <- paste("plot_persingle_high.pdf")
			pdf(string, paper="a4r",width = 0, height = 0)
			barplot(per, col=c("#00ffa5"), main="High peruency clans, one genome", ylab="peruency")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(per, cluster_rows=T, cluster_cols=F, main = "High percentage functions, one genome",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("#ae002b","#00ffa5"))(color_range), 
			cellheight=5,filename ="plot_persingle_high.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}
	} 
	#############################################################################################
	# Low percentage - SINGLE GENOME
	#############################################################################################
	if (dim(per_low_Long)[1] > 1) {

		per <- reshape(per_low_Long, idvar='V8',timevar='V1',v.names='V7', direction="wide")
		rownames(per) <- per[,1]
		per <- as.matrix(per)
		per <- per[,-1, drop=FALSE]
		mode(per) <- 'numeric'
		per[is.na(per)] <- 0

		# Create annotations
		annotationInfo <- colnames(per)
		annotationInfo <- as.data.frame(annotationInfo)
		rownames(annotationInfo) <- colnames(per)
		annotationInfo$annotationInfo <- factor(annotationInfo$annotationInfo, levels=unique(colnames(per)))
		
		generateColours <- brewer.pal(n = 7, name = "RdYlBu")

		annotationColours <- generateColours[1:length(unique(colnames(per)))]
		names(annotationColours) <- unique(colnames(per))
		annotationColours <- list(annotationInfo=annotationColours)

		# Set color range
		if (dim(per_high_Long)[1] > 1) { to <- min(per_high_Long$V7) 
		} else to <-  max(per_low_Long$V7) 
		color_range <- length(seq(min(per_low_Long$V7), to, 0.01))

		# If all percentage are the same, it is not possible to create a heatmap
		if (max(per_low_Long$V7) == min(per_low_Long$V7)) {
			pdf("plot_persingle_low.pdf", paper="a4r",width = 0, height = 0)
			barplot(per, col=c("#ae002b"), main="Low percentage clans, one genome", ylab="percentage")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(per, cluster_rows=T, cluster_cols=F, main = "Low percentage functions, one genome",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white","#ae002b"))(color_range), 
			cellheight=5,filename="plot_persingle_low.pdf",width=widthsingle,treeheight_row= 150,
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
		gomf <- reshape(gomf_Long, idvar='V3',timevar='V1',v.names='V6', direction="wide")
		rownames(gomf) <- gomf[,1]
		gomf <- as.matrix(gomf)
		gomf <- gomf[,-1, drop=FALSE]
		mode(gomf) <- 'numeric'
		gomf[is.na(gomf)] <- 0

		# Create annotations
		annotationInfo <- colnames(gomf)
		annotationInfo <- as.data.frame(annotationInfo)
		rownames(annotationInfo) <- colnames(gomf)
		annotationInfo$annotationInfo <- factor(annotationInfo$annotationInfo, levels=unique(colnames(gomf)))
		generateColours <- brewer.pal(n = 7, name = "RdYlBu")
		annotationColours <- generateColours[1:length(unique(colnames(gomf)))]
		names(annotationColours) <- unique(colnames(gomf))
		annotationColours <- list(annotationInfo=annotationColours)

		# Set color range
		color_range <- length(seq(min(gomf_Long$V6), max(gomf_Long$V6), 1))

		# If all gomfuencies are the same, it is not possible to create a heatmap
		if (max(gomf_Long$V6) == min(gomf_Long$V6)) {
			pdf("plot_gomfsingle_mf.pdf", paper="a4r",width = 0, height = 0)
			barplot(gomf, col=c("#00a0ea"), main="Frequency of Molecular Function GO terms, one genome", ylab="frequency gomf")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(gomf, cluster_rows=T, cluster_cols=F, main = "Frequency of Molecular Function GO terms, one genome",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white","#ea002b", "#00a0ea"))(color_range), 
			cellheight=5,filename ="plot_gomfsingle_mf.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}

	} 
	#############################################################################################
	# Cellular component frequency GO- SINGLE GENOME
	#############################################################################################
	if (dim(gocc_Long)[1] > 1) {

		gocc <- reshape(gocc_Long, idvar='V3',timevar='V1',v.names='V6', direction="wide")
		rownames(gocc) <- gocc[,1]
		gocc <- as.matrix(gocc)
		gocc <- gocc[,-1, drop=FALSE]
		mode(gocc) <- 'numeric'
		gocc[is.na(gocc)] <- 0

		# Create annotations
		annotationInfo <- colnames(gocc)
		annotationInfo <- as.data.frame(annotationInfo)
		rownames(annotationInfo) <- colnames(gocc)
		annotationInfo$annotationInfo <- factor(annotationInfo$annotationInfo, levels=unique(colnames(gocc)))
		generateColours <- brewer.pal(n = 7, name = "RdYlBu")
		annotationColours <- generateColours[1:length(unique(colnames(gocc)))]
		names(annotationColours) <- unique(colnames(gocc))
		annotationColours <- list(annotationInfo=annotationColours)

		# Set color range
		color_range <- length(seq(min(gocc_Long$V6), max(gocc_Long$V6), 1))

		# If all goccuencies are the same, it is not possible to create a heatmap
		if (max(gocc_Long$V6) == min(gocc_Long$V6)) {
			pdf("plot_goccsingle_cc.pdf", paper="a4r",width = 0, height = 0)
			barplot(gocc, col=c("#00a0ea"), main="Frequency of Cellular Component GO terms, one genome", ylab="frequency gocc")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(gocc, cluster_rows=T, cluster_cols=F, main = "Frequency of Cellular Component GO terms, one genome",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white","#bf00ea", "#00a0ea"))(color_range), 
			cellheight=5,filename="plot_goccsingle_cc.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}
	}
	#############################################################################################
	# Biological Process frequency GO- SINGLE GENOME
	#############################################################################################
	if (dim(gobp_Long)[1] > 1) {

		gobp <- reshape(gobp_Long, idvar='V3',timevar='V1',v.names='V6', direction="wide")
		rownames(gobp) <- gobp[,1]
		gobp <- as.matrix(gobp)
		gobp <- gobp[,-1, drop=FALSE]
		mode(gobp) <- 'numeric'
		gobp[is.na(gobp)] <- 0

		# Create annotations
		annotationInfo <- colnames(gobp)
		annotationInfo <- as.data.frame(annotationInfo)
		rownames(annotationInfo) <- colnames(gobp)
		annotationInfo$annotationInfo <- factor(annotationInfo$annotationInfo, levels=unique(colnames(gobp)))
		generateColours <- brewer.pal(n = 7, name = "RdYlBu")
		annotationColours <- generateColours[1:length(unique(colnames(gobp)))]
		names(annotationColours) <- unique(colnames(gobp))
		annotationColours <- list(annotationInfo=annotationColours)

		# Set color range
		color_range <- length(seq(min(gobp_Long$V6), max(gobp_Long$V6), 1))

		# If all frequencies are the same, it is not possible to create a heatmap
		if (max(gobp_Long$V6) == min(gobp_Long$V6)) {
			pdf("plot_gobpsingle_bp.pdf", paper="a4r",width = 0, height = 0)
			barplot(gobp, col=c("#00a0ea"), main="Frequency of Biological Process GO terms, one genome", ylab="frequency gobp")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(gobp, cluster_rows=T, cluster_cols=F, main = "Frequency of Biological Process GO terms, one genome",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white","#ea4a00", "#00a0ea"))(color_range), 
			cellheight=5,filename="plot_gobpsingle_bp.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}
	}
}	



#############################################################################################
# MULTIPLE GENOMES
#############################################################################################
if (length(levels(freq_high_Long$V1)) > 1 & length(levels(freq_low_Long$V1)) > 1) {
	
	#############################################################################################
	# High freq - MULTIPLE GENOMES
	#############################################################################################
	if (dim(freq_high_Long)[1] > 1 ) {
	
		freq <- reshape(freq_high_Long, idvar='V8',timevar='V1',v.names='V6', direction="wide")
		rownames(freq) <- freq[,1]
		freq <- as.matrix(freq)
		freq <- freq[,-1, drop=FALSE]
		mode(freq) <- 'numeric'
		freq[is.na(freq)] <- 0

		# Create annotations
		annotationInfo <- colnames(freq)
		annotationInfo <- as.data.frame(annotationInfo)
		rownames(annotationInfo) <- colnames(freq)
		annotationInfo$annotationInfo <- factor(annotationInfo$annotationInfo, levels=unique(colnames(freq)))
		
		generateColours <- brewer.pal(n = 7, name = "RdYlBu")

		annotationColours <- generateColours[1:length(unique(colnames(freq)))]
		names(annotationColours) <- unique(colnames(freq))
		annotationColours <- list(annotationInfo=annotationColours)

		# Set color range
		if (dim(freq_low_Long)[1] > 1) { from <- max(freq_low_Long$V6) 
		} else from <-  min(freq_high_Long$V6) 
		color_range <- length(seq(from, max(freq_high_Long$V6), 5))


		# If all frequencies are the same, it is not possible to create a heatmap
		if (max(freq_high_Long$V6) == min(freq_high_Long$V6)) {
			string <- paste("plot_freqmulti_high.pdf")
			pdf(string, paper="a4r",width = 0, height = 0)
			barplot(freq, col=c("#00ffa5"), main="High frequency clans, multiple genomes", ylab="Frequency")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(freq, cluster_rows=T, cluster_cols=T, main = "High frequency functions, multiple genomes",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("#ffa500","#00ffa5"))(color_range), 
			cellheight=5,filename ="plot_freqmulti_high.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}
	} 
	#############################################################################################
	# Low freq - MULTIPLE GENOMES
	#############################################################################################
	if (dim(freq_low_Long)[1] > 1) {

		freq <- reshape(freq_low_Long, idvar='V8',timevar='V1',v.names='V6', direction="wide")
		rownames(freq) <- freq[,1]
		freq <- as.matrix(freq)
		freq <- freq[,-1, drop=FALSE]
		mode(freq) <- 'numeric'
		freq[is.na(freq)] <- 0

		# Create annotations
		annotationInfo <- colnames(freq)
		annotationInfo <- as.data.frame(annotationInfo)
		rownames(annotationInfo) <- colnames(freq)
		annotationInfo$annotationInfo <- factor(annotationInfo$annotationInfo, levels=unique(colnames(freq)))
		
		generateColours <- brewer.pal(n = 7, name = "RdYlBu")

		annotationColours <- generateColours[1:length(unique(colnames(freq)))]
		names(annotationColours) <- unique(colnames(freq))
		annotationColours <- list(annotationInfo=annotationColours)

		# Set color range
		if (dim(freq_high_Long)[1] > 1) { to <- min(freq_high_Long$V6) 
		} else to <-  max(freq_low_Long$V6) 
		color_range <- length(seq(min(freq_low_Long$V6), to, 1))

		# If all frequencies are the same, it is not possible to create a heatmap
		if (max(freq_low_Long$V6) == min(freq_low_Long$V6)) {
			string <- paste("plot_freqmulti_low.pdf")
			pdf(string, paper="a4r",width = 0, height = 0)
			barplot(freq, col=c("#ffa500"), main="Low frequency clans, multiple genomes", ylab="Frequency")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(freq, cluster_rows=T, cluster_cols=T, main = "Low frequency functions, multiple genomes",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white","#ffa500"))(color_range), 
			cellheight=5,filename="plot_freqmulti_low.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}
	}
}	

#############################################################################################
# Percentage - MULTIPLE GENOMES
#############################################################################################
if (length(levels(per_high_Long$V1)) > 1 & length(levels(per_low_Long$V1)) > 1) {

	#############################################################################################
	# High percentage - MULTIPLE GENOMES
	#############################################################################################
	if (dim(per_high_Long)[1] > 1 ) {
	
		per <- reshape(per_high_Long, idvar='V8',timevar='V1',v.names='V7', direction="wide")
		rownames(per) <- per[,1]
		per <- as.matrix(per)
		per <- per[,-1, drop=FALSE]
		mode(per) <- 'numeric'
		per[is.na(per)] <- 0

		# Create annotations
		annotationInfo <- colnames(per)
		annotationInfo <- as.data.frame(annotationInfo)
		rownames(annotationInfo) <- colnames(per)
		annotationInfo$annotationInfo <- factor(annotationInfo$annotationInfo, levels=unique(colnames(per)))
		
		generateColours <- brewer.pal(n = 7, name = "RdYlBu")

		annotationColours <- generateColours[1:length(unique(colnames(per)))]
		names(annotationColours) <- unique(colnames(per))
		annotationColours <- list(annotationInfo=annotationColours)

		# Set color range
		if (dim(per_low_Long)[1] > 1) { from <- max(per_low_Long$V7) 
		} else from <-  min(per_high_Long$V7) 
		color_range <- length(seq(from, max(per_high_Long$V7), 0.5))

		# If all peruencies are the same, it is not possible to create a heatmap
		if (max(per_high_Long$V7) == min(per_high_Long$V7)) {
			string <- paste("plot_permulti_high.pdf")
			pdf(string, paper="a4r",width = 0, height = 0)
			barplot(per, col=c("#00ffa5"), main="High peruency clans, multiple genomes", ylab="peruency")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(per, cluster_rows=T, cluster_cols=T, main = "High percentage functions, multiple genomes",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("#ae002b","#00ffa5"))(color_range), 
			cellheight=5,filename ="plot_permulti_high.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}
	} 
	#############################################################################################
	# Low percentage - MULTIPLE GENOMES
	#############################################################################################
	if (dim(per_low_Long)[1] > 1) {

		per <- reshape(per_low_Long, idvar='V8',timevar='V1',v.names='V7', direction="wide")
		rownames(per) <- per[,1]
		per <- as.matrix(per)
		per <- per[,-1, drop=FALSE]
		mode(per) <- 'numeric'
		per[is.na(per)] <- 0

		# Create annotations
		annotationInfo <- colnames(per)
		annotationInfo <- as.data.frame(annotationInfo)
		rownames(annotationInfo) <- colnames(per)
		annotationInfo$annotationInfo <- factor(annotationInfo$annotationInfo, levels=unique(colnames(per)))
		
		generateColours <- brewer.pal(n = 7, name = "RdYlBu")

		annotationColours <- generateColours[1:length(unique(colnames(per)))]
		names(annotationColours) <- unique(colnames(per))
		annotationColours <- list(annotationInfo=annotationColours)

		# Set color range
		if (dim(per_high_Long)[1] > 1) { to <- min(per_high_Long$V7) 
		} else to <-  max(per_low_Long$V7) 
		color_range <- length(seq(min(per_low_Long$V7), to, 0.01))

		# If all percentage are the same, it is not possible to create a heatmap
		if (max(per_low_Long$V7) == min(per_low_Long$V7)) {
			pdf("plot_permulti_low.pdf", paper="a4r",width = 0, height = 0)
			barplot(per, col=c("#ae002b"), main="Low percentage clans, multiple genomes", ylab="percentage")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(per, cluster_rows=T, cluster_cols=T, main = "Low percentage functions, multiple genomes",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white","#ae002b"))(color_range), 
			cellheight=5,filename="plot_permulti_low.pdf",width=widthsingle,treeheight_row= 150,
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
		gomf <- reshape(gomf_Long, idvar='V3',timevar='V1',v.names='V6', direction="wide")
		rownames(gomf) <- gomf[,1]
		gomf <- as.matrix(gomf)
		gomf <- gomf[,-1, drop=FALSE]
		mode(gomf) <- 'numeric'
		gomf[is.na(gomf)] <- 0

		# Create annotations
		annotationInfo <- colnames(gomf)
		annotationInfo <- as.data.frame(annotationInfo)
		rownames(annotationInfo) <- colnames(gomf)
		annotationInfo$annotationInfo <- factor(annotationInfo$annotationInfo, levels=unique(colnames(gomf)))
		generateColours <- brewer.pal(n = 7, name = "RdYlBu")
		annotationColours <- generateColours[1:length(unique(colnames(gomf)))]
		names(annotationColours) <- unique(colnames(gomf))
		annotationColours <- list(annotationInfo=annotationColours)

		# Set color range
		color_range <- length(seq(min(gomf_Long$V6), max(gomf_Long$V6), 1))

		# If all gomfuencies are the same, it is not possible to create a heatmap
		if (max(gomf_Long$V6) == min(gomf_Long$V6)) {
			pdf("plot_gomfmulti_mf.pdf", paper="a4r",width = 0, height = 0)
			barplot(gomf, col=c("#00a0ea"), main="Frequency of Molecular Function GO terms, multiple genomes", ylab="frequency gomf")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(gomf, cluster_rows=T, cluster_cols=T, main = "Frequency of Molecular Function GO terms, multiple genomes",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white","#ea002b", "#00a0ea"))(color_range), 
			cellheight=5,filename ="plot_gomfmulti_mf.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}

	} 
	#############################################################################################
	# Cellular component frequency GO- MULTIPLE GENOMES
	#############################################################################################
	if (dim(gocc_Long)[1] > 1) {

		gocc <- reshape(gocc_Long, idvar='V3',timevar='V1',v.names='V6', direction="wide")
		rownames(gocc) <- gocc[,1]
		gocc <- as.matrix(gocc)
		gocc <- gocc[,-1, drop=FALSE]
		mode(gocc) <- 'numeric'
		gocc[is.na(gocc)] <- 0

		# Create annotations
		annotationInfo <- colnames(gocc)
		annotationInfo <- as.data.frame(annotationInfo)
		rownames(annotationInfo) <- colnames(gocc)
		annotationInfo$annotationInfo <- factor(annotationInfo$annotationInfo, levels=unique(colnames(gocc)))
		generateColours <- brewer.pal(n = 7, name = "RdYlBu")
		annotationColours <- generateColours[1:length(unique(colnames(gocc)))]
		names(annotationColours) <- unique(colnames(gocc))
		annotationColours <- list(annotationInfo=annotationColours)

		# Set color range
		color_range <- length(seq(min(gocc_Long$V6), max(gocc_Long$V6), 1))

		# If all goccuencies are the same, it is not possible to create a heatmap
		if (max(gocc_Long$V6) == min(gocc_Long$V6)) {
			pdf("plot_goccmulti_cc.pdf", paper="a4r",width = 0, height = 0)
			barplot(gocc, col=c("#00a0ea"), main="Frequency of Cellular Component GO terms, multiple genomes", ylab="frequency gocc")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(gocc, cluster_rows=T, cluster_cols=T, main = "Frequency of Cellular Component GO terms, multiple genomes",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white","#bf00ea", "#00a0ea"))(color_range), 
			cellheight=5,filename="plot_goccmulti_cc.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}
	}
	#############################################################################################
	# Biological Process frequency GO- MULTIPLE GENOMES
	#############################################################################################
	if (dim(gobp_Long)[1] > 1) {

		gobp <- reshape(gobp_Long, idvar='V3',timevar='V1',v.names='V6', direction="wide")
		rownames(gobp) <- gobp[,1]
		gobp <- as.matrix(gobp)
		gobp <- gobp[,-1, drop=FALSE]
		mode(gobp) <- 'numeric'
		gobp[is.na(gobp)] <- 0

		# Create annotations
		annotationInfo <- colnames(gobp)
		annotationInfo <- as.data.frame(annotationInfo)
		rownames(annotationInfo) <- colnames(gobp)
		annotationInfo$annotationInfo <- factor(annotationInfo$annotationInfo, levels=unique(colnames(gobp)))
		generateColours <- brewer.pal(n = 7, name = "RdYlBu")
		annotationColours <- generateColours[1:length(unique(colnames(gobp)))]
		names(annotationColours) <- unique(colnames(gobp))
		annotationColours <- list(annotationInfo=annotationColours)

		# Set color range
		color_range <- length(seq(min(gobp_Long$V6), max(gobp_Long$V6), 1))

		# If all frequencies are the same, it is not possible to create a heatmap
		if (max(gobp_Long$V6) == min(gobp_Long$V6)) {
			pdf("plot_gobpmulti_bp.pdf", paper="a4r",width = 0, height = 0)
			barplot(gobp, col=c("#00a0ea"), main="Frequency of Biological Process GO terms, multiple genomes", ylab="frequency gobp")
		} else {
			# Plot heatmap with clustering on column
			pheatmap(gobp, cluster_rows=T, cluster_cols=T, main = "Frequency of Biological Process GO terms, multiple genomes",fontsize_row=5,fontsize=8,
			border_color="white",color = colorRampPalette(c("white","#ea4a00", "#00a0ea"))(color_range), 
			cellheight=5,filename="plot_gobpmulti_bp.pdf",width=widthsingle,treeheight_row= 150,
			annotation=annotationInfo, annotation_colours=annotationColours,show_colnames=F,show_rownames=T)
		}
	}
}	
