#ONLY required arg is the region number
args <- commandArgs(TRUE)
myregion <- as.integer(args[1])

Merger_mgl_output <- function(region=1){

	for (i in region){
		haplotypes <- list.files(full.names=F)[grep(pattern=paste0("\\.",i,".tab.mlg"),list.files(,full.names=F))]
		n_haplo <- length(haplotypes)

		for(j in 1:length(haplotypes)){
			assign(paste(haplotypes[j]), read.delim(toString(haplotypes[j]),header=T))
		}
		
		merge.1 <- merge(eval(as.name(haplotypes[1])),eval(as.name(haplotypes[2])),by="single_motif")
		
		if(n_haplo>2){
			for(k in 3:length(haplotypes)){
				merge.1 <-  merge(merge.1,eval(as.name(haplotypes[k])),by="single_motif")
			}
		}

		clean_merge <- merge.1[,c(1,2,seq(from=4,to=(n_haplo*3+1),by=3))]
		
		colnames(clean_merge) <- c("motif_seq","position",haplotypes)
		out_name <- paste0(region,".",n_haplo,"_haplotypes",".tab.mlg")
		write.table(clean_merge,out_name,sep="\t",quote=F,row.names=F)
	}
}

#Run function above using the Rscript arg[1]
Merger_mgl_output(region=myregion)

