#Author: Arvis Sulovari, PhD
#Main functions required to detect and count repeat units of arbitrary length in FASTA records.   

require(seqinr)
require(Biostrings)
require(stringr)
require(GenomicRanges)


#################
####STEP 1/3#####
#################

##This function finds all pure repeat motifs for a single sequence
singleSeq_motif_finder <- function(seq_input=yor_hsa[[300]],coords=TRUE){
  motif_arr_with_coords <- array(NA,dim=c(1,3))
  motif_arr_no_coords <- array(NA,dim=c(1,1))
  myseq <- seq_input
  myseq_len <- length(myseq)
  i=1
  while(i<myseq_len) {
    first <- myseq[i]
    start_scnd_pointer <- (i+1)
    end_scnd_pointer <- (i+round((myseq_len-i)/2))
    for(j in start_scnd_pointer:end_scnd_pointer){
      second <- myseq[j]
      period=as.numeric(j-i-1)
      if(first == second && (toString(myseq[i:(i+period)]) == toString(myseq[j:(j+period)]))){
        single_motif <- toString(myseq[i:(i+period)])
        dup_start = i
        dup_end = (j+period)
        if(coords){
          motif_row_with_coords <- cbind(single_motif,dup_start,dup_end)
          motif_arr_with_coords <- rbind(motif_row_with_coords,motif_arr_with_coords)
        }
        else{
          motif_arr_no_coords <- rbind(single_motif,motif_arr_no_coords)
        }
      }
      else{
      }
    }
    i=i+1
  }
 if(coords){
   return(na.omit(unique(motif_arr_with_coords)))
 }
  else{
    colnames(motif_arr_no_coords) <- c("single_motif")
    row.names(motif_arr_no_coords) <- NULL
    return(na.omit(unique(motif_arr_no_coords)))
  }
}


#################
####STEP 2/3#####
#################

##GENOTYPER function
#For each region in the dataset above, search for the greatest number of tandemly repeated motifs.
PureRepeatCounter <- function(my_dataset=my_arr,output=out_arr,start=1,end=10,multicounts=T){
  #my_dataset <- unique(as.data.frame(cbind(motif_arr_out[,-c(3,4)],PureRepeatsCounter=rep(0,nrow(motif_arr_out)))))
  output <- as.data.frame(output)
  output$PureRepeatsCounter <- as.numeric(as.character(output$PureRepeatsCounter))
  levels(output$PureRepeatsCounter) <- c(1)
  for(i in start:end){
    #Iterate through the list of sequences, i.e., rows in the dataframe
    slice <- my_dataset[i,]
    myseq <- toString(slice$seq_of_slice)
    mymotif <- toString(slice$single_motif)
    table_tmp_out <- str_locate_all(myseq,mymotif)[[1]]
    mystart <- table_tmp_out[,1]
    myend <- table_tmp_out[,2]
    #Create IRanges object
    #coords <- GRanges(seqnames=rep(c("A"),length(mystart)),ranges=IRanges(start=mystart,end=myend))
    coords <- IRanges(start=mystart,end=myend)
    coords_red <- reduce(coords)
    mx_ir <- as.matrix(findOverlaps(coords,coords_red))
    cond = nrow(mx_ir)
    if(cond>0){
      myval <- as.numeric(table(mx_ir[,2]))
      #Count all occurrence of the motif OR the highest multiplicity (depending on multicounts parameter)
      if(multicounts){
        purerep_count <- as.array(myval)
        purerep_count <- purerep_count[which(purerep_count>1)]
        output[i,]$PureRepeatsCounter <- toString(purerep_count)
      }
      else{
        purerep_count <- max(myval)
        output[i,]$PureRepeatsCounter <- as.numeric(purerep_count)
      }
    }
    else {
      output[i,]$PureRepeatsCounter <- 0
    }
  }
  return(output[start:end,])
}



#################
####STEP 3/3#####
#################

##This function evokes the singleSeq_motif_finder() and PureRepeatCounter() to produce a slice of the final output
multiSeq_motif_finder <- function(fasta_input=yor_hsa){
  tmp_out <- array(NA,dim=c(1,4))
  colnames(tmp_out) <- c("coords_of_slice","single_motif","seq_of_slice","PureRepeatsCounter")
  n_seqs <- length(fasta_input)
  all_names <- names(fasta_input)
  for(i in 1:n_seqs){
    print(i)
    if(length(grep("-",toString(fasta_input[i])))==0){
      #Do not output the coords of each duplicated motif; just the unique motif for each FASTA record
      motifs_of_slice <- singleSeq_motif_finder(seq_input = fasta_input[[i]],coords = FALSE)
      ##Remove start-end of each motif
      motifs_of_slice[,1] <- toupper(gsub(", ","",motifs_of_slice[,1]))
      seq_of_slice <- toString(fasta_input[[i]])
      seq_of_slice <- toupper(gsub(", ","",seq_of_slice))
      coords_of_slice <- all_names[i]
      #Combine into input dataframe:
      bound_slice <- as.data.frame(unique(cbind(coords_of_slice,single_motif=motifs_of_slice,seq_of_slice,PureRepeatsCounter=rep(0,length(motifs_of_slice)))))
      #Count evoking PureRepeatCounter()
      bound_slice_counted <- PureRepeatCounter(my_dataset = bound_slice,output = bound_slice,start = 1,end=nrow(bound_slice),multicounts = F)
      ###
      tmp_out <- rbind(bound_slice_counted,tmp_out)
    }
    else{
    }
  }
  return(na.omit(tmp_out))
}


