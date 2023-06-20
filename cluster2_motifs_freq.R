install.packages('seqinr') # do this once
library(seqinr)
library(ggplot2)
library(gridExtra)

rm(list=ls()) #remove objects from specified envt

processFile = function(multifasta) {
  genes <- c()
  seqs <- c()
  con = file(multifasta, "r")
  sequence=""
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    print(line)
    ## header
    if(length(grep('>', line))){
      g=unlist(strsplit(line, split=" "))[1]
      g= gsub('>','',g) ### Remove '>'
      genes=c(genes,g)
      seqs=c(seqs, sequence)
      sequence=""
    }else{
      ### sequence
      sequence=paste0(sequence,line)
    }

  }
  seqs=c(seqs, sequence) ## last sequence
  #close(con)
  seqs=seqs[2:length(seqs)]
  dat=cbind(genes,seqs)
  return(dat)
}

dat=processFile("allupreg8_cluster2pfpvp_filtered_251021.txt")


cluster2<-c('TTATATAATA','AGGTAA','GCACTA')
for (i in 1:length(cluster2)){
  assign(paste0("frequencyMotif_",toString(i)), 0) ##  create a new varaible
  motif=tolower(cluster2[i])
  l.motif=nchar(motif)
  frequency=c()
  for(j in 1:dim(dat)[1]){
    count=0
    gene_seq=tolower(dat[j,2])
    l.gene_seq=nchar(gene_seq)
    for(ii in 1:((l.gene_seq-l.motif)+1)){
      test=substr(gene_seq, ii, (ii+l.motif-1))
      cat(test , " .... ", nchar(test),"\n")
      cat(motif , " .... ", nchar(motif),"\n")
      if(tolower(motif) == tolower(test)){
        count=count+1
      }
    }
    frequency=c(frequency, count)
  }
assign(paste0("frequencyMotif_",toString(i)), frequency)
}


freq_list <- list(frequencyMotif_1=frequencyMotif_1,
                  frequencyMotif_2=frequencyMotif_2,
                  frequencyMotif_3=frequencyMotif_3)


freqresultlist <- list()
for (f in 1:length(freq_list)) {
  freq_matrix = as.matrix(freq_list[[f]])
  colnames(freq_matrix) = "freq"
  motif_freq = as.data.frame(cbind(dat, freq_matrix))
 # name = paste("dat",freq_list[f],sep="")
  freqresultlist[[f]] <- motif_freq
  
} 

name <- as.character(names(freq_list))
names(freqresultlist) <- name


uni1=as.matrix(unique(freqresultlist[[1]][3]))
uni1=as.numeric(sort(uni1))

  
species = c("PBANKA","PCHAS", "PF3D7","PVX")  
#confirm total freq before separation
# fq=c()
# for ( i in 1:length(uni1)){
#   n=length(which(freqresultlist[[1]][3] == uni1[i]))
#   fq=c(fq,n)
# }

speciesfq1=c()
for (spec in species){
  count=subset(freqresultlist[[1]], startsWith(as.character(genes), spec))
 for (i in 1:length(uni1)){
n1=length(which(count[3] == uni1[i]))
speciesfq1=c(speciesfq1,n1)
}
  }

speciesfqdf1=data.frame(x=rep(uni1,2),y=speciesfq1,variable=rep(species, each = 4))

# plot
p1 <- ggplot(data = speciesfqdf1, aes(x=x, y=y)) + geom_line(aes(colour=variable)) +
  theme_classic() +
  theme(axis.title=element_text(size=10)) +
  ggtitle("TTATATAATA") +
  ylab("Frequency") + xlab("") +
  theme(legend.position = "none") 
 

uni2=as.matrix(unique(freqresultlist[[2]][3]))
uni2=as.numeric(sort(uni2))

#confirm total freq before separation
# fq=c()
# for ( i in 1:length(uni2)){
#   n=length(which(freqresultlist[[2]][3] == uni2[i]))
#   fq=c(fq,n)
# }

speciesfq2=c()
for (spec in species){
  count=subset(freqresultlist[[2]], startsWith(as.character(genes), spec))
  for (i in 1:length(uni2)){
    n2=length(which(count[3] == uni2[i]))
    speciesfq2=c(speciesfq2,n2)
  }
}

speciesfqdf2=data.frame(x=rep(uni2,2),y=speciesfq2,variable=rep(species, each = 4))

# plot
p2 <- ggplot(data = speciesfqdf2, aes(x=x, y=y)) + geom_line(aes(colour=variable)) +
  theme_classic() +
  theme(axis.title=element_text(size=10)) +
  ggtitle("AGGTAA") +
  ylab("") + xlab("Motif count per gene") +
  theme(legend.position = "none")
 


uni3=as.matrix(unique(freqresultlist[[3]][3]))
uni3=as.numeric(sort(uni3))

#confirm total freq before separation
# fq=c()
# for ( i in 1:length(uni3)){
#   n=length(which(freqresultlist[[3]][3] == uni3[i]))
#   fq=c(fq,n)
# }

speciesfq3=c()
for (spec in species){
  count=subset(freqresultlist[[3]], startsWith(as.character(genes), spec))
  for (i in 1:length(uni3)){
    n3=length(which(count[3] == uni3[i]))
    speciesfq3=c(speciesfq3,n3)
  }
}

speciesfqdf3=data.frame(x=rep(uni3,2),y=speciesfq3,variable=rep(species, each = 5))

# plot
p3 <- ggplot(data = speciesfqdf3, aes(x=x, y=y)) + geom_line(aes(colour=variable)) +
  #scale_x_continuous(breaks=c(0,1,2)) +
  theme_classic() +
  theme(axis.title=element_text(size=10)) +
  ggtitle("GCACTA") +
  ylab("") + xlab("") +
  theme(legend.position = "none")
 
 

plist <- list(p1,p2,p3)
motifplots <- grid.arrange(grobs = plist, ncol = 3, top = "Distribution of Cluster2 motif seeds")



