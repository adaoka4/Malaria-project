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

dat=processFile("allupreg8_cluster1pfpvp_filtered_251021.txt")


cluster1<-c("AGACAA","TGTTGTA","AATTTTTAAA","ATTACA","CAACCAA","TGCATGCA")
for (i in 1:length(cluster1)){
  assign(paste0("frequencyMotif_",toString(i)), 0) ##  create a new varaible
  motif=tolower(cluster1[i])
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
                  frequencyMotif_3=frequencyMotif_3,
                  frequencyMotif_4=frequencyMotif_4,
                  frequencyMotif_5=frequencyMotif_5,
                  frequencyMotif_6=frequencyMotif_6)


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
  #scale_x_continuous(breaks=c(0,1,2)) +
  theme_classic() +
  theme(axis.title=element_text(size=10)) +
  ggtitle("AGACAA") +
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

speciesfqdf2=data.frame(x=rep(uni2,2),y=speciesfq2,variable=rep(species, each = 5))

# plot
p2 <- ggplot(data = speciesfqdf2, aes(x=x, y=y)) + geom_line(aes(colour=variable)) +
 # scale_x_continuous(breaks=c(0,1,2)) +
  theme_classic() +
  theme(axis.title=element_text(size=10)) +
  ggtitle("TGTTGTA") +
  ylab("") + xlab("") +
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

speciesfqdf3=data.frame(x=rep(uni3,2),y=speciesfq3,variable=rep(species, each = 3))

# plot
p3 <- ggplot(data = speciesfqdf3, aes(x=x, y=y)) + geom_line(aes(colour=variable)) +
  scale_x_continuous(breaks=c(0,1,2)) +
  theme_classic() +
  theme(axis.title=element_text(size=10)) +
  ggtitle("AATTTTTAAA") +
  ylab("") + xlab("") +
  theme(legend.position = "none")
#labs(colour="Species") +
#theme(legend.title = element_text( size=10), legend.text=element_text(size=8))


uni4=as.matrix(unique(freqresultlist[[4]][3]))
uni4=as.numeric(sort(uni4))

#confirm total freq before separation
# fq=c()
# for ( i in 4:length(uni4)){
#   n=length(which(freqresultlist[[4]][3] == uni4[i]))
#   fq=c(fq,n)
# }

speciesfq4=c()
for (spec in species){
  count=subset(freqresultlist[[4]], startsWith(as.character(genes), spec))
  for (i in 1:length(uni4)){
    n4=length(which(count[3] == uni4[i]))
    speciesfq4=c(speciesfq4,n4)
  }
}

speciesfqdf4=data.frame(x=rep(uni4,2),y=speciesfq4,variable=rep(species, each = 7))

# plot
p4 <- ggplot(data = speciesfqdf4, aes(x=x, y=y)) + geom_line(aes(colour=variable)) +
  #scale_x_continuous(breaks=c(0,1,2)) +
  theme_classic() +
  theme(axis.title=element_text(size=10)) +
  ggtitle("ATTACA") +
  ylab("Frequency") + xlab("") +
  theme(legend.position = "none") 


uni5=as.matrix(unique(freqresultlist[[5]][3]))
uni5=as.numeric(sort(uni5))

#confirm total freq before separation
# fq=c()
# for ( i in 1:length(uni5)){
#   n=length(which(freqresultlist[[5]][3] == uni5[i]))
#   fq=c(fq,n)
# }

speciesfq5=c()
for (spec in species){
  count=subset(freqresultlist[[5]], startsWith(as.character(genes), spec))
  for (i in 1:length(uni5)){
    n5=length(which(count[3] == uni5[i]))
    speciesfq5=c(speciesfq5,n5)
  }
}

speciesfqdf5=data.frame(x=rep(uni5,2),y=speciesfq5,variable=rep(species, each = 4))

# plot
p5 <- ggplot(data = speciesfqdf5, aes(x=x, y=y)) + geom_line(aes(colour=variable)) +
  # scale_x_continuous(breaks=c(0,1,5)) +
  theme_classic() +
  theme(axis.title=element_text(size=10)) +
  ggtitle("CAACCAA") +
  ylab("") + xlab("Motif count per gene") +
  theme(legend.position = "none")



uni6=as.matrix(unique(freqresultlist[[6]][3]))
uni6=as.numeric(sort(uni6))

#confirm total freq before separation
# fq=c()
# for ( i in 1:length(uni6)){
#   n=length(which(freqresultlist[[6]][6] == uni6[i]))
#   fq=c(fq,n)
# }

speciesfq6=c()
for (spec in species){
  count=subset(freqresultlist[[6]], startsWith(as.character(genes), spec))
  for (i in 1:length(uni6)){
    n6=length(which(count[6] == uni6[i]))
    speciesfq6=c(speciesfq6,n6)
  }
}

speciesfqdf6=data.frame(x=rep(uni6,2),y=speciesfq6,variable=rep(species, each = 3))

# plot
p6 <- ggplot(data = speciesfqdf6, aes(x=x, y=y)) + geom_line(aes(colour=variable)) +
  scale_x_continuous(breaks=c(0,1,2)) +
  theme_classic() +
  theme(axis.title=element_text(size=10)) +
  ggtitle("TGCATGCA") +
  ylab("") + xlab("") +
  theme(legend.position = "none")
#labs(colour="Species") +
#theme(legend.title = element_text( size=10), legend.text=element_text(size=8))



plist <- list(p1,p2,p3,p4,p5,p6)
motifplots <- grid.arrange(grobs = plist, ncol = 3, top = "Distribution of Cluster1 motif seeds")



