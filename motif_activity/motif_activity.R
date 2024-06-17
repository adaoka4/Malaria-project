
##Npm  ---> motif counts per promotor
## Eps  ----> gene expression matrix
### "aggtaa_fpkm.csv"  "gcacta_fpkm.csv"  "npm-aggtaa.csv"  "npm-gcacta.csv"   
aggtaa=read.table("aggtaa_fpkm.csv", header=T, sep=',')
aggtaa_npm= read.table("npm-aggtaa.csv", header=T, sep='\t')
npm_aggtaa=as.data.frame(aggtaa_npm[,2])
Eps_aggtaa=aggtaa[,c(3:8)]
stages=colnames(Eps_aggtaa)

Eps_noramlize= function(df){
    row_means=c()
    col_means=c()
    for (i in  1:dim(df)[1]){
       row_means[i]=mean(as.numeric(df[i,])) 
    }
    for (i in  1:dim(df)[2]){
       col_means[i]=mean(as.numeric(df[,i]))
    }
    df_noramlize=df
     for (i in  1:dim(df_noramlize)[1]){
        df_noramlize[i,]=  df_noramlize[i,] - row_means[i]
    }
    for (i in  1:dim(df_noramlize)[2]){
       df_noramlize[,i]= df_noramlize[,i] - col_means[i]
    }

    return(df_noramlize)
}

Npm_noramlize= function(df){
    col_means=c()

    for (i in  1:dim(df)[2]){
       col_means[i]=mean(as.numeric(df[,i]))
    }
    df_noramlize=df
    for (i in  1:dim(df_noramlize)[2]){
       df_noramlize[,i]= df_noramlize[,i] - col_means[i]
    }

    return(df_noramlize)
}
## aggtaa
aggtaa_Eps= Eps_noramlize(Eps_aggtaa)
aggtaa_Npm=Npm_noramlize(npm_aggtaa)
aggtaa_motif_activity=c()
for(i in 1:dim(aggtaa_Eps)[2]){
    model=lm(aggtaa_Npm[,1]~aggtaa_Eps[,i])
    aggtaa_motif_activity[i]=coef(model)[2]

}
aggtaa_motif_mean=mean(aggtaa_motif_activity)
aggtaa_motif_activity_f= aggtaa_motif_activity-aggtaa_motif_mean
## gcacta
### "aggtaa_fpkm.csv"  "gcacta_fpkm.csv"  "npm-aggtaa.csv"  "npm-gcacta.csv"   
gcacta=read.table("gcacta_fpkm.csv", header=T, sep=',')
gcacta_npm= read.table("npm-gcacta.csv", header=T, sep='\t')
npm_gcacta=as.data.frame(gcacta_npm[,2])
Eps_gcacta=gcacta[,c(3:8)]

gcacta_Eps= Eps_noramlize(Eps_gcacta)
gcacta_Npm=Npm_noramlize(npm_gcacta)
gcacta_motif_activity=c()
for(i in 1: dim(gcacta_Eps)[2]){
    model=lm(gcacta_Npm[,1]~gcacta_Eps[,i])
    gcacta_motif_activity[i]=coef(model)[2]

}


gcacta_motif_mean=mean(gcacta_motif_activity)
gcacta_motif_activity_f= gcacta_motif_activity-gcacta_motif_mean

plot(aggtaa_motif_activity_f,type="l",col="red")
    par(new=TRUE)
plot(gcacta_motif_activity_f,type="l", col="green")

final_data= as.matrix(aggtaa_motif_activity_f)
final_data= cbind(final_data, as.matrix(gcacta_motif_activity_f))
colnames(final_data)=c("aggtaa", "gcacta")
png("Motif_activity.png", width = 600, height = 600)
par(mar = c(5, 5, 5, 5) + 0.3,bg = "#f7f7f7")              # Additional space for second y-axis
plot(final_data[,1], type="b", pch = 16, col = 2, ylab="Motif activity 'aggtaa'", 
xaxt = "n", xlab="Stages", lwd=2.0)              # Create first plot
par(new = TRUE)                             # Add new plot
plot(final_data[,2],type="b", pch = 17, col = 3,              # Create second plot without axes
     axes = FALSE, xlab = "", ylab = "", lwd=2.0)
axis(side = 4, at = pretty(range(final_data[,2])))      # Add second axis
mtext("Motif activity 'gcacta'", side = 4, line = 3)  
legend("topleft",
       legend = c("aggtaa", "gcacta"),
       col = 2:3,
       pch = 16:17) 
 axis(1, at = c(1:6),
     labels = stages)    
dev.off()
 rownames(final_data) = stages    

 write.table(final_data, file="motif_activity.txt", sep="\t")
