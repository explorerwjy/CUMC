options(stringsAsFactors=F)

tab <- read.delim("Emily_Kinship_Putative_Table.txt")
gen <- read.table("../../../VCF_joint_calling/FreezeOne.genome", header=T)
kin <- read.delim("../../../VCF_joint_calling/Regeneron_Freeze_One.relatedness2")
sex <- read.table("../../../VCF_joint_calling/FreezeOne.sexcheck", header=T)
kintab <- read.delim("z:/Exome_Seq/resources/kintab.txt")

kinind <- apply(kin[,1:2], 1, function(x){paste(x[order(x)], collapse="")})
whi <- which(!duplicated(kinind))
kin <- kin[whi,]
rownames(kin) <- kinind[whi]

rownames(sex) <- sex[,2]

rownames(kintab) <- kintab[,1]

datind <- apply(tab[,2:3], 1, function(x){paste(x[order(x)], collapse="")})

sexlst <- c("Unidentified", "Male", "Female")

#get putative degree and kinship coefficient
dat <- cbind(tab, kintab[tab[,"Putative_Relationship"],"KinRet"], kin[datind,"RELATEDNESS_PHI"])
colnames(dat)[7:8] <- c("Putative_Degree", "Kinship_Coefficient")
if(any(is.na(dat[,"Putative_Degree"]))) dat[is.na(dat[,"Putative_Degree"]),"Putative_Degree"] <- ""
if(any(is.na(dat[,"Kinship_Coefficient"]))) dat[is.na(dat[,"Kinship_Coefficient"]),"Kinship_Coefficient"] <- ""
whi <- which(dat[,"Kinship_Coefficient"]!="")
# get inferred degree
dat <- cbind(dat, "")
colnames(dat)[9] <- "Inferred_Degree"
kincoef <- as.numeric(dat[whi,"Kinship_Coefficient"])
kincoef[kincoef<0] <- 0
dat[whi,"Inferred_Degree"] <- round(-log2(kincoef),1)-1
#check inferred matches expected
whi <- which(dat[,"Putative_Degree"]!=""&dat[,"Inferred_Degree"]!="")
kindiff <- abs(as.numeric(dat[whi,"Putative_Degree"])-as.numeric(dat[whi,"Inferred_Degree"]))
kindiff[is.nan(kindiff)] <- 0
kincheck <- rep("TRUE", length(whi))
kincheck[which(kindiff>0.45)] <- "?"
kincheck[which(kindiff>0.8)] <- "FALSE"
kincheck[which(dat[whi,"Putative_Degree"]==Inf&dat[whi,"Inferred_Degree"]!=Inf&as.numeric(dat[whi,"Inferred_Degree"])>5)] <- "?"
kincheck[which(dat[whi,"Putative_Degree"]==Inf&dat[whi,"Inferred_Degree"]!=Inf&as.numeric(dat[whi,"Inferred_Degree"])>9)] <- "TRUE"
dat <- cbind(dat, "")
colnames(dat)[10] <- "Check_Kinship"
dat[whi,"Check_Kinship"] <- kincheck

#get sexes and check them
dat <- cbind(dat, sexlst[sex[dat[,"Indv1"],"SNPSEX"]+1], sexlst[sex[dat[,"Indv2"],"SNPSEX"]+1])
colnames(dat)[11:12] <- c("Inferred_Sex1", "Inferred_Sex2")
dat <- cbind(dat, dat[,"Inferred_Sex1"]==dat[,"Sex1"], dat[,"Inferred_Sex2"]==dat[,"Sex2"])
colnames(dat)[13:14] <- c("Check_Sex1", "Check_Sex2")
if(any(dat[,"Sex1"]=="")) dat[dat[,"Sex1"]=="","Check_Sex1"] <- ""
if(any(dat[,"Sex2"]=="")) dat[dat[,"Sex2"]=="","Check_Sex2"] <- ""
for(i in 11:14){
  if(any(is.na(dat[,i]))) dat[is.na(dat[,i]),i] <- ""
}
if(any(dat[,"Inferred_Sex1"]=="Unidentified"))dat[dat[,"Inferred_Sex1"]=="Unidentified","Check_Sex1"] <- "Unknown"
if(any(dat[,"Inferred_Sex2"]=="Unidentified")) dat[dat[,"Inferred_Sex2"]=="Unidentified","Check_Sex2"] <- "Unknown"
dat[dat[,"Inferred_Degree"]==Inf,"Inferred_Degree"] <- "Unrelated"
write.table(dat, "$KinTab", col.names=T, row.names=F, quote=F, sep="\t")

#for mismatches get the plink IBD results
whi <- which(dat[,"Check_Kinship"]=="FALSE")
if(length(whi)>0){
  outgen <- vector()
  for(i in unique(dat[whi,"Family"])){
    mism <- unique(c(dat[dat[,"Family"]==i,"Indv1"],dat[dat[,"Family"]==i,"Indv2"]))
    outgen <- rbind(outgen, gen[gen[,1]%in%mism&gen[,3]%in%mism,])
  }
  mismind <- apply(dat[whi,c("Indv1", "Indv2")], 1, function(x){paste(x[order(x)], collapse="")})
  outgenind <- apply(outgen[,c(1,3)], 1, function(x){paste(x[order(x)], collapse="")})
  msim <- vector()
  msim[outgenind%in%mismind] <- "Mismatch"
  outgen <- cbind(outgen, msim)
  outgen[is.na(outgen[,"msim"]),"msim"] <- ""
  colnames(outgen)[15] <- "Mismatch"
  write.table(outgen, "MismatchedKinship_plinkIDB.txt", col.names=T, row.names=F, quote=F, sep="\t")
}
##chck the sex and output plots for mismatches
whi1 <- which(dat[,"Check_Sex1"]=="FALSE")
whi2 <- which(dat[,"Check_Sex2"]=="FALSE")
nsex <- matrix(nc=2, nr=0)
if(length(whi1)>0) nsex <- dat[whi1,c("Indv1", "Sex1")]
colnames(nsex) <- c("Indv", "ExpectedSex")
if(length(whi2)>0){
  temp <- dat[whi2,c("Indv2", "Sex2")]
  colnames(temp) <- c("Indv", "ExpectedSex")
  nsex <- rbind(nsex, temp)
}
nsex <- unique(nsex)
if(nrow(nsex)>0){
  for(i in 1:nrow(nsex)){
    pdf(paste("IncorrectSex_", nsex[i,1], ".pdf", sep=""))
     hdt <- hist(sex[,"F"], main=paste("Sample",nsex[i,1], "expected:", nsex[i,2]), xlab="", breaks=20)
     abline(v=sex[nsex[i,1],"F"], col="red")
     halfb <- floor(length(hdt\$breaks)/2)
     mtext("<<<------------------- Female", 1, 3, at=hdt\$breaks[halfb+3], adj=1)
     mtext("Male ------------------->>>", 1, 4, at=hdt\$breaks[halfb-1], adj=0)
     mtext("F-Statistic", 1, 2)
     mtext(nsex[i,1], 3, 0, at=sex[nsex[i,1],"F"], col="red")
    dev.off()
  }
}
