df1 <- read.table("gencode.v21.chr_patch_hapl_scaff.annotation.txt.gz", header=T)
head(df1)
df2 <- read.table("featureCount.Length.txt", header=T)
head(df2)

# clean df2
df_location <- df2
colnames(df_location)[1] = "Gene"
df_location$Chr <- gsub("\\;.+$", "", df_location$Chr)  # remove duplicated parts
df_location$Strand <- gsub("\\;.+$", "", df_location$Strand)  # remove duplicated parts
df_location$Start <- gsub("\\;.+$", "", df_location$Start)  # remove duplicated parts
df_location$End <- gsub("^.+\\;", "", df_location$End)  # remove duplicated parts
head(df_location)

# merge
df_out <- merge(df1, df_location, by=1)
head(df1)
head(df2)
head(df_out)
dim(df1)
dim(df2)
dim(df_out)
write.table(df_out, "gencode.v21.primary_assembly.anno.txt", row.names=F, quote=F)
