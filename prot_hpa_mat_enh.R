#Creat mateix from HPAdataset
#source("/Users/dokada/Dropbox/analysis/2022.5/prot_hpa_mat_enh.R") #Done, 2023.6.17
out_path <- "/Users/dokada/Desktop/work/prot_hpa_mat_enh/"
if(file.exists(out_path)){
    unlink(out_path, recursive=TRUE)
    dir.create(out_path)
}else{
    dir.create(out_path)
}
nt <- read.table("/Users/dokada/Desktop/work/aging_proteome_data/normal_tissue.tsv", header = TRUE, stringsAsFactors = FALSE, sep="\t",quote="")

#Edit nt matrix
nt_app <- nt[nt[,"Reliability"] %in% c("Enhanced"),]
level_num <- rep(NA, nrow(nt_app))
level_num[nt_app$Level=="Not detected"] <- 0
level_num[nt_app$Level=="Low"] <- 1
level_num[nt_app$Level=="Medium"] <- 2
level_num[nt_app$Level=="High"] <- 3
level_num[nt_app$Level %in% c("Ascending","Descending", "Not representative")] <- 1
nt_app$Level <- level_num

#Calculate cell types by genes expression matrix
unq_ct <- unique(nt_app[,"Cell.type"])
unq_genes <- unique(nt_app[,"Gene"])
mat <- matrix(NA, length(unq_ct), length(unq_genes))
rownames(mat) <- unq_ct
colnames(mat) <- unq_genes
for(i in 1:length(unq_ct)){
    tmp_ct <- unq_ct[i]
    for(j in 1:length(unq_genes)){
        tmp_gene <- unq_genes[j]
        idxs <- intersect(which(nt_app[,"Cell.type"]==tmp_ct), which(nt_app[,"Gene"]==tmp_gene))
        if(length(idxs)==0){
            mat[i,j] <- -1
        }else{
            mat[i,j] <- max(nt_app[idxs, "Level"])
        }
    }
}
write.csv(mat, file=paste0(out_path,"mat_ct.csv"))

#Calculate Tissues by genes expression matrix
unq_tissues <- unique(nt_app[,"Tissue"])
unq_genes <- unique(nt_app[,"Gene"])
mat_tissue <- matrix(NA, length(unq_tissues), length(unq_genes))
rownames(mat_tissue) <- unq_tissues
colnames(mat_tissue) <- unq_genes
for(i in 1:length(unq_tissues)){
    tmp_tissue <- unq_tissues[i]
    for(j in 1:length(unq_genes)){
        tmp_gene <- unq_genes[j]
        idxs <- intersect(which(nt_app[,"Tissue"]==tmp_tissue), which(nt_app[,"Gene"]==tmp_gene))
        if(length(idxs)==0){
            mat_tissue[i,j] <- -1
        }else{
            mat_tissue[i,j] <- max(nt_app[idxs, "Level"])
        }
    }
}
write.csv(mat_tissue, file=paste0(out_path,"mat_ts.csv"))

