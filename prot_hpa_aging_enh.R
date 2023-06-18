#2023.6.16
#source("/Users/dokada/Dropbox/analysis/2022.5/prot_hpa_aging_enh.R") #Done, 2023.6.17
out_path <- "/Users/dokada/Desktop/work/prot_hpa_aging_enh/"
if(file.exists(out_path)){
    unlink(out_path, recursive=TRUE)
    dir.create(out_path)
}else{
    dir.create(out_path)
}


#data read
lins <- read.csv("/Users/dokada/Desktop/work/aging_proteome_data/leha_protein_linear.csv", header = TRUE,row.names=1, stringsAsFactors = FALSE)
annot <- read.csv("/Users/dokada/Desktop/work/aging_proteome_data/leha_protein_annot.csv", header=TRUE, row.names=1, stringsAsFactors = FALSE)
all(sort(rownames(annot))==sort(rownames(lins)))
nt <- read.table("/Users/dokada/Desktop/work/aging_proteome_data/normal_tissue.tsv", header = TRUE, stringsAsFactors = FALSE, sep="\t",quote="")
id_table <- read.table("/Users/dokada/Desktop/work/aging_proteome_data/human_ens_ent3.txt", header = TRUE, stringsAsFactors = FALSE, sep="\t")
all_prots <- rownames(annot)
n_prot <- length(all_prots)

#ID table
res <- list()
for(i in 1:nrow(annot)){
    tmp_prot <- all_prots[i]
    tmp_ids <- annot[tmp_prot,"EntrezGeneID"]
    if(grepl(",", tmp_ids)){
        tmp_ids2 <- as.numeric(strsplit(tmp_ids,",")[[1]])
    }else{
        tmp_ids2 <- as.numeric(strsplit(tmp_ids," ")[[1]])
    }
    tmp <- NULL
    for(j in 1:length(tmp_ids2)){
        tmp_ens <- id_table[which(id_table[,2]==tmp_ids2[j]),"Gene.stable.ID"]
        if(length(tmp_ens)>0){
            tmp <- c(tmp, tmp_ens)
        }
    }
    res[[i]] <- tmp
}
names(res) <- all_prots
#res[["MBL2.3000.66.1"]]

#aging proteins for linear and nonlinear (tota; 1379, same with original paper)
lin_posi <- rownames(lins)[(lins[,"q.Age"] < 0.05) & (lins[,"Coefficient.Age"]>0)]
lin_nega <- rownames(lins)[(lins[,"q.Age"] < 0.05) & (lins[,"Coefficient.Age"]<0)]
geneset_posi <- rep(0, n_prot)
names(geneset_posi) <- all_prots
geneset_posi[lin_posi] <- 1
geneset_nega <- rep(0, n_prot)
names(geneset_nega) <- all_prots
geneset_nega[lin_nega] <- 1
lin <- geneset_posi + geneset_nega

#Calculate cell types by genes expression matrix
mat_ct <- read.csv("/Users/dokada/Desktop/work/prot_hpa_mat_enh/mat_ct.csv", header=TRUE, row.names=1, stringsAsFactors = FALSE)
mat_ts <- read.csv("/Users/dokada/Desktop/work/prot_hpa_mat_enh/mat_ts.csv", header=TRUE, row.names=1, stringsAsFactors = FALSE)
mat_list <- list(mat_ct, mat_ts)
tspt_all <- matrix(NA, n_prot, 2)
names(mat_list) <- c("ct", "ts")
thresh <- 1
for(m in 1:2){
    mat <- mat_list[[m]]
    unq_cat <- unique(rownames(mat))
    unq_genes <- unique(colnames(mat))
    tsig <- rep(NA, length(unq_genes))
    names(tsig) <- unq_genes
    tsig <- rep(NA, length(unq_genes))
    names(tsig) <- unq_genes
    for(i in 1:length(unq_genes)){
        vec <- mat[,unq_genes[i]]
        idxs_high <- which(vec==3)
        if(length(idxs_high)==1){
            if(all(vec[-idxs_high] < thresh)){
                tsig[i] <- unq_cat[idxs_high]
            }else{
                tsig[i] <- "Not_tissue_specific"
            }
        }else{
            tsig[i] <- "Not_tissue_specific"
        }
    }
    table(tsig)

    #mapping aging proteomic resul
    ts_info_prot <- rep(NA, n_prot)
    names(ts_info_prot) <- all_prots
    for(i in 1:n_prot){
        tmp_prot <- all_prots[i]
        tmp_ens <- res[[tmp_prot]]
        tmp_ens_ts <- tsig[tmp_ens]
        unq_ts_lab <- unique(tmp_ens_ts[!is.na(tmp_ens_ts)])
        if(length(unq_ts_lab)==0){
            ts_info_prot[i] <- "nohit"
        }else if(length(unq_ts_lab)==1){
            ts_info_prot[i] <- unq_ts_lab
        }else if(any(unq_ts_lab=="Not_tissue_specific")){
            ts_info_prot[i] <- "Not_tissue_specific"
        }else{
            cat("WARNING", i, k, tmp_prot, tmp_ens, tmp_ens_ts, "\n")
        }
    }
    table(tsig)
    tspt_all[,m] <- ts_info_prot
}
colnames(tspt_all) <- c("cell.type", "tissue")

#output
annot2 <- cbind(annot, geneset_posi, geneset_nega, lin, tspt_all)
write.csv(annot2, file=paste0(out_path, "annot_enh.csv"))
