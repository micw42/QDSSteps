library(dplyr)
library(tidyverse)
library(sctransform)
library(Seurat)
library(sva)
library(DESeq2)
library(purrr)

select = dplyr::select
rename = dplyr::rename
reduce = purrr::reduce

normal_deseq = function(full_df, cond_df) {
  dds <- DESeqDataSetFromMatrix(countData = full_df, colData = cond_df, design = ~ 1)
  vsd <- DESeq2::varianceStabilizingTransformation(dds)
  vsd = as.data.frame(assay(vsd))
  return(vsd)
}

#V2 norm
v2_sct = function(df, variable.features=3000) {
  df = df %>% as.matrix() %>%
    CreateSeuratObject() %>%
    SCTransform(vst.flavor = "v2", variable.features.n=variable.features) %>%
    GetAssayData(slot="scale.data") %>% 
    as.data.frame()
  return(df)
}

filter_cells = function(df, nUMI_filt=500, nGene_filt=250, log10_filt=0.75, mito_filt=0.2) {
  df = df %>% get_seurat_obj()
  filtered_seurat <- subset(x = df, 
                            subset= (nUMI >= nUMI_filt) & 
                              (nGene >= nGene_filt) & 
                              (log10GenesPerUMI > log10_filt) & 
                              (mitoRatio < mito_filt))
  return(filtered_seurat)
}

filter_pipeline = function(df, nUMI_filt=500, nGene_filt=250, log10_filt=0.75, mito_filt=0.2) {
  df = df %>% 
    filter_cells(nUMI_filt = nUMI_filt, nGene_filt = nGene_filt, log10_filt = log10_filt, mito_filt = mito_filt) %>%
    GetAssayData() %>%
    as.data.frame()
  return(df)
}

get_seurat_obj = function(df) {
  df = df %>% 
    as.matrix() %>% CreateSeuratObject()
  df$log10GenesPerUMI <- log10(df$nFeature_RNA) / log10(df$nCount_RNA)
  df$mitoRatio <- PercentageFeatureSet(object = df, pattern = "^MT-")
  df$mitoRatio <- df@meta.data$mitoRatio / 100
  metadata <- df@meta.data
  metadata <- metadata %>%
    dplyr::rename(sample = orig.ident,
                  nUMI = nCount_RNA,
                  nGene = nFeature_RNA)
  df@meta.data <- metadata
  return(df)
}

filter_cells = function(df, nUMI_filt=500, nGene_filt=250, log10_filt=0.75, mito_filt=0.2) {
  df = df %>% get_seurat_obj()
  filtered_seurat <- subset(x = df, 
                            subset= (nUMI >= nUMI_filt) & 
                              (nGene >= nGene_filt) & 
                              (log10GenesPerUMI > log10_filt) & 
                              (mitoRatio < mito_filt))
  return(filtered_seurat)
}

filter_genes = function(counts) {
  # Output a logical vector for every gene on whether the more than zero counts per cell
  nonzero <- counts > 0
  
  # Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
  keep_genes <- Matrix::rowSums(nonzero) >= 10
  
  # Only keeping those genes expressed in more than 10 cells
  filtered_counts <- counts[keep_genes, ]
  
  return(filtered_counts)
}

bulk_preproc = function(df) {
	good_genes = df %>% cpm_filter() %>% rownames()
	cond_df = data.frame(row.names=colnames(df), Condition=rep(1, ncol(df)))
	df = (df %>% normal_deseq(cond_df = cond_df)) %>%
		rownames_to_column(var="Geneid") %>%
		filter(Geneid %in% good_genes)
	return(df)
}

get_good_genes = function(counts, prop=0.1) {
  thresh = (ncol(counts)*prop) %>% round()  
  nonzero <- counts > 0
  keep_genes <- Matrix::rowSums(nonzero) >= thresh
  good_genes = rownames(counts)[keep_genes]
  return(good_genes)
}


convert_species = function(full_df, hm) {
  hm = hm %>% filter(Rat.homology.type=="ortholog_one2one" & Rat.orthology.confidence..0.low..1.high.==1) %>%
	select(Gene.stable.ID, Rat.gene.stable.ID)		
  musc_df = inner_join(full_df, hm, by=c("Geneid"="Gene.stable.ID")) %>%
    select(-Geneid) %>%
    distinct(Rat.gene.stable.ID, .keep_all = T) %>%
    dplyr::rename(Geneid=Rat.gene.stable.ID) %>%
    relocate(Geneid)
  return(musc_df)
}


# convert gene symbol to ensembl ID
id_convert = function(gene_map, full_df) {
  gene_map = gene_map %>%
    select(converted_alias, initial_alias) %>%
    filter(converted_alias!="None") %>%
    distinct(initial_alias, .keep_all=T)
  
  full_df = inner_join(full_df, gene_map, by=c("Geneid"="initial_alias")) %>%
    select(-Geneid) %>%
    dplyr::rename(Geneid=converted_alias) %>%
    relocate(Geneid) %>%
    distinct(Geneid, .keep_all=T)
  
  return(full_df)
}

reformat_TMS = function(TMS_file) {
  musc_df = readRDS(TMS_file) %>%
    column_to_rownames(var="index") %>%
    t() %>%
    as.data.frame() 
  musc_df = musc_df %>% rownames_to_column(var="Geneid")
  musc_df$Geneid = musc_df$Geneid %>% sapply(FUN=function(x) {
    if (startsWith(x, "X") & endsWith(x, "Rik")) {
      return(substring(x, 2))
    } else {
      return(x)
    }
  })
  saveRDS(musc_df, TMS_file)
}

filter_TMS = function(musc_df, gene_map, nUMI_filt=500, nGene_filt=250, 
                      log10_filt=0.75, mito_filt=0.2) {
  musc_df = musc_df %>% 
    column_to_rownames(var="Geneid") %>%
    filter_cells(nUMI_filt=nUMI_filt, nGene_filt=nGene_filt, 
                                     log10_filt=log10_filt, mito_filt=mito_filt) %>%
    GetAssayData() %>%
    as.data.frame() %>%
    rownames_to_column(var="Geneid") %>% 
    relocate(Geneid) %>%
    mutate(Geneid=toupper(Geneid))
  musc_df = id_convert(gene_map, musc_df) %>%
    distinct(Geneid, .keep_all = T) %>%
    column_to_rownames(var="Geneid")
  return(musc_df)
}

run_ComBat = function(train_df, test_df, batch) {
    n_train = ncol(train_df) - 1
    n_test = ncol(test_df) - 1
    train_batch = rep("train", n_train)
    test_batch = rep("test", n_test)
    batch = c(train_batch, test_batch)
    
    full_df = inner_join(train_df, test_df, by="Geneid") %>%
        column_to_rownames(var="Geneid")
    combat_edata = ComBat(full_df, batch=batch, mean.only = F)
    train_df = as.data.frame(combat_edata[,1:n_train] %>% t())
    test_df = as.data.frame(combat_edata[,(n_train+1):ncol(combat_edata)] %>% t())
    out = setNames(list(train_df, test_df), c("train", "test"))
    return(out)
}

cpm_filter = function(df, cell_thresh=2, count_thresh=0.5) {
  cpm_df = cpm_norm(df)
  idx = rowSums( cpm_df >= count_thresh ) >= cell_thresh
  df = df[idx,]
  return(df)
}

cpm_norm = function(df) {
  cpm_df = df %>% scale(center = FALSE,
                        scale = colSums(df)/10e6)
  return(cpm_df)
}

y_weights <- function(xy){
  cts <- xy %>% dplyr::count(y)
  tot <- sum(cts$n)
  cts$Weight <- 1 - cts$n / tot; cts$n <- NULL; names(cts) <- c("y", "Weight")

  wdf <- inner_join(data.frame(y=xy$y), cts, by="y")
  weights <- wdf$Weight

  return(weights)
}

build_model <- function(xy, alpha, nfolds = 5, family="binomial", type.measure="class", weighted=T, standardize=T){
  xy <- xy %>% arrange(y)
  x <- as.matrix(xy %>% select(-y))
  y <- as.matrix(xy %>% select(y))

  if (weighted){weights <- y_weights(xy)}
  else{weights <- NULL}

  fit <- cv.glmnet(x, y, nfolds=nfolds, family=family, type.measure=type.measure, weights=weights, alpha=alpha, standardize=standardize)

  return(fit)
}
                   
alpha_test = function(xy, n_alphas, nfolds=NULL, family="gaussian", type.measure="mse") {
    if (is.null(nfolds)) {
        nfolds=nrow(xy)
    }
    alphas = seq(0, 1, length.out=n_alphas)
    fit_list = list()
    err_vec = c()
    for (i in 1:length(alphas)) {
        alpha = alphas[i]
        fit = LOOCV(xy, alpha=alpha, nfolds=nfolds, family=family, type.measure=type.measure)
        fit_list[[i]] = fit
        mce = min(fit$cvm)
        err_vec = c(err_vec, mce)
    }
    err_df = data.frame(Alpha=alphas, Errors=err_vec)
    out = setNames(list(err_df, fit_list), c("err", "fit"))
    return(out)
}
                   
                   
