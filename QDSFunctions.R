suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(sctransform))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(glmnet))

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
  df = df %>%
    SCTransform(vst.flavor = "v2", variable.features.n=variable.features) %>%
    GetAssayData(slot="scale.data") %>% 
    as.data.frame()
  return(df)
}

filter_cells = function(df, nUMI_filt=500, nGene_filt=250, log10_filt=0.75, mito_filt=0.2) {
  filtered_seurat <- subset(x = df, 
                            subset= (nUMI >= nUMI_filt) & 
                              (nGene >= nGene_filt) & 
                              (log10GenesPerUMI > log10_filt) & 
                              (mitoRatio < mito_filt))
  return(filtered_seurat)
}

sc_preproc = function(df, prop=0.3, 
                      nUMI_filt=500, nGene_filt=250, log10_filt=0.75, mito_filt=0.2) {
    print("Before filter:")
    print(df)
    df = df %>% filter_cells(nUMI_filt=nUMI_filt, nGene_filt=nGene_filt,
                            log10_filt=log10_filt, mito_filt=mito_filt)
    print("After filter:")
    print(df)
    good_genes = df %>% get_good_genes(prop=prop)
    df = df %>% v2_sct() %>%
        rownames_to_column(var="Geneid") %>%
        filter(Geneid %in% good_genes)
    return(df) }

get_seurat_obj = function(df, mt.pattern="^MT-") {
  df = df %>% 
    as.matrix() %>% CreateSeuratObject()
  df$log10GenesPerUMI <- log10(df$nFeature_RNA) / log10(df$nCount_RNA)
  df$mitoRatio <- PercentageFeatureSet(object = df, pattern = mt.pattern)
  df$mitoRatio <- df@meta.data$mitoRatio / 100
  metadata <- df@meta.data
  metadata <- metadata %>%
    dplyr::rename(sample = orig.ident,
                  nUMI = nCount_RNA,
                  nGene = nFeature_RNA)
  df@meta.data <- metadata
  return(df)
}

make_hist = function(obj, metric, grouping=NULL, ann_df=NULL) {
    meta = obj@meta.data
    if (!is.null(ann_df)) {
        meta = meta %>% rownames_to_column(var="sample_id") %>%
            inner_join(ann_df, by=c("sample_id"=colnames(ann_df)[1])) %>%
            column_to_rownames(var="sample_id")
    }
    vals = meta[[metric]]
    ub = quantile(vals)[3]+1.5*IQR(vals)
    p = meta %>%
        ggplot(aes_string(x=metric, color = grouping, fill=grouping)) +
        geom_density(alpha = 0.2) +
        theme_classic() +
        xlim(NA, ub) +
      ggtitle(paste(metric, "\ngrouped by", grouping))
    return(p)
    }
                       

bulk_preproc = function(df) {
	good_genes = df %>% cpm_filter() %>% rownames()
	cond_df = data.frame(row.names=colnames(df), Condition=rep(1, ncol(df)))
	df = (df %>% normal_deseq(cond_df = cond_df)) %>%
		rownames_to_column(var="Geneid") %>%
		filter(Geneid %in% good_genes)
	return(df)
}

get_good_genes = function(obj, prop=0.1) {
  counts = obj %>% GetAssayData()  
  thresh = (ncol(counts)*prop) %>% round()  
  nonzero <- counts > 0
  keep_genes <- Matrix::rowSums(nonzero) >= thresh
  good_genes = rownames(counts)[keep_genes]
  return(good_genes)
}


convert_species = function(full_df, hm, hom_type, conf, old_id, new_id) {
  hm = hm %>% filter(!!sym(hom_type)=="ortholog_one2one" & !!sym(conf)==1) %>%
	select(c(old_id, new_id)) %>%
	distinct(!!sym(old_id), .keep_all=T) 		
  conv_df = inner_join(full_df, hm, by=c("Geneid"=old_id)) %>%
    select(-Geneid) %>%
    distinct(!!sym(new_id), .keep_all = T) %>%
    dplyr::rename(Geneid=!!sym(new_id)) %>%
    relocate(Geneid)
  return(conv_df)
}


# convert gene symbol to ensembl ID
convert_symbol = function(gene_map, full_df) {
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


run_ComBat = function(train_df, test_df) {
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
                   
alpha_test = function(x, y, n_alphas, nfolds=NULL, family="gaussian", type.measure="mse") {
    if (is.null(nfolds)) {
        nfolds=nrow(x)
    }
    x["y"] = y
    alphas = seq(0, 1, length.out=n_alphas)
    fit_list = list()
    err_vec = c()
    for (i in 1:length(alphas)) {
        alpha = alphas[i]
        fit = build_model(x, alpha=alpha, nfolds=nfolds, family=family, type.measure=type.measure)
        fit_list[[i]] = fit
        mce = min(fit$cvm)
        err_vec = c(err_vec, mce)
    }
    err_df = data.frame(Alpha=alphas, Errors=err_vec)
    out = setNames(list(err_df, fit_list), c("err", "fit"))
    return(out)
}

make_boxplot = function(pred_df, ann_df, grouping, title=NULL) {
  pred_df = inner_join(pred_df, ann_df, by=c(sample_id=colnames(ann_df)[1])) 
  p = ggplot(pred_df, aes_string(x=grouping, y="QDS")) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(shape=16, height=0, width=0.2) +
    ggtitle(title)
  return(p)
}

make_grouped_hist = function(pred_df, ann_df, grouping1, grouping2, title=NULL) {
    pred_df = inner_join(pred_df, ann_df, by=c(sample_id=colnames(ann_df)[1])) 
    modes = pred_df %>% group_by(!!sym(grouping1), !!sym(grouping2)) %>%
        group_split() %>%
        lapply(FUN=function(x) {
            cond = x %>% pull(!!sym(grouping1)) %>% unique()
            mode_idx = which.max(density(x$QDS)$y)
            mode = density(x$QDS)$x[mode_idx]
            df = data.frame(Condition=cond, QDS_mode=mode)
            colnames(df)[1] = grouping1
            return(df)
        })
    modes_df = bind_rows(modes) %>% group_by(!!sym(grouping1)) %>% summarise(QDS_mod_med=median(QDS_mode))

    p = ggplot(pred_df, aes(x=QDS, fill=NULL, group=!!sym(grouping2)))+
        geom_density(adjust=1.5, alpha=.4) +
        geom_vline(data=modes_df, mapping=aes(xintercept=QDS_mod_med, color="red")) +
        facet_wrap(as.formula(paste("~", grouping1)), ncol=1)
    return(p)
}
                   
                   
