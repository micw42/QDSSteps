# QDS Workflow
## Build a linear regression model to predict quiescence depth
In this workflow, we build a linear regression model using an RNA-seq dataset with samples of known quiescence depth, then use the model to predict the quiescence depth of another RNA-seq dataset. The training data is a bulk RNA-seq dataset of rat embryonic fibroblasts (REF) serum starved for a range of 2 to 16 days. The test data is an scRNA-seq dataset of human lung cancer cells in various stages of development. We start with the train and test as raw counts, and we preprocess the data by filtering lowly expressed genes and (in the case of scRNA-seq) cells with low expression levels. Then, we normalize and VST the data. After that, we merge the train and test data and remove batch effects between them. Lastly, we build a linear regression model on the train data and use it to make predictions on the test data.

### Train data preprocessing (bulk RNA-seq)
Load in source file
```
source("QDSFunctions.R")
```
To preprocess the training data, subset the genes to those with CPM > 0.5 in at least 2 samples. Then, normalize and VST the data using DESeq2.
```
#Read in train data (raw counts)
train_raw = read.delim("/path/to/train/data") 

#Filter lowly expressed genes
train_df_norm = train_raw %>% distinct(Geneid, .keep_all=T) %>%
    column_to_rownames(var="Geneid") %>%
    bulk_preproc() 
```

### Test data preprocessing (scRNA-seq)
To process the test data, filter out cells and genes that have low expression values. Then, normalize and VST the data using Seurat.
```
#Read in test data (raw counts)
test_raw = read.delim("/path/to/test/data")

#Filter cells
test_filt = test_raw %>% distinct(Geneid, .keep_all=T) %>%
    column_to_rownames(var="Geneid") %>%
    filter_pipeline()

#Filter genes
good_genes = test_norm %>% get_good_genes(prop=0.3)
test_norm = test_filt %>% v2_sct() %>%
    rownames_to_column(var="Geneid") %>%
    filter(Geneid %in% good_genes)
```

### Merge train and test 
The test dataset uses human gene names, while the train dataset uses rat Ensembl IDs. Convert the human gene names in the test dataset to human Ensembl IDs using a gene map file from gProfiler, then convert the human Ensembl IDs to rat Ensembl IDs using a homology file from Ensembl.
```
gene_map = read.csv("/path/to/gene/map")
hm = read.delim("/path/to/homology/file")

test_conv = test_norm %>% id_convert(gene_map = gene_map) %>%
    convert_species(hm=hm)
```
Inner join the train and test dataset and remove batch effects using the sva package.
```
train_df = readRDS("/xdisk/guangyao/michellewei/arrayValTest3/REF_rat_d2d16_QDSInput.rds")
test_df = readRDS("/xdisk/guangyao/michellewei/epithelial_cells_conv.rds")

out = run_ComBat(train_df, test_df)
```

### Build linear regression model
Add a column of y values to the train data.
```
train = out$train
test = out$test

train$y = rep(c(2:4, seq(6, 16, by=2)), each=3)
```
Test a sequence of 20 alpha values to find the optimal one.
```
alpha_out = alpha_test(train, 20, family="gaussian", type.measure="mse")
alpha_df = out$err
fit_list = out$fit
min_idx = which.min(alpha_df$Errors)
opt_fit = fit_list[[min_idx]]
```
Use the model with the optimal alpha to predict the QDS of the test samples.
```
pred_df = data.frame(sample_id=rownames(test),
        pred=predict(opt_fit, newx=as.matrix(test), s=opt_fit$lambda.min, type="response")) %>%
        rename(QDS=s1)
```

### Visualize results
Merge the prediction table with sample metadata table.
```
ann_df = read.delim("/path/to/sample/metadata/")
pred_df = inner_join(pred_df, ann_df, by="sample_id")
```
Make a boxplot of test sample QDS values (replace Condition with relevant variable name).
```
ggplot(pred_df, aes(x = Condition, y = QDS)) + geom_boxplot(outlier.shape = NA) + 
        geom_jitter(shape = 16, height = 0, width = 0.2) 
```








