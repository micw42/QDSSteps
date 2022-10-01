# QDS Workflow
## Build a linear regression model to predict quiescence depth
In this workflow, we build a linear regression model using an RNA-seq dataset with samples of known quiescence depth, then use the model to predict the quiescence depth of another RNA-seq dataset. The training data is a bulk RNA-seq dataset of rat embryonic fibroblasts (REF) serum starved for a range of 2 to 16 days. The test data is an scRNA-seq dataset of human lung cancer cells in various stages of development. We start with the train and test as raw counts, and we preprocess the data by filtering lowly expressed genes and (in the case of scRNA-seq) cells with low expression levels. Then, we normalize and VST the data. After that, we merge the train and test data and remove batch effects between them. Lastly, we build a linear regression model on the train data and use it to make predictions on the test data.

Load in the source file.
```
source("/xdisk/guangyao/michellewei/QDSWorkflow/QDSFunctions.R")
```

### Train data preprocessing (bulk RNA-seq)
Start with a table of raw counts, where rows are genes and columns are samples. Use the `bulk_preproc()` function to filter, normalize, and VST the data.
```
#Read in train data (raw counts)
train_raw = read.delim("/path/to/train/data") 

#Filter lowly expressed genes, normalize, and VST
train_df_norm = train_raw %>% bulk_preproc() 
```

### Test data preprocessing (scRNA-seq)
Start with a table of raw counts, where rows are genes and columns are samples. Use the `get_seurat_obj()` function to create a Seurat object from the table and calculate QC metrics. To process the test data, filter out cells and genes 
```
#Read in test data (raw counts)
test_raw = read.delim("/path/to/test/data")

#Create a Seurat object
obj = test_raw %>% get_seurat_obj()
```
Use the `make_hist()` function to make histograms of the QC metrics. Use the histograms to decide the QC metric thresholds to use to filter the cells. `ann_df` is an optional table with additional sample metadata that can be used to group the histograms. The first column of `ann_df` should be the sample names. To group the histograms, set `grouping` to the name of the desired `ann_df` column.
```
#Make histograms of QC metrics (nUMI, nGene, log10GenesPerUMI, and mitoRatio)
make_hist(obj, ann_df = NULL, grouping=NULL, metric="nUMI")
make_hist(obj, ann_df = NULL, grouping=NULL, metric="nGene")
make_hist(obj, ann_df = NULL, grouping=NULL, metric="log10GenesPerUMI")
make_hist(obj, ann_df = NULL, grouping=NULL, metric="mitoRatio")

```
Use the histograms to decide the QC metric thresholds to use to filter the cells. Use the `sc_preproc()` function to filter the cells by QC metric, filter the genes by raw count, and normalize and VST the data. `prop` is the minimum proportion of cells with raw count > 0 that a gene must have in order to pass the filter.
```
test_norm = obj %>% sc_preproc(prop=0.3, nUMI_filt=500, nGene_filt=250, log10_filt=0.75, mito_filt=0.2)

```

### Merge train and test 
To merge the train and test datasets, their gene IDs must match. If the train and test dataset use different species Ensembl IDs, convert the IDs to the same species using the `convert_species()` function. Specify the names of the homology type, confidence, original ID, and new ID columns of the homology table using the`hom_type`, `conf`, `old_id`, and `new_id` parameters, respectively. If a dataset uses gene symbols instead of IDs, convert the gene symbols to Ensembl IDs using the `convert_symbol()` function before converting to the correct species. For steps to get the `gene_map` and `hm` tables, see the chapter. 
```
gene_map = read.csv("/path/to/gene/map")
hm = read.delim("/path/to/homology/file")

test_conv = test_norm %>% convert_symbol(gene_map = gene_map) %>%
    convert_species(hm=hm, hom_type = "Rat.homology.type",
                conf="Rat.orthology.confidence..0.low..1.high.",
                old_id="Gene.stable.ID", new_id="Rat.gene.stable.ID")
```
Use the `run_ComBat()` function to inner join the train and test dataset and remove batch effects. This will output a named list with the new train and test tables.
```
out = run_ComBat(train_df, test_df)
```

### Build linear regression model
Use the `alpha_test()` function to test a sequence of 20 alpha values and get the corresponding cross-validation error. Set `y` as a vector of y values corresponding to the training samples. Y values should be numeric and correspond to the quiescence depth of the samples. For example, they can be the number of serum starvation days.
```
train = out$train
test = out$test
alpha_out = alpha_test(x=train, y=c("vector", "of", "y", "values"),
                        20, family="gaussian", type.measure="mse")
```
Get the model with the alpha value that gives the lowest cross-validation error. Use this model to predict the QDS of the test samples.
```
alpha_df = alpha_out$err
fit_list = alpha_out$fit
min_idx = which.min(alpha_df$Errors)
opt_fit = fit_list[[min_idx]]
pred_df = data.frame(sample_id=rownames(test),
        pred=predict(opt_fit, newx=as.matrix(test), s=opt_fit$lambda.min, type="response")) %>%
        rename(QDS=s1)
```

### Visualize results
Use the `make_boxplot()` function to make a boxplot of test sample QDS values. `ann_df` is a sample metadata table, where the first column contains sample names. Set `grouping` as the column of `ann_df` to use to group the boxplot.
```
make_boxplot(pred_df, ann_df, grouping="variable_to_group_by")
```
Use the `make_grouped_hist()` function to make a histogram of test sample QDS values grouped by 2 variables. This will output a series of histograms grouped by `grouping1`, with the density curves on each individual histogram grouped by `grouping2`.
```
make_grouped_hist(pred_df, ann_df, grouping1="1st_grouping_variable", grouping2="2nd_grouping_variable")
```








