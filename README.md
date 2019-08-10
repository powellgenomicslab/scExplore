# scExplore

R package for single-cell data visualization

## Instalation

```r
devtools::install_github("powellgenomicslab/scExplore")
```

## Examples

Plot merged gene expression signal to detect *CD8+ T cells*

```r
plotMerge(pbmc, gene1 = "CD3D", gene2 = "CD8A", size = 0.2)
```

![Example](misc/example.png)


Plot divided by cluster information

```r
plotMerge(pbmc, gene1 = "CD3D", gene2 = "CD8A", group = "seurat_clusters", size = 0.2)
```

![Example](misc/example_clusters.png)


Plot expression of `CD3D` in pca

```r
plotFeature(pbmc, "CD3D", reduction = "pca")
```

![Example](misc/example_expression.png)


Plot clusters

```r
plotFeature(pbmc, "seurat_clusters")
```

![Example](misc/example_clusters.png)
