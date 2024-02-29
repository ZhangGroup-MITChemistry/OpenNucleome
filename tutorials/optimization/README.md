## Optimization code

When developing the nucleus model, we optimized the interaction parameters between chromosomes and those between chromosome-nuclear landmarks separately. These optimizations make use of experimental data for HFF cells, and the resulting model is cell-type specific. When applying the model to a different cell type, one may consider reoptimizing the parameters with corresponding experimental data. The two folders here provide details on how to carry out parameter optimization.  

- chr_NL_optimization: parameter optimization for chromosome-nuclear landmarks interactions

- chr_chr_optimization: parameter optimization for chromosome-chromosome interactions
