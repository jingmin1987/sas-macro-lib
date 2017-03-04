# sas-macro-lib
A general purpose macr library for variable selection.

# Usage
macro library for: General purpose           
macro prefix     : gp_                       
macro list       : gp_info_value             
                   gp_stability_chart        
                   gp_stability_index        
                   gp_var_split              
                   gp_univariate_c           
                   gp_cluster_corr           
                   gp_cluster_size           
                                             
use %gp_help(macro_name) for more info. 

# Examples
## gp_cluster_corr
```sas
%gp_cluster_corr(
    data_table  = train.diag_3,
    out_table   = check,
    agg_funcs   = mean p50 p75 min max std,
    delete_flag = 1
)
```

### screenshot of the result table
![alt text][corr]
[corr]: https://github.com/jingmin1987/sas-macro-lib/blob/master/examples/gp_cluster_corr.png 'Cluster Correlation'

## gp_cluster_size
```sas
%gp_cluster_size(
    data_table = train.diag_3,
    delete_flag = 1,
    plot_flag = 1,
    log_flag = 1,
    plot_clus_list = 10 20 40 80 160 320 640,
    plot_num_cols = 4,
    plot_width = 5in,
    plot_height = 5in
)
```
### screenshot of histogram plot
![alt text][size_hist]
[size_hist]: https://github.com/jingmin1987/sas-macro-lib/blob/master/examples/gp_cluster_size_histogram.png 'Histogram Plot'

### screenshot of histogram plot
![alt text][size_den]
[size_den]: https://github.com/jingmin1987/sas-macro-lib/blob/master/examples/gp_cluster_size_density.png 'Density Plot'
=======
Handy SAS macros created for variable selection

# Purpose
This is part of my SAS library that I intended to showcase how a standardized SAS environment could look like. Each time when a new project is created, one should call %include "" to initialize the standard macro libraries similar to the ways in other general-purpose programming languages.
