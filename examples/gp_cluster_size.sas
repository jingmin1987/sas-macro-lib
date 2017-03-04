*********************************************************************************************;
***                               Demo of gp_cluster_size                                 ***;
*********************************************************************************************;

* Initialization;
%let train = /data2/users/jzhan162/POP_Consumer_USA_2016/2016/Var_Prescreen/Medium_FICO;
%include '/data2/users/jzhan162/macro_libs/ml_general.sas';

libname train "&train";

* train.diag_3 is an output of OUTSTATS from Proc Varclus;

* Default behavior with plot;
%gp_cluster_size(
    data_table = train.diag_3,
    delete_flag = 1,
    plot_flag = 1
)

* Log scale on cluster size (recommended);
%gp_cluster_size(
    data_table = train.diag_3,
    log_flag = 1,
    delete_flag = 1,
    plot_flag = 1
)

* Customize at what number to plot;
%gp_cluster_size(
    data_table = train.diag_3,
    delete_flag = 1,
    plot_flag = 1,
    log_flag = 1,
    plot_clus_list = 10 20 40 80 160 320 640
)

* Customize subplot layout and plot size, actual plot in github;
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

