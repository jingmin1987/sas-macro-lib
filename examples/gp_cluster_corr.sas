*********************************************************************************************;
***                               Demo of gp_cluster_corr                                 ***;
*********************************************************************************************;

* Initialization;
%let train = /data2/users/jzhan162/POP_Consumer_USA_2016/2016/Var_Prescreen/Medium_FICO;
%include '/data2/users/jzhan162/macro_libs/ml_general.sas';

libname train "&train";

* train.diag_3 is an output of OUTSTATS from Proc Varclus;

* Default behavior;
%gp_cluster_corr(
    data_table  = train.diag_3,
    out_table   = check,
    delete_flag = 1
)

* Add other aggregating functions;
%gp_cluster_corr(
    data_table  = train.diag_3,
    out_table   = check,
    agg_funcs   = mean p50 p75 min max std,
    delete_flag = 1
)