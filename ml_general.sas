options nosource;
%put ****************************************************************************************;
%put ***     macro library for: General purpose                                           ***;   
%put ***     macro prefix     : gp_                                                       ***;
%put ***     macro list       : gp_info_value                                             ***;
%put ***                        gp_stability_chart                                        ***;
%put ***                        gp_stability_index                                        ***;
%put ***                        gp_var_split                                              ***;
%put ***                        gp_univariate_c                                           ***;
%put ***                        gp_cluster_corr                                           ***;
%put ***                        gp_cluster_size                                           ***;
%put ***                        gp_cluster_struct                                         ***;
%put ***                        gp_cluster_add_vars                                       ***;        
%put ***                                                                                  ***;
%put ***     use macro gp_help(macro_name) for more info.                                 ***;
%put ****************************************************************************************;

*********************************************************************************************;
***                                    Public Macros                                      ***;
*********************************************************************************************;
%macro gp_info_value(
                    data_table  = /* Source table name                                      */,
                    out_table   = /* Output table name                                      */,
                    n_threshold = /* Nlevel threshold to determine a variable if continuous */,
                    perf_var    = /* Performance variable such AS ODEFIND                   */,
                    date_var    = /* If given, also calculate IV over time by date_var      */,
                    graph_out   = /* Graph output if date_var is given. 1 for on and        */
                                  /* missing or 0 for off                                   */,
                    delete_flag = /* 1 or 0(miss). Delete all temporary data after running  */,
                    special_woe = 9.999
                                  /* Special value of WOE when sum(good) or sum(bad) = 0    */,
                    trans_data  = /* Dataset containing monotonically transformed data by   */
                                  /* proc transreg.                                         */

                        /* Purpose: calculates IV value based on WOE method                 */
                        /* Example: gp_info_value(data_table  = sample,                     */
                        /*                         out_table   = example,                   */
                        /*                         n_threshold = 20                         */
                        /*                         perf_var    = ODEFIND                    */
                        /*          )                                                       */
);

    %local func_num;

    * Initialization;
    %if %length(&date_var) = 0 %then %do;
        %let graph_out = 0;
    %end;
    %let discrete_vars =;
    %let continuous_vars =;
    %let transvars =;

    %let perf_var = %lowcase(&perf_var);
    %let date_var = %lowcase(&date_var);

    * Prefix table names for better clarification;
    %let func_num = gp_01;

    proc datasets lib = work nolist;
        delete &func_num.:;
    quit;

    * Make a small dataset for initial screening;
    data &func_num._sample;
        set &data_table;
        if ranuni(999) < 0.01; /* 1% sample */
    run;

    * Find if nlevles > n_threshold;
    ods select nlevels;
    ods output nlevels = &func_num._var_levels;

    proc freq data = &func_num._sample nlevels; run;

    proc sql noprint;
        SELECT tablevar INTO: discrete_vars SEPARATED BY ' '
            FROM &func_num._var_levels
            WHERE nlevels <= &n_threshold
        ;
    quit;

    ods select nlevels;
    ods output nlevels = &func_num._var_levels_u;

    proc freq data = &data_table.(keep = &discrete_vars) nlevels; run;

    * Find variable types num or char;
    ods select variables;
    ods output variables  = &func_num._var_type;

    proc contents data = &func_num._sample; run;

    ods select all;

    proc sort data = &func_num._var_levels; by tablevar; run;
    proc sort data = &func_num._var_levels_u; by tablevar; run;

    data &func_num._var_levels;
        update &func_num._var_levels &func_num._var_levels_u;
        by tablevar;
        rename tablevar = var_name;
    run;

    proc sql noprint;
        CREATE TABLE &func_num._var_check AS
            SELECT a.var_name,
                   a.nlevels,
                   b.type
            FROM &func_num._var_levels AS a
            FULL JOIN
                 &func_num._var_type AS b
            ON a.var_name = b.variable
        ;

        SELECT var_name INTO: discrete_vars SEPARATED BY ' '
            FROM &func_num._var_check
            WHERE (nlevels <= &n_threshold OR lowcase(type) = 'char')
                  AND lowcase(var_name) NOT IN ("&perf_var", "&date_var")
        ;

        SELECT var_name INTO: continuous_vars SEPARATED BY ' '
            FROM &func_num._var_check
            WHERE nlevels > &n_threshold AND lowcase(type) = 'num'
                  AND lowcase(var_name) NOT IN ("&perf_var", "&date_var")
        ;
    quit;

    data &func_num._disc_data(keep = &discrete_vars &perf_var &date_var)
         &func_num._cont_data(keep = &continuous_vars &perf_var &date_var)
         &func_num._perf_data(keep = &perf_var &date_var);
        set &data_table end = flag;

        good + (&perf_var = 0);
        bad + (&perf_var = 1);

        if flag then do;
            call symput('total_good', good);
            call symput('total_bad', bad);
        end;
    run;

    * Smooth data monotonically;
    %if %length(&continuous_vars) %then %do;
        %if %length(&trans_data) %then %do;
            %let var_keep = T%sysfunc(tranwrd(&continuous_vars, %str( ), %str( T)));
            data &func_num._cont_data_T;
                set &trans_data;
                keep &perf_var &var_keep;
            run;
        %end;
        %else %do;
            proc transreg data = &func_num._cont_data noprint;
                model identity(&perf_var) = monotone(&continuous_vars);
                output out = &func_num._cont_data_T;
            run;
        %end;

        * Select variables that are successfully transformed or in the trans_data;
        ods select variables;
        ods output variables  = &func_num._tvars;

        proc contents data = &func_num._cont_data_T; run;

        ods select all;

        proc sql noprint;
            SELECT variable INTO: transvars SEPARATED BY ' '
                FROM &func_num._tvars
                WHERE variable LIKE 'T%'
            ;
        quit;
    %end;


    * In-memory loading;
    %if %length(&discrete_vars) %then %do;
        sasfile &func_num._disc_data open;
    %end;
    %if %length(&continuous_vars) %then %do;
        sasfile &func_num._cont_data_T open;
    %end;

    * Calculate base_ratio that is used to calculate WOE for each var;
    proc sql noprint exec;
        SELECT sum(&perf_var = 0) / sum(&perf_var = 1) INTO: base_ratio
            FROM &func_num._perf_data
        ;

        * For discrete vars;
        %let i = 1;
        %do %while(%length(%scan(&discrete_vars, &i)) > 0);
            %let var = %scan(&discrete_vars, &i);
            %put Running WOE and IV calculation for var: &var;

            CREATE TABLE &func_num._woe_d_&i as
                SELECT DISTINCT
                       case when (sum(&perf_var = 0) = 0) then -&special_woe
                            when (sum(&perf_var = 1) = 0) then +&special_woe
                            else log(sum(&perf_var = 0) / sum(&perf_var = 1)) - log(&base_ratio) 
                       end AS WOE,
                       sum(&perf_var = 0) / &total_good AS dist_good,
                       sum(&perf_var = 1) / &total_bad AS dist_bad,
                       calculated dist_good - calculated dist_bad AS dist_diff,
                       calculated dist_diff * calculated WOE AS IV,
                       &var,
                       "&var" AS var_name length = 36
                FROM &func_num._disc_data
                group by &var
            ;
            %let i = %eval(&i + 1);
        %end;

        * For continuous vars;
        %let i = 1;
        %do %while(%length(%scan(&continuous_vars, &i)) > 0 
                   and %index(&transvars, T%scan(&continuous_vars, &i)));
            %let var = %scan(&continuous_vars, &i);
            %put Running WOE and IV calculation for var: &var;

            CREATE TABLE &func_num._woe_c_&i as
                SELECT DISTINCT
                       case when (sum(&perf_var = 0) = 0) then -&special_woe
                            when (sum(&perf_var = 1) = 0) then +&special_woe
                            else log(sum(&perf_var = 0) / sum(&perf_var = 1)) - log(&base_ratio) 
                       end AS WOE,
                       sum(&perf_var = 0) / &total_good AS dist_good,
                       sum(&perf_var = 1) / &total_bad AS dist_bad,
                       calculated dist_good - calculated dist_bad AS dist_diff,
                       calculated dist_diff * calculated WOE AS IV,
                       T&var,
                       "&var" AS var_name length = 36
                FROM &func_num._cont_data_T
                group by T&var
            ;
            %let i = %eval(&i + 1);
        %end;
    quit;

    * In-memory unloading;
    %if %length(&discrete_vars) %then %do;
        sasfile &func_num._disc_data close;
    %end;
    %if %length(&continuous_vars) %then %do;
        sasfile &func_num._cont_data_T close;
    %end;

    data &func_num._WOE;
        set &func_num._woe_:;
    run;

    proc means data = &func_num._WOE nway noprint;
        class var_name;
        var IV;
        output out = &func_num._iv(drop = _:) sum =;
    run;

    * If date_var is provided, then calculate IV by date_var AS well;
    %if %length(&date_var) %then %do;
        proc datasets lib = work nolist;
            delete &func_num._woe_:;
        quit;

        %if %length(&continuous_vars) %then %do;
            data &func_num._cont_data_T;
                set &func_num._cont_data_T;
                set &func_num._cont_data(keep = &date_var);
            run;
        %end;

        * In-memory loading;
        %if %length(&discrete_vars) %then %do;
            sasfile &func_num._disc_data open;
        %end;
        %if %length(&continuous_vars) %then %do;
            sasfile &func_num._cont_data_T open;
        %end;

        proc sql noprint exec;
            CREATE TABLE &func_num._count_by_date AS
                SELECT sum(&perf_var = 0) AS total_good,
                       sum(&perf_var = 1) AS total_bad,
                       (calculated total_good + 1e-3) / (calculated total_bad + 1e-3) AS base_ratio,
                       &date_var
                FROM &func_num._perf_data
                GROUP BY &date_var
                ORDER BY &date_var
            ;

            * For discrete vars;
            %let i = 1;
            %do %while(%length(%scan(&discrete_vars, &i)) > 0
                       and %scan(&discrete_vars, &i) ne &date_var);

                %let var = %scan(&discrete_vars, &i);
                %put Running WOE and IV calculation for var: &var;

                CREATE TABLE &func_num._woe_d_&i AS
                    SELECT DISTINCT
                           CASE WHEN (sum(&perf_var = 0) = 0) THEN -&special_woe
                                WHEN (sum(&perf_var = 1) = 0) THEN +&special_woe
                                ELSE log(sum(&perf_var = 0) / sum(&perf_var = 1)) - log(base_ratio) 
                           END AS WOE,
                           sum(&perf_var = 0) / total_good AS dist_good,
                           sum(&perf_var = 1) / total_bad AS dist_bad,
                           CALCULATED dist_good - CALCULATED dist_bad AS dist_diff,
                           CALCULATED dist_diff * CALCULATED WOE AS IV,
                           a.&date_var,
                           &var,
                           "&var" AS var_name length = 36
                    FROM &func_num._disc_data AS a,
                         &func_num._count_by_date AS b
                    WHERE a.&date_var = b.&date_var
                    GROUP BY a.&date_var, &var
                    ORDER BY a.&date_var, &var
                ;
                %let i = %eval(&i + 1);
            %end;

            * For continuous vars;
            %let i = 1;
            %do %while(%length(%scan(&continuous_vars, &i)) > 0 
                       and %index(&transvars, T%scan(&continuous_vars, &i))
                       and %scan(&continuous_vars, &i) ne &date_var);

                %let var = %scan(&continuous_vars, &i);
                %put Running WOE and IV calculation for var: &var;

                CREATE TABLE &func_num._woe_c_&i as
                    SELECT DISTINCT
                           CASE WHEN (sum(&perf_var = 0) = 0) THEN -&special_woe
                                WHEN (sum(&perf_var = 1) = 0) THEN +&special_woe
                                ELSE log(sum(&perf_var = 0) / sum(&perf_var = 1)) - log(base_ratio) 
                           END AS WOE,
                           sum(&perf_var = 0) / total_good AS dist_good,
                           sum(&perf_var = 1) / total_bad AS dist_bad,
                           CALCULATED dist_good - CALCULATED dist_bad AS dist_diff,
                           CALCULATED dist_diff * CALCULATED WOE AS IV,
                           a.&date_var,
                           T&var,
                           "&var" AS var_name length = 36
                    FROM &func_num._cont_data_T AS a,
                         &func_num._count_by_date AS b
                    WHERE a.&date_var = b.&date_var
                    GROUP BY a.&date_var, T&var
                    ORDER BY a.&date_var, T&var
                ;
                %let i = %eval(&i + 1);
            %end;
        quit;

        * In-memory unloading;
        %if %length(&discrete_vars) %then %do;
            sasfile &func_num._disc_data close;
        %end;
        %if %length(&continuous_vars) %then %do;
            sasfile &func_num._cont_data_T close;
        %end;

        data &func_num._WOE_by_date;
            set &func_num._woe_:;
        run;

        proc means data = &func_num._WOE_by_date nway noprint;
            class &date_var var_name;
            var IV;
            output out = &func_num._iv_by_date(drop = _:) sum =;
        run;

        data &func_num._iv;
            set &func_num._iv &func_num._iv_by_date;
        run;

        %if &graph_out = 1 %then %do;
            ods graphics on/height = 240px width = 640px;

            proc sql noprint;
                SELECT DISTINCT var_name INTO: plot_vars SEPARATED BY ' '
                    FROM &func_num._iv
                ;
            quit;

            %let i = 1;
            %do %while(%length(%scan(&plot_vars, &i)) > 0);
                %let var = %scan(&plot_vars, &i);

                data &func_num._plot_&i;
                    set &func_num._iv end = _end;

                    retain _min 0 _max 0.1;
                    _min = min(iv, _min);
                    _max = max(iv, _max);

                    if missing(&date_var) then do;
                        call symput('base_iv', iv);
                    end;
                    else output;
                    
                    if _end then do;
                        call symput('y_min', _min);
                        call symput('y_max', _max);
                    end;

                    where var_name = "&var";
                    drop _:;
                run;

                proc sgplot data = &func_num._plot_&i;
                    series x = &date_var
                           y = iv/markers lineattrs = (color = bigb)
                                          legendlabel = 'Information Value';
                    xaxis type = discrete display = (nolabel);
                    yaxis label = 'Information Value' min = &y_min max = &y_max;
                    refline &base_iv/transparency = 0.5 
                                     lineattrs = (pattern = dash)
                                     label = 'Overal IV';
                    refline 0.1/transparency = 0.5;
                    title "IV Chart for Variable = &var";
                run;
                %let i = %eval(&i + 1);
            %end;
        %end;
    %end;

    data &out_table;
        set &func_num._iv;
    run;

    %if &delete_flag = 1 %then %do;
        proc datasets lib = work nolist;
            delete &func_num.:;
        quit;
    %end;

    title;
    ods graphics on/height = 480px width = 640px; /* Default size */
%mend gp_info_value;

%macro gp_stability_chart(
                    data_table = /* Source table name                                       */,
                    var_list   = /* Variable list for stability chart. SPACE separated      */, 
                    date_var   = /* A variable by which the each var will be group          */,
                    data_type  = /* Data type of var_list. Char or Num (missing)            */

                           /* Purpose: print raw stability chart for given variables        */
                           /* Example: gp_stability_chart(data_table = sample,              */
                           /*                             var_list   = 20                   */
                           /*                             date_var   = QTR                  */
                           /*          )                                                    */
);

    %local func_num;

    * Prefix table names for better clarification;
    %let func_num = gp_02;

    proc datasets lib = work nolist;
        delete &func_num.:;
    quit;

    %if &data_type ne char %then %do;
        * Bin numeric data and generate format;
        %_gp_num_fmt(&func_num, &data_table, &var_list)

        proc format cntlin = &func_num._fmt; run;
    %end;   

    * Create a copy of external data/view - the PC way;
    proc sql noprint;
        CREATE TABLE &func_num._data AS
            SELECT *
            FROM &data_table
            ORDER BY &date_var
        ;
    quit;

    * In-memory loading;
    sasfile &func_num._data open;

    %let i = 0;
    %do %while (%length(%scan(&var_list, %eval(&i + 1))) > 0);
        %let i = %eval(&i + 1);
        %let var = %scan(&var_list, &i);
        %put Current variable in calculation: &var;

        proc freq data = &func_num._data noprint;
            by &date_var;
            table &var./out = &func_num._freq missing;
            %if &data_type ne char %then %do;
                format &var &var._.;
            %end;    
        run;
        
        proc sgplot data = &func_num._freq;
            vbar &date_var./response =  percent group =  &var groupdisplay = stack missing;
            xaxis discreteorder = data display = (nolabel);
            yaxis grid values = (0 to 100 by 10);
            title "Raw Stability Chart for Variable = &var";
        run;
    %end;

    sasfile &func_num._data close;
    title;
%mend gp_stability_chart;

%macro gp_stability_index(
                    data_table = /* Source table name                                       */,
                    out_table  = /* Output table name. Ignored if missing                   */,
                    var_list   = /* Variable list for stability index. SPACE separated      */,
                    date_var   = /* A variable by which the each var will be group          */,
                    graph_out  = /* Graph output. 1 for on and missing or 0 for off         */,
                    norm_var   = /* A variable by which the distribution is normalized.     */
                                 /* e.g. what if norm_var kept the same over the period     */
                                 /* NOTE: norm_var has to be in the var_list                */,
                    compare    = /* Compare raw index/chart to normalized ones. 1 for on    */
                                 /* and missing or 0 for off                                */,
                    data_type  = /* Data type of var_list. Char or Num (missing)            */
  
                           /* Purpose: calculate stability index for given variables. The   */
                           /*          first time period is assumed to be the expected      */
                           /*          distribution                                         */
                           /* Example: gp_stability_index(data_table = sample,              */
                           /*                             var_list   = 20                   */
                           /*                             date_var   = QTR                  */
                           /*          )                                                    */
);

    %local func_num;

    * Initialization;
    %if %length(&out_table) %then %do;
        data &out_table; set _null_; run;
    %end;

    %if %length(&norm_var) = 0 %then %do;
        %let compare = 1;
    %end;

    %let data_type = %lowcase(&data_type);
    %let date_var  = %lowcase(&date_var);
    %let norm_var  = %lowcase(&norm_var);

    * Prefix table names for better clarification;
    %let func_num = gp_03;

    proc datasets lib = work nolist;
        delete &func_num.:;
    quit;

    %if &data_type ne char %then %do;
        * Bin numeric data and generate format;
        %_gp_num_fmt(&func_num, &data_table, &var_list)

        proc format cntlin = &func_num._fmt; run;
    %end;  

    * Create a copy of external data/view - the PC way;
    proc sql noprint;
        CREATE TABLE &func_num._data AS
            SELECT *
            FROM &data_table
            ORDER BY &date_var
        ;
    quit;

    * In-memory loading;
    sasfile &func_num._data open;

    %let i = 0;
    %do %while (%length(%scan(&var_list, %eval(&i + 1))) > 0);
        %let i = %eval(&i + 1);
        %let var = %lowcase(%scan(&var_list, &i));
        %put Current variable in calculation: &var;

        * Calculate raw PCTs;
        proc freq data = &func_num._data noprint;
            by &date_var;
            table &var./out = &func_num._freq missing;
            %if &data_type ne char %then %do;
                format &var &var._.;
            %end;
        run;

        * Output the expected PCT;
        %if &data_type ne char %then %do;
            data &func_num._freq &func_num._expected(rename = (percent = pct_exp));
                set &func_num._freq(rename = (&var = _drop));
                by &date_var;
                length &var. $36;
                &var = put(_drop, &var._.);

                retain _first;
                if _N_ = 1 then _first = &date_var;
                if &date_var = _first then output &func_num._expected;
                output &func_num._freq;

                drop _: count;
            run;
        %end;
        %else %do;
            data &func_num._freq &func_num._expected(rename = (percent = pct_exp));
                set &func_num._freq;
                by &date_var;

                retain _first;
                if _N_ = 1 then _first = &date_var;
                if &date_var = _first then output &func_num._expected;
                output &func_num._freq;

                drop _: count;
            run;
        %end;

        %if &compare = 1 %then %do; 
            * Calculate the raw SI;
            data &func_num._freq;
                if _N_ = 1 then do;
                    if 0 then set &func_num._expected;
                    declare hash map(dataset: "&func_num._expected");
                    map.definekey("&var");
                    map.definedata('pct_exp');
                    map.definedone();
                end;

                set &func_num._freq;

                rc = map.find();
                if rc ne 0 then pct_exp = 0;

                percent = max(percent, 1e-3);
                pct_exp = max(pct_exp, 1e-3);

                stab_idx = (percent - pct_exp) * log(percent / pct_exp) / 100;
            run;

            proc means data = &func_num._freq nway noprint;
                class &date_var;
                var stab_idx;
                output out = &func_num._si(drop = _:) sum =;
            run;

            data &func_num._si;
                set &func_num._si;

                stab_idx_diff = sum(stab_idx - lag1(stab_idx), 0);
            run;

            %if %length(&out_table) %then %do;
                data &func_num._si;
                    set &func_num._si;
                    length type $10 var_name $36;

                    type = 'raw';
                    var_name = "&var";
                run;

                data &out_table;
                    set &out_table &func_num._si;
                run;
            %end;
        %end;

        * Normalization;
        %if %length(&norm_var) and (&norm_var ne &var) %then %do;
            proc freq data = &func_num._data noprint;
                by &date_var;
                table &var. * &norm_var./outpct out = &func_num._freq_cross_raw missing;
                %if &data_type ne char %then %do;
                    format &var &var._.
                           &norm_var &norm_var._.;
                %end;
            run;
    
            %if &data_type ne char %then %do;
                data &func_num._freq_cross_raw 
                     &func_num._freq_cross_base(rename = (
                                                    pct_row = pct_row_exp 
                                                    percent = pct_exp
                                               ));
                    set &func_num._freq_cross_raw(rename = (
                                                    &var = _drop1
                                                    &norm_var = _drop2
                                                 ));
                    by &date_var;
                    length &var &norm_var $36;
                    &var = put(_drop1, &var._.);
                    &norm_var = put(_drop2, &norm_var._.);

                    retain _first;
                    if _N_ = 1 then _first = &date_var;

                    if &date_var = _first then output &func_num._freq_cross_base;
                    output &func_num._freq_cross_raw;

                    drop _: pct_col count;
                run;
            %end;
            %else %do;
                data &func_num._freq_cross_raw 
                     &func_num._freq_cross_base(rename = (
                                                    pct_row = pct_row_exp
                                                    percent = pct_exp
                                               ));
                    set &func_num._freq_cross_raw;
                    by &date_var;

                    retain _first;
                    if _N_ = 1 then _first = &date_var;

                    if &date_var = _first then output &func_num._freq_cross_base;
                    output &func_num._freq_cross_raw;

                    drop _: pct_col count;
                run;
            %end;
            
            * Calculated normalized percent;
            data &func_num._freq_cross_raw;
                if _N_ = 1 then do;
                    if 0 then set &func_num._freq_cross_base;
                    declare hash map(dataset: "&func_num._freq_cross_base");
                    map.definekey("&var", "&norm_var");
                    map.definedata('pct_row_exp', 'pct_exp');
                    map.definedone();
                end;
                set &func_num._freq_cross_raw;

                rc = map.find();
                if rc ne 0 then do;
                    pct_exp = 0;
                    pct_row_exp = pct_row; /* if a new category, then no adjustment */
                end;
                
                * Preventing 0/0, num/0 cases and make 0/0 = 1;
                percent = max(percent, 1e-3);
                pct_exp = max(pct_exp, 1e-3);
                pct_row = max(pct_row, 1e-3);
                pct_row_exp = max(pct_row_exp, 1e-3);

                percent = percent * pct_row_exp / pct_row;
            run;

            proc means data = &func_num._freq_cross_raw nway noprint;
                class &date_var;
                var percent;
                output out = &func_num._freq_cross_sum(drop = _:) sum = pct_sum;
            run;

            * Calculate normalized SIs;
            data &func_num._freq_cross;
                if _N_ = 1 then do;
                    if 0 then set &func_num._freq_cross_sum;
                    declare hash map(dataset: "&func_num._freq_cross_sum");
                    map.definekey("&date_var");
                    map.definedata('pct_sum');
                    map.definedone();
                end;
                set &func_num._freq_cross_raw;
                
                by &date_var;
                
                rc = map.find();
                if rc ne 0 then put 'gp_stability_index: How is this possible';
                percent = percent / pct_sum * 100; /* Enforce sum(percent) = 100 */
                
                stab_idx = (percent - pct_exp) * log(percent / pct_exp) / 100;
            run;
            
            proc means data = &func_num._freq_cross nway noprint;
                class &date_var;
                var stab_idx;
                output out = &func_num._si_norm(drop = _:) sum =;
            run;

            data &func_num._si_norm;
                set &func_num._si_norm;

                stab_idx_diff = sum(stab_idx - lag1(stab_idx), 0);
            run;

            %if %length(&out_table) %then %do;
                data &func_num._si_norm;
                    set &func_num._si_norm;
                    length type $10 var_name $36;

                    type = 'normalized';
                    var_name = "&var";
                run;

                data &out_table;
                    set &out_table &func_num._si_norm;
                run;
            %end;
        %end;
            
        %if &graph_out = 1 %then %do;
            * Chart output layout depending on if there is a comparison;
            %if not (&compare = 1) or %length(&norm_var) = 0 or &var = &norm_var %then %do;
                %if &compare = 1 or &var = &norm_var %then %do; /* Raw SIs */
                    %let table_1 = &func_num._freq;
                    %let table_2 = &func_num._si;
                    %let title_1 = Raw Stability Chart for var = &var;
                    %let title_2 = Raw Stability Index for var = &var;
                %end;
                %else %if &var ne &norm_var %then %do; /* Normalized SIs */
                    %let table_1 = &func_num._freq_cross;
                    %let table_2 = &func_num._si_norm;
                    %let title_1 = Normalized Stability Chart for var = &var./norm_var = &norm_var;
                    %let title_2 = Normalized Stability Index for var = &var./norm_var = &norm_var;
                %end;
                ods graphics on/height = 480px width = 640px;
                proc sgplot data = &table_1;
                    vbar &date_var./response = percent group = &var groupdisplay = stack missing;
                    xaxis discreteorder = data display = (nolabel);
                    yaxis grid values = (0 to 100 by 10);
                    title "&title_1";
                run;

                ods graphics on/height = 240px width = 640px;
                proc sgplot data = &table_2;
                    series x = &date_var
                           y = stab_idx/markers lineattrs = (color = bigb)
                                                legendlabel = 'Stability Index';
                    series x = &date_var
                           y = stab_idx_diff/ lineattrs = (pattern = dash color = bigb)
                                              legendlabel = 'Stability Index Delta'
                                              transparency = 0.4;
                    xaxis type = discrete display = (nolabel);
                    yaxis max = 1 label = 'Stability Index';
                    refline 0.2/transparency = 0.6 lineattrs = (pattern = dash) label = 'Alert';
                    refline 0.4/transparency = 0.3 lineattrs = (pattern = dash) label = 'Warning';
                    title "&title_2";
                run;
            %end;
            %else %if &var ne &norm_var %then %do;
                data &func_num._freq_total;
                    set &func_num._freq_cross(in = a)
                        &func_num._freq(in = b);
                    length type $10;

                    if a then type = 'Normalized';
                    if b then type = 'Raw';
                run;

                data &func_num._si_total;
                    set &func_num._si(rename = (
                                        stab_idx = stab_idx_raw
                                        stab_idx_diff = stab_idx_raw_diff
                                     ));
                    set &func_num._si_norm(rename = (
                                            stab_idx = stab_idx_norm
                                            stab_idx_diff = stab_idx_norm_diff
                                          ));
                run;

                ods graphics on/height = 400px width = 800px;
                proc sgpanel data = &func_num._freq_total;
                    panelby type;
                    vbar &date_var./response = percent group = &var groupdisplay = stack missing;
                    colaxis display = (nolabel);
                    title "Stability Chart Comparison for var = &var";
                run;
                    
                ods graphics on/height = 200px width = 800px;
                proc sgplot data = &func_num._si_total;
                    series x = &date_var
                           y = stab_idx_raw/markers lineattrs = (color = bigb)
                                                    legendlabel = 'Raw Stability Index';
                    series x = &date_var
                           y = stab_idx_norm/markers lineattrs = (color = crimson)
                                                     legendlabel = 'Normalized Stability Index';
                    series x = &date_var
                           y = stab_idx_raw_diff/ lineattrs = (pattern = dash color = bigb) 
                                                  legendlabel = 'Raw Stability Index Delta'
                                                  transparency = 0.4;
                    series x = &date_var
                           y = stab_idx_norm_diff/ lineattrs = (pattern = dash color = crimson) 
                                                   legendlabel = 'Normalized Stability Index Delta'
                                                   transparency = 0.4;
                    xaxis type = discrete display = (nolabel);
                    yaxis max = 1 label = 'Stability Index';
                    refline 0.2/transparency = 0.6 lineattrs = (pattern = dash) label = 'Alert';
                    refline 0.4/transparency = 0.3 lineattrs = (pattern = dash) label = 'Warning';
                    title "Stability Index Comparison for var = &var./norm_var = &norm_var";
                run;
            %end;
        %end;    
    %end;

    * In-memory unloading;
    sasfile &func_num._data close;
    title;
    ods graphics on/height = 480px width = 640px; /* Default size */
%mend gp_stability_index;


%macro gp_var_split(
                    out_table = /* Output table name.                                       */,
                    var_list  = /* Variable list for stability chart. SPACE separated       */,
                    blk_size  = /* Number of variables for each block                       */
  
                           /* Purpose: split var_list into smaller blocks whose size is     */
                           /*          specified in the blk_size argument                   */
                           /* Example: gp_var_split(out_table = my_blk,                     */
                           /*                       var_list  = var_1 var_2 var_3,          */
                           /*                       blk_size  = 300                         */
                           /*          )                                                    */
);

    %local func_num;

    * Prefix table names for better clarification;
    %let func_num = gp_04;

    proc datasets lib = work nolist;
        delete &func_num.:;
    quit;
    
    %if %length(&out_table) %then %do;
        * Cannot use cards in macro so no input @@ trick;
        /*
            data &out_table;
                input @@;
                _infile_ = "&var_list";
                input var_name $ @@;
                cards; 
            run;
        */
        
        %let i = 0;
        data &out_table;
            length var_name $36;
            _count = 0;
            block  = 1;

            %do %while (%length(%scan(&var_list, %eval(&i + 1))));
                %let i = %eval(&i + 1);
                %let var = %scan(&var_list, &i);
                var_name = "&var";
                _count = _count + 1;
                if _count > &blk_size then do;
                    _count = 1;
                    block  = block + 1;
                end;
                output;
            %end;
            drop _:;
        run;
    %end;
%mend gp_var_split;

%macro gp_univariate_c(
                    data_table  = /* Source table name                                      */,
                    out_table   = /* Output table name. Ignored if missing                  */,
                    perf_var    = /* Performance variable such AS ODEFIND                   */,
                    delete_flag = /* 1 or 0(miss). Delete all temporary data after run      */,
                    trans_data  = /* Dataset containing monotonically transformed data by   */
                                  /* proc transreg.                                         */

                        /* Purpose: calculate concordance value for each variable based on  */
                        /*          univariate logistic regression                          */
                        /* Example: gp_univariate_c(data_table = sample,                    */
                        /*                          out_table  = example,                   */
                        /*                          perf_var   = ODEFIND                    */
                        /*          )                                                       */
);

    %local func_num;

    * Prefix table names for better clarification;
    %let func_num = gp_05;

    proc datasets lib = work nolist;
        delete &func_num.:;
    quit;

    * Initialization;
    %let num_vars  =;
    %let char_vars =;
    %let perf_var  = %lowcase(&perf_var);

    data &func_num._concord;
        set _null_;
    run;
    
    * Obtain variable names and types;
    ods select variables;
    ods output variables = &func_num._var_type;

    proc contents data = &data_table; run;

    ods select all;

    proc sql noprint;
        SELECT variable INTO: num_vars SEPARATED BY ' '
            FROM &func_num._var_type
            WHERElowcase(variable) ne "&perf_var" and lowcase(type) = 'num'
        ;

        SELECT variable INTO: char_vars SEPARATED BY ' '
            FROM &func_num._var_type
            WHERElowcase(variable) ne "&perf_var" and lowcase(type) = 'char'
        ;
    quit;

    data &func_num._char_data(keep = &char_vars &perf_var)
         &func_num._num_data(keep = &num_vars &perf_var);
        set &data_table;
    run;

    * Monotonic transformation;
    %if %length(&num_vars) > 0 %then %do;
        %if %length(&trans_data) %then %do;
            %let var_keep = T%sysfunc(tranwrd(&num_vars, %str( ), %str( T)));
            data &func_num._num_data_T;
                set &trans_data;
                keep &perf_var &var_keep;
            run;
        %end;
        %else %do;
            proc transreg data = &func_num._num_data noprint;
                model identity(&perf_var) = monotone(&num_vars);
                output out = &func_num._num_data_T;
            run;
        %end;
    %end;

    * In-memory loading;
    %if %length(&char_vars) %then %do;
        sasfile &func_num._char_data open;
    %end;

    %let i = 1;
    %do %while(%length(%scan(&char_vars, &i)));
        %let var  = %scan(&char_vars, &i);
        
        * Calculate the condordance value;
        ods select Association; 
        ods output Association = &func_num._association;
        

        proc logistic data = &func_num._char_data;
            class &var;
            model &perf_var = &var;
        run;

        ods select all;

        data &func_num._association;
            set &func_num._association;
            length var_name $36;
            var_name = "&var";
            keep nValue2 var_name;
            rename nValue2 = concord;
            where Label2 = 'c';
        run;

        data &func_num._concord;
            set &func_num._concord &func_num._association;
        run;
        %let i = %eval(&i + 1);
    %end;

    * In-memory unloading;
    %if %length(&char_vars) %then %do;
        sasfile &func_num._char_data close;
    %end;

    * In-memory loading;
    %if %length(&num_vars) %then %do;
        sasfile &func_num._num_data_T open;
    %end;

    %let i = 1;
    %do %while(%length(%scan(&num_vars, &i)));
        %let var  = %scan(&num_vars, &i);

        ods select Association; 
        ods output Association = &func_num._association;
        
        proc logistic data = &func_num._num_data_T;
            model &perf_var = T&var;
        run;

        ods select all;

        data &func_num._association;
            set &func_num._association;
            length var_name $36;
            var_name = "&var";
            keep nValue2 var_name;
            rename nValue2 = concord;
            where Label2 = 'c';
        run;

        data &func_num._concord;
            set &func_num._concord &func_num._association;
        run;
        %let i = %eval(&i + 1);
    %end;

    * In-memory unloading;
    %if %length(&num_vars) %then %do;
        sasfile &func_num._num_data_T close;
    %end;

    %if %length(&out_table) %then %do;
        data &out_table;
            set &func_num._concord;
        run;
    %end;

    %if &delete_flag = 1 %then %do;
        proc datasets lib = work nolist;
            delete &func_num.:;
        quit;
    %end;
%mend gp_univariate_c;

%macro gp_fast_transreg(
                    data_table  = /* Source table name                                      */,
                    out_table   = /* Output table name. Ignored if missing                  */,
                    perf_var    = /* Performance variable such AS ODEFIND                   */,
                    var_list    = /* List of interested numeric variables. SPACE separated  */, 
                    delete_flag = /* 1 or 0(miss). Delete all temporary data after run      */,
                    blk_size    = 100
                                  /* block size for each subset of variables                */    
                      
                        /* Purpose: accelarate transreg speed by splitting calculationg     */
                        /*          into smaller pieces                                     */
                        /* Example: gp_fast_transreg(data_table = sample,                   */
                        /*                           out_table   = example,                 */
                        /*                           perf_var    = ODEFIND                  */
                        /*                           var_list    = A B C D                  */
                        /*          )                                                       */
);

    %local func_num;

    * Prefix table names for better clarification;
    %let func_num = gp_06;

    proc datasets lib = work nolist;
        delete &func_num.:;
    quit;

    * Initialization;
    %let perf_var  = %lowcase(&perf_var);

    %gp_var_split(out_table = &func_num._var_blk, 
                  var_list  = &var_list, 
                  blk_size  = &blk_size)
    
    data &func_num._data;
        set &data_table;
        keep &var_list &perf_var;
    run;

    * In memory loading;
    sasfile &func_num._data open;

    %let i = 1;
    proc sql noprint;
        SELECT var_name INTO: curr_vars SEPARATED BY ' '
            FROM &func_num._var_blk
            WHERE block = &i
        ;
    quit;

    %do %while(%length(&curr_vars));
        proc transreg data = &func_num._data(keep = &curr_vars &perf_var) noprint;
            model identity(&perf_var) = monotone(&curr_vars);
            output out = &func_num._data_&i;
        run;

        %let i = %eval(&i + 1);
        %let curr_vars =;
        proc sql noprint;
            SELECT var_name INTO: curr_vars SEPARATED BY ' '
                FROM &func_num._var_blk
                WHERE block = &i
            ;
        quit;
    %end;

    %if %length(&out_table) %then %do;
        data &out_table;
            merge &func_num._data_:;
        run;
    %end;

    %if &delete_flag = 1 %then %do;
        proc datasets lib = work nolist;
            delete &func_num.:;
        quit;
    %end;

    sasfile &func_num._data close;

%mend gp_fast_transreg;

%macro gp_cluster_corr(
                    data_table  = /* Source table name, should be the OUTSTATS FROM Varclus */,
                    out_table   = /* Output table name. Ignored if missing                  */,
                    agg_funcs   = mean
                                  /* List of aggregating function(s). Default is mean       */,
                    delete_flag = /* 1 or 0(miss). Delete all temporary data after run      */

                        /* Purpose: calculate statistics of clusters correlation by each    */
                        /*          step of decomposition                                   */
                        /* Example: gp_cluster_corr(data_table = sample,                    */
                        /*                          out_table  = example,                   */
                        /*                          agg_func   = p75                        */
                        /*          )                                                       */
);

    %local func_num;

    * Prefix table names for better clarification;
    %let func_num = gp_07;

    proc datasets lib = work nolist;
        delete &func_num.:;
    quit;

    ods select Position;
    ods output Position = &func_num._position;
    proc contents data = &data_table varnum; run;
    ods select all;

    * Supposedly, the first 3 variables are _NCL_ _TYPE_ _NAME_;
    data _null_;
        set &func_num._position end = _end;
        
        if Num = 4 then do;
            call symput('start_var', Variable);
        end;

        if _end then do;
            call symput('end_var', Variable);
            call symput('end_num', Num);
        end;
    run;

    %let dim = %eval(&end_num - 3);
    
    data &func_num._ccorr;
        set &data_table;
        where _TYPE_ = 'CCORR';
    run;

    data &func_num._ccorr;
        length clus_num 4;
        set &func_num._ccorr;
        clus_num = input(substr(_NAME_, 5), best4.);

        array old_clus{*} &start_var -- &end_var;
        array new_clus{*} clus_1 - clus_&dim;

        do i = 1 to dim(old_clus);
            if i < clus_num then do;
                new_clus{i} = old_clus{i};
            end;
        end;

        drop WOE_: i _NAME_ _TYPE_;
        rename _NCL_ = num_clusters;
    run;

    proc transpose data = &func_num._ccorr out = &func_num._ccorr_t;
        by num_clusters clus_num;
    run;

    data &func_num._ccorr_t;
        set &func_num._ccorr_t;
        COL1 = abs(COL1);
        where COL1 ne .;
    run;

    proc means data = &func_num._ccorr_t nway noprint;
        class num_clusters;
        var COL1;
        output out = &func_num._out(drop = _:) 
            %let i = 1;
            %do %while (%length(%scan(&agg_funcs, &i)));
                %let agg_func = %scan(&agg_funcs, &i);
                &agg_func = &agg_func._cluster_corr
                %let i = %eval(&i + 1);
            %end;
            ;
    run;

    %if %length(&out_table) %then %do;
        data &out_table;
            set &func_num._out;
        run;
    %end;

    %if &delete_flag = 1 %then %do;
        proc datasets lib = work nolist;
            delete &func_num.:;
        quit;
    %end;

%mend gp_cluster_corr;

%macro gp_cluster_size(
                    data_table  = /* Source table name, should be the OUTSTATS FROM Varclus */,
                    out_table   = /* Output table of KDE. Ignored if missing                */,
                    log_flag    = /* 1 or 0(miss). If log transform the size of clusters    */,
                    plot_flag   = /* 1 or 0(miss). If plot the histogram and density        */,
                    delete_flag = /* 1 or 0(miss). Delete all temporary data after run      */,
                    plot_clus_list =
                                  /* A list of clusters represented by the number to be     */
                                  /* plotted. If not given, a plot of 21 subplots based on  */
                                  /* log interval will be generated                         */,
                    plot_hist_xmax = 20
                                  /* Max value of xaxis for the histogram plot              */,
                    plot_num_cols = 7
                                  /* Number of columns in the plots                         */,
                    plot_width    = 6in
                                  /* Plot width passed to ods graphics width option         */,
                    plot_height   = 8in
                                  /* Plot height passed to ods graphics height option       */

                        /* Purpose: calculate cluster size distribution by number of        */
                        /*          clusters. Histogram and KDE plots are optional          */
                        /* Example: gp_cluster_size(data_table = sample,                    */
                        /*                          out_table  = kde_out,                   */
                        /*                          plot_flag  = 1,                         */
                        /*                          plot_clus_list = 1 3 5 8 10 11 125      */
                        /*          )                                                       */
);

    %local func_num;

    * Prefix table names for better clarification;
    %let func_num = gp_08;

    * Initialization of variables;
    %if &log_flag = 1 %then %let var_name = log_num;
    %else %let var_name = COL1;

    %if %length(&plot_clus_list) %then %do;
        %let plot_flag = 1;
    %end;
        
    proc datasets lib = work nolist;
        delete &func_num.:;
    quit;

    * Obtain member information;
    data &func_num._members;
        set &data_table;
        where _TYPE_ = 'MEMBERS';
        rename _NCL_ = num_clusters;
        drop _NAME_ _TYPE_; 
    run;

    proc transpose data = &func_num._members out = &func_num._members_t(drop = _:);
        by num_clusters;
    run;

    data &func_num._members_t;
        set &func_num._members_t end = _end;
        log_num = log(COL1);
        where COL1 ne .;

        if _end then call symput('num_clus', num_clusters);
        drop _:;
    run;

    %if &log_flag = 1 %then %let gridu = 6;
    %else %let gridu = &num_clus; /* Just an approximation of ymax */

    * Density estimation;
    proc kde data = &func_num._members_t;
        by num_clusters;
        univar &var_name (gridl=1 gridu=&gridu) / ngrid = 21 plots = none
        out = &func_num._den noprint;
        where num_clusters > 1; /* Avoid error of insufficient data */
    run;
    
    %if %length(&out_table) %then %do;
        data &out_table;
            set &func_num._den;
        run;
    %end;

    * Plot based on specification;
    %if &plot_flag = 1 %then %do;
        data &func_num._den_plot;
            set &func_num._den;
            by num_clusters;

            low = 0;

            * Clean density data;
            if density < 1e-3 then density = 0;

            %if %length(&plot_clus_list) %then %do;
                where num_clusters in (&plot_clus_list);
            %end;
            %else %do;
                * Make a 3 x 7 plot;
                retain i is_plot 0;
                interval = log(&num_clus) / 22;
                if first.num_clusters then is_plot = 0;
                if log(num_clusters) >= i * interval then do;
                    i = i + 1;
                    is_plot = 1;
                end; 
            
                drop i is_plot interval;
                if is_plot;
            %end;
        run;

    ods graphics / reset width = &plot_width height = &plot_height;
    title 'Histogram Plot of Number of Members for Each Cluster by Total Cluster Number';

    proc sgpanel data = &func_num._den_plot nocycleattrs;
        panelby num_clusters / layout = panel onepanel novarname noborder
                               colheaderpos = bottom columns = &plot_num_cols;
        highlow y = value high = count low = low / lineattrs = (thickness = 4) 
                                                   transparency = 0.5;
        rowaxis label = 'Log(Num of Members in Cluster)' grid;
        colaxis display = (nolabel) valueattrs = (size = 6) max = &plot_hist_xmax;
    run;
    
    /* 
        There is a bug where if one has ods graphcis / reset at the beginning and before every
        following plot, all plots will just be the last plot.
    */

    title 'Density Plot of Number of Members for Each Cluster by Total Cluster Number';

    proc sgpanel data = &func_num._den_plot nocycleattrs;
        panelby num_clusters / layout = panel onepanel novarname noborder
                               colheaderpos = bottom columns = &plot_num_cols;
        highlow y = value high = density low = low / lineattrs = (thickness = 4) 
                                                     transparency = 0.5;
        rowaxis label = 'Log(Num of Members in Cluster)' grid;
        colaxis display = (nolabel) valueattrs = (size = 6);
    run;

    title '';
    ods graphics / reset;    

    %end;

    %if &delete_flag = 1 %then %do;
        proc datasets lib = work nolist;
            delete &func_num.:;
        quit;
    %end;

%mend gp_cluster_size;

%macro gp_cluster_struct(
                    data_table  = /* Source table name, should be the OUTSTATS FROM Varclus */,
                    out_table   = /* Output table of KDE. Ignored if missing                */,
                    select_clus = /* Number representing which cluster is selected          */,
                    delete_flag = /* 1 or 0(miss). Delete all temporary data after run      */

                        /* Purpose: Show the cluster structure and RSquared Ratios          */

);

    %local func_num;

    * Prefix table names for better clarification;
    %let func_num = gp_09;

    proc datasets lib = work nolist;
        delete &func_num.:;
    quit;

    * Pad the number with 0 so it could be sorted correctly;
    data &func_num._structure;
        set &data_table;

        if _NAME_ = 'Clus1' then _NAME_ = 'Clus01';
        if _NAME_ = 'Clus2' then _NAME_ = 'Clus02';
        if _NAME_ = 'Clus3' then _NAME_ = 'Clus03';
        if _NAME_ = 'Clus4' then _NAME_ = 'Clus04';
        if _NAME_ = 'Clus5' then _NAME_ = 'Clus05';
        if _NAME_ = 'Clus6' then _NAME_ = 'Clus06';
        if _NAME_ = 'Clus7' then _NAME_ = 'Clus07';
        if _NAME_ = 'Clus8' then _NAME_ = 'Clus08';
        if _NAME_ = 'Clus9' then _NAME_ = 'Clus09';

        rename _NAME_ = cluster;

        where _NCL_ = &select_clus and _TYPE_ = 'STRUCTUR';
    run;

    proc transpose data = &func_num._structure out = &func_num._structure_t;
        by cluster;
        var WOE_:;
    run;
    
    * Calculate rsquared;
    data &func_num._structure_t;
        set &func_num._structure_t;
        COL1 = COL1 ** 2;
        rename _NAME_ = var_name
               COL1 = rsquared
        ;
        drop _LABEL_;
    run;

    proc sort data = &func_num._structure_t;
        by var_name rsquared;
    run;

    * Output rsquare ratios;
    data &func_num._out;
        set &func_num._structure_t;
        by var_name rsquared;

        retain rsquared_next clus_next;
        
        * Initialize to missing;
        if first.var_name then rsquared_next = .;
        
        if last.var_name then do;
            rsquared_ratio = (1 - rsquared) / (1 - rsquared_next);
            output;
        end;

        rsquared_next = rsquared;
        clus_next = cluster;
    run;

    proc sort data = &func_num._out;
        by cluster rsquared_ratio;
    run;

    %if %length(&out_table) %then %do;
        data &out_table;
            set &func_num._out;
        run;
    %end;

    %if &delete_flag = 1 %then %do;
        proc datasets lib = work nolist;
            delete &func_num.:;
        quit;
    %end;

%mend gp_cluster_struct;

%macro gp_cluster_add_vars(
                    data_table  = /* Data source table                                      */,
                    score_table = /* OUTSTATS FROM Varclus                                  */,
                    out_table   = /* Output table of KDE. Ignored if missing                */,
                    select_clus = /* Number representing which cluster is selected          */,
                    select_vars = /* Selected variables to add to the existing structure. If*/
                                  /* missing, all variables in data_table will be added     */,
                    delete_flag = /* 1 or 0(miss). Delete all temporary data after run      */

                        /* Purpose: Show the cluster structure and RSquared Ratios          */

);

    %local func_num;

    * Prefix table names for better clarification;
    %let func_num = gp_10;

    proc datasets lib = work nolist;
        delete &func_num.:;
    quit;

    * Re-construct the cluster components;
    %_gp_clus_score(&func_num, &data_table, &score_table, &func_num._cmp, &select_clus)

    data &func_num._cmp;
        set &func_num._cmp;
        set &data_table.(keep = &select_vars);
    run;

    proc corr data = &func_num._cmp outp = &func_num._corr noprint;
        var woe:;
        with clus:;
    run;

    data &func_num._corr;
        set &func_num._corr;

        if _NAME_ = 'Clus1' then _NAME_ = 'Clus01';
        if _NAME_ = 'Clus2' then _NAME_ = 'Clus02';
        if _NAME_ = 'Clus3' then _NAME_ = 'Clus03';
        if _NAME_ = 'Clus4' then _NAME_ = 'Clus04';
        if _NAME_ = 'Clus5' then _NAME_ = 'Clus05';
        if _NAME_ = 'Clus6' then _NAME_ = 'Clus06';
        if _NAME_ = 'Clus7' then _NAME_ = 'Clus07';
        if _NAME_ = 'Clus8' then _NAME_ = 'Clus08';
        if _NAME_ = 'Clus9' then _NAME_ = 'Clus09';

        rename _NAME_ = cluster;
        where _NAME_ is not missing;
        drop _TYPE_;
    run;

    proc transpose data = &func_num._corr out = &func_num._corr_t;
        by cluster;
    run;

    data &func_num._corr_t;
        set &func_num._corr_t;

        COL1 = COL1 ** 2;

        rename _NAME_ = var_name
               COL1 = rsquared;
        drop _LABEL_;
    run;

    proc sort data = &func_num._corr_t;
        by var_name rsquared;
    run;

    * Output rsquare ratios;
    data &func_num._out;
        set &func_num._corr_t;
        by var_name rsquared;

        retain rsquared_next clus_next;
        
        * Initialize to missing;
        if first.var_name then rsquared_next = .;
        
        if last.var_name then do;
            rsquared_ratio = (1 - rsquared) / (1 - rsquared_next);
            output;
        end;

        rsquared_next = rsquared;
        clus_next = cluster;
    run;

    proc sort data = &func_num._out;
        by cluster rsquared_ratio;
    run;

    %if %length(&out_table) %then %do;
        data &out_table;
            set &func_num._out;
        run;
    %end;

    %if &delete_flag = 1 %then %do;
        proc datasets lib = work nolist;
            delete &func_num.:;
        quit;
    %end;

%mend gp_cluster_add_vars;

%macro gp_template(
                    data_table  = /* Source table name                                      */,
                    out_table   = /* Output table name. Ignored if missing                  */,
                    delete_flag = /* 1 or 0(miss). Delete all temporary data after run      */

                        /* Purpose: a dummy macro AS a coding style reference               */

);

    %local func_num;

    * Prefix table names for better clarification;
    %let func_num = gp_07;

    proc datasets lib = work nolist;
        delete &func_num.:;
    quit;

    /*
        Some code here
    */

    %if %length(&out_table) %then %do;
        data &out_table;
            set &func_num._out;
        run;
    %end;

    %if &delete_flag = 1 %then %do;
        proc datasets lib = work nolist;
            delete &func_num.:;
        quit;
    %end;

%mend gp_template;

*********************************************************************************************;
***                                       Help Macro                                      ***;    
*********************************************************************************************;

%macro gp_help(m_name);
%put ****************************************************************************************;
    %if %lowcase(%trim(&m_name)) = gp_info_value %then %do;
        %put macro: gp_info_value;
        %put ;
        %put purpose: calculates IV value based on WOE method.;
        %put ;
        %put arguments:;
        %put        data_table  = : Source table name;                                        
        %put        out_table   = : Output table name;                                        
        %put        n_threshold = : Nlevel threshold to determine a variable if continuous;   
        %put        perf_var    = : Performance variable such as ODEFIND;                     
        %put        date_var    = : If given, also calculate IV over time by date_var;        
        %put        graph_out   = : Graph output if date_var is given. 1 for on and;        
        %put                        missing or 0 for off;                                     
        %put        delete_flag     1 or 0(miss). Delete all temporary data after running;    
        %put        special_woe = : Special value of WOE when sum(good) or sum(bad) = 0;
        %put ;
        %put example:;
        %put          gp_info_value(data_table  = sample,; 
        %put                        out_table   = example,;         
        %put                        n_threshold = 20;     
        %put                        perf_var    = ODEFIND; 
        %put          );                    
    %end;
    %else %if %lowcase(%trim(&m_name)) = gp_stability_chart %then %do;
        %put macro: gp_stability_chart;
        %put ;
        %put purpose: print raw stability chart for given variables;
        %put ;
        %put arguments:;
        %put        data_table = : Source table name;                                           
        %put        var_list   = : Variable list for stability chart. SPACE separated;    
        %put        date_var   = : A variable by which the each var will be group; 
        %put        data_type  = : Data type of var_list. Char or Num (missing);
        %put ;
        %put example:;
        %put          gp_stability_chart(data_table = sample,;   
        %put                             var_list   = OFFR_DT FICO_SCR_NB,;         
        %put                             date_var   = QTR;        
        %put          );                       
    %end;
    %else %if %lowcase(%trim(&m_name)) = gp_stability_index %then %do;
        %put macro: gp_stability_index;
        %put ;
        %put purpose: calculate stability index for given variables. The;
        %put          first time period is assumed to be the expected distribution;
        %put ;
        %put arguments:;
        %put        data_table = : Source table name;                                     
        %put        out_table  = : Output table name. If missing then ignored;            
        %put        var_list   = : Variable list for stability chart. SPACE separated;    
        %put        date_var   = : A variable by which the each var will be group;        
        %put        graph_out  = : Graph output. 1 for on and missing or 0 for off;       
        %put        norm_var   = : A variable by which the distribution is normalized.; 
        %put                       e.g. what if norm_var kept the same over the period;   
        %put                       NOTE: norm_var has to be in the var_list;              
        %put        compare    = : Compare raw index/chart to normalized ones. 1 for on;  
        %put                       and missing or 0 for off;   
        %put        data_type  = : Data type of var_list. Char or Num (missing); 
        %put ;
        %put example:;
        %put          gp_stability_index(data_table = sample,;   
        %put                             var_list   = OFFR_DT FICO_SCR_NB,;         
        %put                             date_var   = QTR,;
        %put                             norm_var   = FICO_SCR_NB; 
        %put          );                       
    %end;
    %else %if %lowcase(%trim(&m_name)) = gp_var_split %then %do;
        %put macro: gp_var_split;
        %put ;
        %put purpose: split var_list into smaller blocks whose size is;
        %put          specified in the blk_size argument;             
        %put ;
        %put arguments:;
        %put        out_table = : Output table name;                                    
        %put        var_list  = : Variable list for stability chart. SPACE separated;   
        %put        blk_size  = : Number of variables for each block;                  
        %put ;
        %put example:;
        %put          gp_var_split(out_table = my_blk,;                 
        %put                       var_list  = var_1 var_2 var_3,;               
        %put                       blk_size  = 300;                     
        %put          );                                   
    %end;
    %else %if %lowcase(%trim(&m_name)) = gp_univariate_c %then %do;
        %put macro: gp_univariate_c;
        %put ;
        %put purpose: calculate concordance value for each variable based on;
        %put          sunivariate logistic regression;                        
        %put ;
        %put arguments:;
        %put        data_table  = : Source table name;                                    
        %put        out_table   = : Output table name. Ignored if missing;   
        %put        perf_var    = : Performance variable such as ODEFIND;                
        %put        delete_flag = : 1 or 0(miss). Delete all temporary data after run;
        %put ;
        %put example:;
        %put          gp_univariate_c(data_table = sample,;       
        %put                          out_table  = example,;               
        %put                          perf_var   = ODEFIND;       
        %put          );                           
    %end;
    %else %if %lowcase(%trim(&m_name)) = gp_fast_transreg %then %do;
        %put macro: gp_fast_transreg;
        %put ;
        %put purpose: accelarate transreg speed by splitting calculationg;
        %put          into smaller pieces;                                
        %put ;
        %put arguments:;
        %put        data_table  = : Source table name;                                        
        %put        out_table   = : Output table name. Ignored if missing;                
        %put        perf_var    = : Performance variable such as ODEFIND;                    
        %put        var_list    = : List of interested numeric variables. SPACE separated;
        %put        delete_flag = : 1 or 0(miss). Delete all temporary data after run;    
        %put        blk_size    = : block size for each subset of variables;
        %put ;                             
        %put example:;
        %put          gp_fast_transreg(data_table = sample,;        
        %put                           out_table  = example,;               
        %put                           perf_var   = ODEFIND;       
        %put                           var_list   = A B C D; 
        %put          );                     
    %end;
    %else %if %lowcase(%trim(&m_name)) = gp_cluster_corr %then %do;
        %put macro: gp_cluster_corr;
        %put ;
        %put purpose: calculate cluster size distribution by number of; 
        %put          clusters. Histogram and KDE plots are optional;   
        %put ;
        %put arguments:;
        %put        data_table  = : Source table name, should be the OUTSTATS FROM Varclus;
        %put        out_table   = : Output table name. Ignored if missing;   
        %put        agg_funcs   = : List of aggregating functions. Default is mean;   
        %put        delete_flag = : 1 or 0(miss). Delete all temporary data after run;
        %put ;
        %put example:;
        %put          gp_cluster_corr(data_table = sample,;       
        %put                          out_table  = example,;               
        %put                          agg_funcs  = mean p50 p75 max min;           
        %put          );                          
    %end;
    %else %if %lowcase(%trim(&m_name)) = gp_cluster_size %then %do;
        %put macro: gp_cluster_size;
        %put ;
        %put purpose: calculate statistics of clusters correlation by each;
        %put          step of decomposition;                               
        %put ;
        %put arguments:;
        %put        data_table  = : Source table name, should be the OUTSTATS FROM Varclus;
        %put        out_table   = : Output table of KDE. Ignored if missing;               
        %put        log_flag    = : 1 or 0(miss). If log transform the size of clusters;   
        %put        plot_flag   = : 1 or 0(miss). If plot the histogram and density;       
        %put        delete_flag = : 1 or 0(miss). Delete all temporary data after run;     
        %put        plot_clus_list = : A list of clusters represented by the number to be;
        %put                           plotted. If not given, a plot of 21 subplots based on;
        %put                           log interval will be generated;
        %put        plot_hist_xmax = : Max value of xaxis for the histogram plot. Default 20;                
        %put        plot_num_cols  = : Number of columns in the plots. Default 7;                       
        %put        plot_width     = : Plot width passed to ods graphics width option.;
        %put                           Default 6in;
        %put        plot_height    = : Plot height passed to ods graphics height option.;
        %put                           Default 8in;
        %put ;
        %put example:;
        %put          gp_cluster_size(data_table = sample,;              
        %put                          out_table  = kde_out,;               
        %put                          plot_flag  = 1,;                               
        %put                          plot_clus_list = 1 3 5 8 10 11 125;
        %put          );
    %end;
    %else %do;
        %put please input a valid macro name;
        %put current macros in the ordw library:;
        %put gp_info_value;
        %put gp_stability_chart;
        %put gp_stability_index;
        %put gp_var_split;
        %put gp_univariate_c;
        %put gp_cluster_corr;
        %put gp_cluster_size;
    %end;
%put ****************************************************************************************;
%mend gp_help;

*********************************************************************************************;
***                                   Private Macros                                      ***;    
*********************************************************************************************;

%macro _gp_num_fmt(func_num, data_table, var_list);
    proc univariate data = &data_table noprint;
        var &var_list;
        output out = &func_num._fmt 
               pctlpts = (1, 5 to 95 by 10, 99)
               pctlpre = %sysfunc(tranwrd(%sysfunc(compbl(&var_list))%str( ), %str( ), %str(_ )));
    run;

    proc transpose data = &func_num._fmt out = &func_num._fmt; run;
    
    data &func_num._fmt;
        set &func_num._fmt;
        length fmtname $36;
        fmtname = strip(scan(_label_, -1, ',')) || '_';
        value   = round(col1, .001);
        keep fmtname value;
    run;

    proc sort data = &func_num._fmt nodup; by fmtname value; run;

    * Create the format for the variables using percentile bins;
    data &func_num._fmt;
        length start end last_val $16 label $36;
        retain last_val;
        set &func_num._fmt;
        by fmtname;

        if first.fmtname then do;
            * Label for count of missings. Always executed;
            start = '.';
            label = 'missing';
            last_val = '';
            output; /* FMT for missing */

            if value ne . then do;
                start = 'low';
                end   = strip(put(value, best16.));
                label = catx(' ', '<=', end);
                last_val = end;
                output; /* FMT for <= min */
            end;
        end;
        else do; /* Supposedly because of prior sorting, remaining obs should be non-miss */
            end = strip(put(value, best16.));

            if missing(last_val) then do; /* Last obs was missing */
                start = 'low';
                label = catx(' ', '<=', end);
            end;
            else do;
                start = last_val;
                label = catx(' - ', start, end);
            end;

            last_val = end;
            output; /* FMT for val_1 ~~ val_2 */

            if last.fmtname then do;
                start = end;
                end   = 'high';
                label = catx(' ', '>', start);
                output; /* FMT for > max */
            end;
        end;
    run;
%mend _gp_num_fmt;

%macro _gp_diag9(data_table, out_table, clus_table);
    data &out_table;
        set _null_;
    run;

    proc sql noprint;
        SELECT DISTINCT clus INTO: clusters SEPARATED BY ' '
            FROM &clus_table
        ;

        SELECT DISTINCT _var_name INTO: var_list SEPARATED BY ' '
            FROM &clus_table
        ;

    quit;

    data _diag9_data;
        set &data_table;
        keep &var_list;
    run;

    sasfile _diag9_data open;

    %let i = 1;
    %do %while(%length(%scan(&clusters, &i)));
        %let clus = %scan(&clusters, &i);
        proc sql noprint;
            SELECT _var_name INTO: clus_vars SEPARATED BY ' '
                FROM &clus_table
                WHERE clus = &clus
            ;
        quit;

        proc corr data = &data_table.(keep = &clus_vars)
                  outp = _diag9
                  noprint;
        run;

        proc sort data = _diag9; by _name_; run;

        proc transpose data = _diag9 out = _diag9_t;
            by _name_;
            where _type_ = 'CORR';
        run;

        data _diag9_t;
            set _diag9_t;
            col1 = abs(col1);
        run;

        proc sort data = _diag9_t; by col1; run;

        data _diag9_t;
            set _diag9_t;
            clus = &clus;
            if _N_ = 1;
            drop _name_;
        run;

        data &out_table;
            set &out_table
                _diag9_t;
        run;
        %let i = %eval(&i + 1);
    %end;

    proc datasets lib = work nolist;
        delete _diag9:;
    quit;

    data &out_table;
        set &out_table;
        rename col1 = min_abs_corr;
    run;

    sasfile _diag9_data close;
%mend _gp_diag9;

%macro _gp_clus_score(func_num, data_table, oustats_table, out_table, select_clus);
data &func_num._score;
    set &oustats_table;
    where _NCL_ = &select_clus or _NCL_ = .;
    drop _NCL_;
run;

proc score data = &data_table
    score = &func_num._score
    out = &out_table.(keep = clus:);
run;


%mend _gp_clus_score;

options source;
