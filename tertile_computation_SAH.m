%% Compute MAP Reductions (run in SAH Sub Analysis Folder

main_dir = cd;
cd('Analysis\Epoch 60 min')
master_data = readtable('Intervention_Analysis_Epoch_60.csv');
cd(main_dir)
cohort = readtable('SAH_cohort_SELECTED.csv');

mean_map_reductions = []; % indices of this array matches cohort file

%%
for pt = 1:height(cohort) % iterate thru selected patients
    
    MR = cohort.MR(pt,1);
    
    admin_indexes = find(strcmp(master_data.MR,MR) );
    
    % calculate the mean map reduction across all administrations for this
    % patient
    mean_map_reductions(pt) = mean(master_data.MEAN_ABP_AFTER(admin_indexes) - master_data.MEAN_ABP_BEFORE(admin_indexes));
    
    
end



%% Calculate tertiles

min_map_reduction = min(mean_map_reductions)
first_breakpoint = quantile(mean_map_reductions,1/3)
second_breakpoint = quantile(mean_map_reductions,2/3)
max_map_reduction = max(mean_map_reductions)

first_tertile_inds = find(mean_map_reductions <= first_breakpoint); % GREATEST MAP REDUCTION
second_tertile_inds = find(mean_map_reductions > first_breakpoint & mean_map_reductions <= second_breakpoint);
third_tertile_inds = find(mean_map_reductions > second_breakpoint); % LEAST MAP REDUCTION ( SOME ARE POSITIVE )

%% Assign tertiles to patients in cohort file

assignment = zeros(height(cohort),1); % 1 x N array where N = number of patients. Each entry is the tertile assignment (1,2,3) for the patient 

assignment(first_tertile_inds) = 1; 
assignment(second_tertile_inds) = 2;
assignment(third_tertile_inds) = 3;

cohort = [cohort table(assignment)] % add assignments to cohort table

%% Assign tertiles to administrations in master_data file

admin_assignments = zeros(height(master_data),1); % assign each administration to a tertile based on patient assignment

for pt = 1:height(cohort) % iterate thru selected patients
    
    MR = cohort.MR(pt,1); % get MR for this patient
    assign = assignment(pt); % get tertile assignment for this patient
    
    admin_indexes = find(strcmp(master_data.MR,MR) ); % find indices of pt in master_data file
   
    
    admin_assignments(admin_indexes) = assign;
    
    
end

%% Results calculated from summary stats file


vars_of_interest_admins = {'PERCENT_BELOW_BEFORE', 'PERCENT_ABOVE_BEFORE', 'PERCENT_BELOW_AFTER', ...
                    'PERCENT_ABOVE_AFTER', 'MEAN_ABP_BEFORE',  'MEAN_ABP_AFTER',  'MEAN_COX_BEFORE', ...
                    'MEAN_COX_AFTER',  'MEAN_MAPopt_BEFORE', 'MEAN_MAPopt_AFTER', 'delta_MAP_AFTER', ...
                    'delta_MAP_BEFORE',...
                    'MEAN_NIRS_ON_BEFORE', 'MEAN_NIRS_ON_AFTER','MEAN_NIRS_OFF_BEFORE', 'MEAN_NIRS_OFF_AFTER',...
                    'MEAN_NIRS_ON_BELOW','MEAN_NIRS_ON_NOTBELOW','MEAN_NIRS_OFF_BELOW','MEAN_NIRS_OFF_NOTBELOW'};
                

   
data_package = cell2table(cell(0));

for tertile = 1:3
    
    
   
    row = {};
    vars = {};
    tertile_admins = master_data(admin_assignments == tertile,:);
    
    for var = 1:length(vars_of_interest_admins)
        
       % CALCULATING SUMMARY STATS
       var_name = vars_of_interest_admins{var};
       delta_calc = tertile_admins.(var_name); % all data for this var across relevant administrations
       row = [row {mean(delta_calc,'omitnan')}  {std(delta_calc,'omitnan')}]; % calcualte the mean and std for each var in vars_of_interest
       vars = [vars {[var_name,'_MEAN']} {[var_name,'_STD']}];
       
      
    end
    
  
    
    data_package = vertcat(data_package,row);
   

end


data_package.Properties.VariableNames = vars;

      
%% p-values :(

vars_of_interest = {'PERCENT_BELOW', 'PERCENT_ABOVE','MEAN_ABP','MEAN_COX','MEAN_MAPopt', 'delta_MAP',...
                                'MEAN_PRX','MEAN_NIRS_ON','MEAN_NIRS_OFF'};                



p_package = cell2table(cell(0,length(vars_of_interest))); p_package.Properties.VariableNames = vars_of_interest;

for tertile = 1:3

    tertile_admins = master_data(admin_assignments == tertile,:);
    p_values = {};
    
    for var = 1:length(vars_of_interest)
        
        before_data = tertile_admins.([vars_of_interest{var} '_BEFORE']);
        after_data = tertile_admins.([vars_of_interest{var} '_AFTER']);
        
        test_stat = after_data-before_data; [h,p] = adtest(test_stat); % null=normal dist!
        figure; histogram(test_stat); title(vars_of_interest{var})
        
        if h % data is normally distributed
           
            [h,p,ci,stats] = ttest(before_data,after_data); % t-test
            display([vars_of_interest{var} 'is normally distributed'])
            
        else
            
            [p,h,stats] = signrank(before_data,after_data);
            display([vars_of_interest{var} 'is NOT normally distributed'])
        end

        
        p_values = [p_values {p}];

    end
    
     p_package = vertcat(p_package,p_values);

  
end


%% OTHER p-values (not before vs after)

group_1 = {'MEAN_NIRS_ON_BELOW','MEAN_NIRS_OFF_BELOW'}
group_2 = {'MEAN_NIRS_ON_NOTBELOW','MEAN_NIRS_OFF_NOTBELOW'}
table_vars = {'MEAN_NIRS_ON_COMPARE','MEAN_NIRS_OFF_COMPARE'}


p_package_addendum = cell2table(cell(0,length(table_vars))); p_package_addendum.Properties.VariableNames = table_vars;

for tertile = 1:3

    tertile_admins = master_data(admin_assignments == tertile,:);
    p_values = {};
    
    for comparison =1:length(group_1)
   
        group1 = tertile_admins.(group_1{comparison}); % get data
        group2 = tertile_admins.(group_2{comparison});
        
        test_stat = group1-group2; [h,p] = adtest(test_stat); % null=normal dist!
%         figure; histogram(test_stat); title(vars_of_interest{var})
        
        if h % data is normally distributed
           
            [h,p,ci,stats] = ttest(group1,group2); % t-test
            display([table_vars{comparison} 'is normally distributed'])
            
        else
            
            [p,h,stats] = signrank(group1,group2);
            display([table_vars{comparison} 'is NOT normally distributed'])
        end

        
        p_values = [p_values {p}];

    end
    
    p_package_addendum = vertcat(p_package_addendum,p_values);

  
end



writetable([p_package p_package_addendum],"P_Results.csv") % each tertile is a new row. 1st row = 1st tertile. 
                            
                        


%% Cohort stuff

% based on tertile assignment, get indexes of cohort file for each tertile.
cohort_low = cohort(assignment == 3,:);
cohort_med = cohort(assignment == 2,:);
cohort_high = cohort(assignment == 1,:);

percent_below_time = cell(0,2);
mRS_fracs = cell(0,7);

for tertile = 1:3
    

   % fraction of viable monitoring time spent below
   cohort_set = cohort(assignment == tertile,:);
   calc = (cohort_set.total_time_below_sec ./ cohort_set.viable_monitoring_time_sec) * 100;
   percent_below = [{mean(calc,'omitnan')} {std(calc,'omitnan')}]; % proportion of viable monitoring time spent below LLA

   % mRS fractions
   N = height(cohort_set); % # of patients in this tertile cohort
   num_0 = sum(cohort_set.mRS == 0);
   num_1 = sum(cohort_set.mRS == 1);
   num_2 = sum(cohort_set.mRS == 2);
   num_3 = sum(cohort_set.mRS == 3);
   num_4 = sum(cohort_set.mRS == 4);
   num_5 = sum(cohort_set.mRS == 5);
   num_6 = sum(cohort_set.mRS == 6);


   mRS_row =  [ {(num_0 / N)*100} {(num_1 / N)*100} {(num_2 / N)*100} {(num_3 / N)*100} {(num_4 / N)*100} {(num_5 / N)*100} {(num_6 / N)*100} ];


    percent_below_time = [percent_below_time ; percent_below];
    mRS_fracs = [mRS_fracs ; mRS_row];
end

percent_below_table = cell2table(percent_below_time); percent_below_table.Properties.VariableNames = {'Percent_Below_Time_mean','Percent_Below_Time_Std'};
mRS_fracs = cell2table(mRS_fracs); mRS_fracs.Properties.VariableNames = {'percent_0','percent_1','percent_2','percent_3','percent_4','percent_5','percent_6',};
 
 
data_package = [data_package percent_below_table mRS_fracs]

writetable(data_package,"Results.csv") % each tertile is a new row. 1st row = 1st tertile. 
                            

%% anova

calc = (cohort.total_time_below_sec ./ cohort.viable_monitoring_time_sec) * 100;
[p,tbl,stats] = anova1(calc,cohort.assignment)

