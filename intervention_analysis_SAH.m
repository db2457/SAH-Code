%% Revision History
% 11/23: Added code to search for MAP or ABP. MAP exists in hemosphere
% files. ABP exists in ICM+ files.


cd 'C:\Users\db2457\OneDrive - Yale University\Desktop\Projects\Intervention Analysis\SAH Sub Analysis\SAH-Code'
code_folder = cd; % directory where code files are stored
cd ..

allowed_medications = { ...                 
'HYDRALAZINE 10 MG TABLET',...
'HYDRALAZINE 20 MG/ML INJECTION SOLUTION',...
'HYDRALAZINE 25 MG TABLET',...
'HYDRALAZINE 50 MG TABLET',...
'LABETALOL (NORMODYNE) 50 MG HALFTAB',...
'LABETALOL 100 MG TABLET',...
'LABETALOL 200 MG TABLET',...
'LABETALOL 300 MG TABLET',...
'LABETALOL 5 MG/ML INTRAVENOUS SOLUTION',...
'NIMODIPINE 30 MG CAPSULE',...
'NIMODIPINE 30 MG/10ML ORAL SOLUTION',...
'NIMODIPINE 30 MG/5ML ORAL SYRINGE (FOR ORAL USE ONLY)',...
'NIMODIPINE 60 MG/20ML ORAL SOLUTION'}; %EPIC_MED_NAME

viable_threshold = 0.80; %what proportion of data must be viable for window analysis?
below_threshold = 0.5; % min percent below of epoch time threhsold to be considered critical hypoperfusion

jdat_filename = "JDAT_SAH_only.csv";
cohort_filename = "SAH_cohort_ALL.xlsx";
eeg_overlap_filename = "EEG_overlap.csv";

epochs = [15, 30, 60, 90, 120]; % +/- mins

jdat = readtable(jdat_filename); jdat.MAR_TAKEN_TIME = datetime(jdat.MAR_TAKEN_TIME,'ConvertFrom','excel');
cohort = readtable(cohort_filename);
overlap = readtable(eeg_overlap_filename);
fs = 0.1;
off_injury = 0; % 0 = analyze injured sides. 1 = analyze uninjured sides
                
covariates_kim = {'event_id','mrn','date','PERCENT_BELOW_AFTER','mean_before','mean_after','median_before','median_after','min_before',...
                     'min_after','max_before','max_after'};  % DR. KIM CODE, ADDED 9/28
%               
% event_interest = overlap.event_id; % DR. KIM CODE, ADDED 9/28

temp_cohort = readtable('SAH_cohort_SELECTED.xlsx'); temp_cohort = temp_cohort.mr;



side = 1; % Which side is being analyzed. 0 = LEFT, 1 = RIGHT
%% 

kim_analysis = 1; % perform Dr. Kim's analyses


%% IMPORT ICM+ DATA

cd .. 
cd('CSV_Files\Cleaned_SAH')

ALL_DATA = {}; % cell array containing data of interest. each entry is another patient.
files = dir( fullfile(cd,'*.csv') );   %# list all .csv files in directory
files = {files.name}';                      %'# get file names
warning('off','MATLAB:table:ModifiedAndSavedVarnames') % turn annoying import warnings off

for patient = 1:height(cohort) % iterate through patients in cohort file 
  
    
     MR = cohort.mrn{patient};
     MR_files = split(files,'_'); % MRs of all CSV files
        
   
      
      if ~isnan(cohort.exclude_intervention(patient)) %   % if patient is excluded, don't bother importing... 
        reason = cohort.Status(patient);
        reason = reason{1};
        display([MR, ' excluded: ',reason])
        continue % skip this patient
      end 
        
     file_ind =  find(strcmp(MR, MR_files(:,1))); % index of filename that corresponds to this patient
      
      
      if isempty(file_ind) % if the ICM+ CSV file does not exist in file directory
        display([MR, ' does not have a .CSV file'])
        continue % skip this patient
        
      else % import it 
          
          data_table = readtable(files{file_ind});
          
      end

    
    % if actual CSV file is empty
    if isempty(data_table)
        disp([filename, ' is an empty CSV file'])
        ALL_DATA{patient} = []; % mark with empty spot
    end
    
    
    
    ALL_DATA{patient} = data_table; % append to cell array
           
end


%% ANALYSIS
cd(code_folder)
cd ..
mkdir('Analysis') % make new folder for cleaned files
cd('Analysis') % change directory to that folder
analysis_folder = cd;

analyzed_MRs = {};
dates = {};



for epoch = 3
    
    newdir_name = ['Epoch ', num2str(epochs(epoch)),' min'];
    mkdir(newdir_name)
    cd(newdir_name)    
    newdir = cd;
   
    
    if kim_analysis
        epoch_file_kim  = cell2table(cell(0,length(covariates_kim)));  epoch_file_kim.Properties.VariableNames = covariates_kim; % initialize empty epoch file
    end
    
    
% %    
    
    for patient = 1:length(ALL_DATA)
        
        % skip if empty
        if isempty(ALL_DATA{patient})
            continue
        end
        
        
        
        % collect information
        data = ALL_DATA{patient};
        time = data.DateTime;
        MR = cohort.mrn{patient};
        
        
        cohort_ind =  find(strcmp(cohort.mrn, MR)); % index of MR in cohort file
    


         
        % Pull NIRS and ICP-derived limits of autoregulation
        
        if any(contains(data.Properties.VariableNames,'LLA_PRx')) && any(contains(data.Properties.VariableNames,'ULA_PRx')) % do ICP-derived limits exist? 
            
            lla_icp = data.LLA_PRx;
            ula_icp = data.ULA_PRx;
      
            
        else
                             
            lla_icp = [];
            ula_icp = [];
            
            
        end
        
        
        if any(contains(data.Properties.VariableNames,'LLA_L')) && any(contains(data.Properties.VariableNames,'LLA_R')) % do NIRS-derived limits exist (just checks left)
            
 
            
            if side
                
                lla_nirs = data.LLA_R; 
                ula_nirs = data.ULA_R;
      
                
            else
                
                lla_nirs = data.LLA_L; 
                ula_nirs = data.ULA_L;
                
            end
         
 
        else
            
            lla_nirs = [];
            ula_nirs = [];
            cox = [];
   
            
        end
        
        % Choose whether we use NIRS or ICP-derived autoregulatory data
        % (ULA,LLA,MAPopt)
        
        if isempty(lla_icp) || isempty(ula_icp) % if no ICP data...
            
            if isempty(lla_nirs) || isempty(ula_nirs) % AND no NIRS data,,,
                
                display([MR, ' has no ICP OR NIRS-derived autoregulation data. This patient should be excluded from analysis'])
                continue
            end
                
           % use NIRS data
            upper_mean = table(ula_nirs); upper_mean.Properties.VariableNames = {'upper'};
            lower_mean = table(lla_nirs); lower_mean.Properties.VariableNames = {'lower'};


            if any(contains(data.Properties.VariableNames,'MAPopt_flex_R_mmHg_'))
                
                if side

                    mapopt_mean = table(data.MAPopt_flex_R_mmHg_); mapopt_mean.Properties.VariableNames = {'mapopt'};
                    
                    
                else
                    
                    mapopt_mean = table(data.MAPopt_flex_L); mapopt_mean.Properties.VariableNames = {'mapopt'};
                end

             

            else

                 if side

                    mapopt_mean = table(data.MAPopt_flex_R); mapopt_mean.Properties.VariableNames = {'mapopt'};
                    
                    
                else
                    
                    mapopt_mean = table(data.MAPopt_flex_L); mapopt_mean.Properties.VariableNames = {'mapopt'};
                 end
                
            end

            data = [data lower_mean upper_mean mapopt_mean];
            nirs_data_used_indicator = 1;
                
            
            
        elseif isempty(lla_nirs) || isempty(ula_nirs) % if we have ICP data but no NIRS..
            
            % use ICP data
             upper_limit = table(data.ULA_PRx); upper_limit.Properties.VariableNames = {'upper'};
             lower_limit = table(data.LLA_PRx); lower_limit.Properties.VariableNames = {'lower'};
             mapopt_icp = table(data.MAPopt_PRx); mapopt_icp.Properties.VariableNames = {'mapopt'};
             data = [data upper_limit lower_limit mapopt_icp];
            
        else % we have ICP and NIRS derived data
            
            % compare which is better
            
           
            nirs_quality_index = sum(isnan(lla_nirs)) + sum(isnan(ula_nirs)); % higher the number, the more missing NIRS data...
            icp_quality_index = sum(isnan(lla_icp)) + sum(isnan(ula_icp));

            
            if nirs_quality_index > icp_quality_index % if NIRS signal is worse...
               
                  % use ICP data
                 upper_limit = table(data.ULA_PRx); upper_limit.Properties.VariableNames = {'upper'};
                 lower_limit = table(data.LLA_PRx); lower_limit.Properties.VariableNames = {'lower'};
                 mapopt_icp = table(data.MAPopt_PRx); mapopt_icp.Properties.VariableNames = {'mapopt'};
                 data = [data upper_limit lower_limit mapopt_icp];

               
                
            else
                
                % use NIRS data
                upper_mean = table(ula_nirs); upper_mean.Properties.VariableNames = {'upper'};
                lower_mean = table(lla_nirs); lower_mean.Properties.VariableNames = {'lower'};
              
                
                if any(contains(data.Properties.VariableNames,'MAPopt_flex_R_mmHg_'))
                
                    if side

                        mapopt_mean = table(data.MAPopt_flex_R_mmHg_); mapopt_mean.Properties.VariableNames = {'mapopt'};


                    else

                        mapopt_mean = table(data.MAPopt_flex_L); mapopt_mean.Properties.VariableNames = {'mapopt'};
                    end

             

                else

                    if side

                        mapopt_mean = table(data.MAPopt_flex_R); mapopt_mean.Properties.VariableNames = {'mapopt'};


                    else

                        mapopt_mean = table(data.MAPopt_flex_L); mapopt_mean.Properties.VariableNames = {'mapopt'};
                    end
                
                end

                data = [data lower_mean upper_mean mapopt_mean];
                nirs_data_used_indicator = 1;

                
            end
            

        end
        
        
         
        % define ABP variable b.c. Hemosphere files have MAP not ABP :/
        if any(contains(data.Properties.VariableNames,'ABP'))
            ABP = 'ABP';
        elseif any(contains(data.Properties.VariableNames,'MAP'))
            ABP = 'MAP';
        end
        
        % raw NIRS signal
        NIRS_signal_exists = any(contains(data.Properties.VariableNames,'rSO2R')) && any(contains(data.Properties.VariableNames,'rSO2L'));
         if NIRS_signal_exists % if NIRS exist...
             
             if side
                 
                 nirs_mean = table(data.rSO2R); nirs_mean.Properties.VariableNames = {'nirs'};
                 cox_mean = table(data.COxr); cox_mean.Properties.VariableNames = {'cox'};
                 
             else
                 
                 nirs_mean = table(data.rSO2L); nirs_mean.Properties.VariableNames = {'nirs'};
                 cox_mean = table(data.COxl); cox_mean.Properties.VariableNames = {'cox'};
             end
             
        
              
         end
        
       

        
        % match with JDAT
        jdat_indexes = find(strcmp(jdat.MRN,MR) );
        admins = jdat(jdat_indexes,:); % slices JDAT data for just one patient
        admins = admins(contains(admins.EPIC_MED_NAME,allowed_medications),:); % select only those admins that are in the allowed medication list

   
        
        for admin = 1:height(admins) % ITERATE THROUGH MEDICATION ADMINISTRATIONS FOR THIS PATIENT--------------------------
            
           
            
           % collect information
            row = admins(admin,:);
            med = row.EPIC_MED_NAME;
            time_taken = row.MAR_TAKEN_TIME;
            dose = row.DOS;
            event_id = row.EVENT_ID;
            
            % Pull time indexes from ICM data
            start_index = find(abs(time-time_taken) <= seconds(5)); %  time index of dose
            
            if isempty(start_index) % if dose time is out of bounds, skip to next admin
                continue
            else
                start_index = start_index(1);
            end
        
             % window data into before and after epochs
            cushion = round(fs * (epochs(epoch) * 60)); % index size of cushion
          
            
            
            if start_index + cushion > height(data) % check for out of right bounds
                
                continue; % skip this admin because this epoch is inappropriate
            
            elseif start_index - cushion <= 1 % check for out of left bounds
             
                continue; % skip this admin because this epoch is inappropriate
            
                             
            else
                
                 win_before = data(start_index - cushion : start_index,:);
                 win_after = data(start_index : start_index + cushion,:);
             
            end
            
%              11.22.2021: Dr. Kim only wanted the data requirement to
%               apply to the half hour before, rather than a full hour
%              before. 

          
                
            win_before_half = win_before(round((height(win_before))/2):end,:);
            viable_before = (sum(~any(isnan([win_before_half.lower,win_before_half.upper,win_before_half.(ABP)]),2)) / fs) / (height(win_before_half)/fs);
                
         
           
            viable_after = (sum(~any(isnan([win_after.lower,win_after.upper,win_after.(ABP)]),2)) / fs) / (height(win_after)/fs); 
            
           
            if viable_before < viable_threshold || viable_after < viable_threshold % if >90% of data is availabe on both sides

%                 fprintf(['\n',MR, ', Event ID: ', num2str(event_id),', Side: ',num2str(side),', Epoch: ', num2str(epoch),' is not viable for analysis (',num2str(viable_before*100),',',num2str(viable_after*100),')']);

                continue % skip this med admin
            end
            
          % -CALCULATE OVARIATES -----------------------------
      
            
         
            
            cd(code_folder) % to access needed .m files
            
    
            

            
            [PERCENT_BELOW_BEFORE, PERCENT_ABOVE_BEFORE] = outside_limits(win_before,epochs(epoch)*60);
            [PERCENT_BELOW_AFTER, PERCENT_ABOVE_AFTER,TIME_BELOW] = outside_limits(win_after,epochs(epoch)*60);
           
           

            MEAN_ABP_BEFORE = mean(win_before.(ABP),'omitnan');
            MEDIAN_ABP_BEFORE = median(win_before.(ABP),'omitnan');
            MAX_ABP_BEFORE = max(win_before.(ABP));
            MIN_ABP_BEFORE = min(win_before.(ABP));
            
            MEAN_ABP_AFTER = mean(win_after.(ABP),'omitnan');
            MEDIAN_ABP_AFTER = median(win_after.(ABP),'omitnan');
            MAX_ABP_AFTER = max(win_after.(ABP));
            MIN_ABP_AFTER = min(win_after.(ABP));

  
             
           cd(newdir)
              %% DR. KIM CODE, ADDED 9/28
           if epoch == 3 && kim_analysis
               
               
               % package requested covariatses
                covar_kim_data = {event_id,MR,time_taken,PERCENT_BELOW_AFTER*100, MEAN_ABP_BEFORE, MEAN_ABP_AFTER, MEDIAN_ABP_BEFORE, MEDIAN_ABP_AFTER,...
                                  MIN_ABP_BEFORE,MIN_ABP_AFTER,MAX_ABP_BEFORE,MAX_ABP_AFTER};
                epoch_file_kim = [epoch_file_kim ; covar_kim_data];
                
                % create and save requested plots
                concat_data = [win_before ; win_after]; % vert concat before and after data windows
                time_mod = concat_data.DateTime - time_taken;    time_mod.Format = 'm';
                x_ticks = time_mod(1):seconds(5*60):time_mod(end);
                
                time_mod.Format = 'm'; x_ticks = round(x_ticks,'minutes');
               
            
                
                fig = figure('Renderer', 'painters', 'Position', [10 10 1000 500],'visible','off');
                plot(time_mod,concat_data.(ABP),'black','LineWidth',1.5); hold on;
                plot(time_mod,concat_data.lower,'Color',[0 0.4470 0.7410]);
                plot(time_mod,concat_data.upper,'Color',[0 0.4470 0.7410]);
                xticks(x_ticks)
                
                xlabel('Time (min)'); ylabel('MAP (mmHg)'); title(['Event: ', num2str(event_id),', Time taken: ',datestr(time_taken)])
                saveas(fig,[num2str(event_id),'.jpg'])
     
                    
                time_mod_table = table(time_mod);  time_mod_table.Properties.VariableNames = {['time']};
                abp = table(concat_data.(ABP)); abp.Properties.VariableNames = {['abp']};
                lla = table(concat_data.lower); lla.Properties.VariableNames = {['lla']};
                ula = table(concat_data.upper); ula.Properties.VariableNames = {['ula']};
                    
                master_table = [time_mod_table,abp,lla,ula];
                writetable(master_table,[num2str(event_id),'.csv'])
                
           end
          
   
        end
           
           % PACKAGE DATA
           

       
        
  
    end
    
    

       
    
        
    writetable(epoch_file_kim,['Intervention_Analysis_KIM_Epoch_', num2str(epochs(epoch)),'.csv']) % DR. KIM CODE, ADDED 9/28

  
    display(["Just packaged Epoch ", num2str(epochs(epoch))])
    cd(analysis_folder)
    
end


cd(code_folder)





