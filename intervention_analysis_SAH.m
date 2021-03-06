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

jdat_filename = "JDAT_SAH_only_EDIT.csv";
cohort_filename = "SAH_cohort_ALL.xlsx";
eeg_overlap_filename = "EEG_overlap.csv";

epochs = [15, 30, 60, 90, 120]; % +/- mins

jdat = readtable(jdat_filename); jdat.MAR_TAKEN_TIME = datetime(jdat.MAR_TAKEN_TIME,'ConvertFrom','excel');
cohort = readtable(cohort_filename);
overlap = readtable(eeg_overlap_filename);
fs = 0.1;


covariates = {'event_id', 'MR', 'med', 'dose', 'PERCENT_BELOW_BEFORE', 'PERCENT_ABOVE_BEFORE', 'PERCENT_BELOW_AFTER', ...
                    'PERCENT_ABOVE_AFTER', 'MEAN_ABP_BEFORE',  'MEAN_ABP_AFTER',  'MEAN_COX_BEFORE', ...
                    'MEAN_COX_AFTER',  'MEAN_MAPopt_BEFORE', 'MEAN_MAPopt_AFTER', 'delta_MAP_AFTER','delta_MAP_BEFORE','MEAN_PRX_BEFORE', 'MEAN_PRX_AFTER',...
                    'MEAN_NIRS_BEFORE', 'MEAN_NIRS_AFTER',...
                    'MEAN_LLA_BEFORE',   'MEAN_LLA_AFTER',...
                    'MEAN_ULA_BEFORE','MEAN_ULA_AFTER',...
                    'time_taken'};
                
                 
    
                
covariates_kim = {'event_id','mrn','date','PERCENT_BELOW_AFTER','mean_before','mean_after','median_before','median_after','min_before',...
                     'min_after','max_before','max_after'};  % DR. KIM CODE, ADDED 9/28
%               
% event_interest = overlap.event_id; % DR. KIM CODE, ADDED 9/28


temp_cohort = readtable('SAH_cohort_SELECTED.xlsx'); temp_cohort = temp_cohort.mr;

num_events = zeros(1,length(temp_cohort))';
num_responses = zeros(1,length(temp_cohort))';
recording_time = zeros(1,length(temp_cohort))'; % total recroding time
viable_recording_time = zeros(1,length(temp_cohort))'; % viable recroding time (ALL signals exist)
nirs_data_used = zeros(1,length(temp_cohort))'; % =1 if NIRS-derived autoregulatory data was used. =0 for ICP-derived. 
times_below = cell(1,length(temp_cohort))'; 
map_drops = cell(1,length(temp_cohort))'; 
overall_time_below = zeros(1,length(temp_cohort))';
overall_time_above = zeros(1,length(temp_cohort))';
mean_LA_range = zeros(1,length(temp_cohort))';
mean_LLA = zeros(1,length(temp_cohort))';
mean_ULA = zeros(1,length(temp_cohort))';
mean_TOx = zeros(1,length(temp_cohort))';
mean_MAPopt = zeros(1,length(temp_cohort))';
mean_ABP = zeros(1,length(temp_cohort))';
mean_PRx = zeros(1,length(temp_cohort))';



%% 
selected_analysis = 1; % calculate summary statistics for a selected cohort
kim_analysis = 0; % perform Dr. Kim's analyses


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

cd(code_folder)

%% ANALYSIS
cd ..
mkdir('Analysis') % make new folder for cleaned files
cd('Analysis') % change directory to that folder
analysis_folder = cd;

analyzed_MRs = {};
dates = {};

NIRS_on_series_1 = [];
NIRS_on_series_2 = []; 
NIRS_on_series_3 = [];

NIRS_off_series_1 = [];
NIRS_off_series_2 = [];
NIRS_off_series_3 = [];

ABP_series_1 = [];
ABP_series_2 = [];
ABP_series_3 = [];

nirs_data_used_indicator = 0; % switches to 1 if nirs data is used

for epoch = 3
    
    newdir_name = ['Epoch ', num2str(epochs(epoch)),' min'];
    mkdir(newdir_name)
    cd(newdir_name)    
    newdir = cd;
   
    epoch_file = cell2table(cell(0,length(covariates)));  epoch_file.Properties.VariableNames = covariates; % initialize empty epoch file
    
    if kim_analysis
        epoch_file_kim  = cell2table(cell(0,length(covariates_kim)));  epoch_file_kim.Properties.VariableNames = covariates_kim; % initialize empty epoch file
    end
    
    
% %    
    
    for patient = 1:length(ALL_DATA)
        
        nirs_data_used_indicator = 0; % switches to 1 if nirs data is used
        
        % skip if empty
        if isempty(ALL_DATA{patient})
            continue
        end
        
        
        
        % collect information
        data = ALL_DATA{patient};
        time = data.DateTime;
        MR = cohort.mrn{patient};
        
        
        cohort_ind =  find(strcmp(cohort.mrn, MR)); % index of MR in cohort file
    
        
        if selected_analysis 
            
            sub_ind = regexpi(MR,'\([\d]\)'); % match (/d)

            if ~isempty(sub_ind) % if a (/d) exists in the MR

                MR = MR(1: sub_ind - 1);

            end
        
            MR_index = find(strcmp(temp_cohort,MR)); % index of MR in temp_cohort
            recording_time(MR_index) = height(data) / fs; % total monitoring time (s)
       
            

        
            if ~any( strcmp(temp_cohort,MR) ) % if MR does not exist in temp cohort...
                continue %skip patinet
            end
            
        end

        
  

        % Pull NIRS and ICP-derived limits of autoregulation
        
        if any(contains(data.Properties.VariableNames,'LLA_PRx')) && any(contains(data.Properties.VariableNames,'ULA_PRx')) % do ICP-derived limits exist? 
            
            lla_icp = data.LLA_PRx;
            ula_icp = data.ULA_PRx;
      
            
        else
                             
            lla_icp = [];
            ula_icp = [];
            
            
        end
        
        
        if any(contains(data.Properties.VariableNames,'LLA_L')) && any(contains(data.Properties.VariableNames,'LLA_R')) % do NIRS-derived limits exist (just checks left)
            
            % average left and right sides NIRS..
            lla_nirs = mean([data.LLA_L , data.LLA_R],2,'omitnan');
            ula_nirs = mean([data.ULA_L , data.ULA_R],2,'omitnan');
      
 
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
                
 

                mapopt_mean = table(mean([data.MAPopt_flex_L, data.MAPopt_flex_R_mmHg_],2,'omitnan')); mapopt_mean.Properties.VariableNames = {'mapopt'};

            else

                mapopt_mean = table(mean([data.MAPopt_flex_L, data.MAPopt_flex_R],2,'omitnan')); mapopt_mean.Properties.VariableNames = {'mapopt'};
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
                    
                    mapopt_mean = table(mean([data.MAPopt_flex_L, data.MAPopt_flex_R_mmHg_],2,'omitnan')); mapopt_mean.Properties.VariableNames = {'mapopt'};
                    
                else

                    mapopt_mean = table(mean([data.MAPopt_flex_L, data.MAPopt_flex_R],2,'omitnan')); mapopt_mean.Properties.VariableNames = {'mapopt'};
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
        
        % average left and right sides NIRS
        NIRS_signal_exists = any(contains(data.Properties.VariableNames,'rSO2R')) && any(contains(data.Properties.VariableNames,'rSO2L'));
         if NIRS_signal_exists % if NIRS exist...
        
              nirs_mean = table(mean([data.rSO2R , data.rSO2L],2,'omitnan')); nirs_mean.Properties.VariableNames = {'nirs'};
              data = [data nirs_mean];
              
              cox_mean = table(mean([data.COxl, data.COxr],2,'omitnan')); cox_mean.Properties.VariableNames = {'cox'};
              data = [data cox_mean];
              
         end
        
       

        
        % match with JDAT
        sub_ind = regexpi(MR,'\([\d]\)'); % match (/d)
        
        if ~isempty(sub_ind) % if a (/d) exists in the MR
           
            MR = MR(1: sub_ind - 1);
            
        end
        
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

            if kim_analysis
                
                win_before_half = win_before(round((height(win_before))/2):end,:);
                viable_before = (sum(~any(isnan([win_before_half.lower,win_before_half.upper,win_before_half.(ABP)]),2)) / fs) / (height(win_before_half)/fs);
                
            else
                
                %              % check data viability (%)
                viable_before = (sum(~any(isnan([win_before.lower,win_before.upper,win_before.(ABP)]),2)) / fs) / (height(win_before)/fs);
                
            end
            
            
           
            viable_after = (sum(~any(isnan([win_after.lower,win_after.upper,win_after.(ABP)]),2)) / fs) / (height(win_after)/fs); 
            
           
            if viable_before < viable_threshold || viable_after < viable_threshold % if >90% of data is availabe on both sides

%                 fprintf(['\n',MR, ', Event ID: ', num2str(event_id),', Side: ',num2str(side),', Epoch: ', num2str(epoch),' is not viable for analysis (',num2str(viable_before*100),',',num2str(viable_after*100),')']);

                continue % skip this med admin
            end
            
          % -CALCULATE OVARIATES -----------------------------
      
            
         
            
            cd(code_folder) % to access needed .m files
            
    
            
        
            
            [PERCENT_BELOW_BEFORE, PERCENT_ABOVE_BEFORE] = outside_limits(win_before,epochs(epoch)*60);
            [PERCENT_BELOW_AFTER, PERCENT_ABOVE_AFTER,TIME_BELOW] = outside_limits(win_after,epochs(epoch)*60);
       
            % don't really need to calcualte effects manually like this...
%             MEAN_NIRS_ON_LA_EFFECT = MEAN_NIRS_ON_BELOW - MEAN_NIRS_ON_NOTBELOW;
%             MEAN_NIRS_OFF_LA_EFFECT = MEAN_NIRS_OFF_BELOW - MEAN_NIRS_OFF_NOTBELOW;
            
           
         
           
           
           

            MEAN_ABP_BEFORE = mean(win_before.(ABP),'omitnan');
            MEDIAN_ABP_BEFORE = median(win_before.(ABP),'omitnan');
            MAX_ABP_BEFORE = max(win_before.(ABP));
            MIN_ABP_BEFORE = min(win_before.(ABP));
            
            MEAN_ABP_AFTER = mean(win_after.(ABP),'omitnan');
            MEDIAN_ABP_AFTER = median(win_after.(ABP),'omitnan');
            MAX_ABP_AFTER = max(win_after.(ABP));
            MIN_ABP_AFTER = min(win_after.(ABP));

            MEAN_LLA_BEFORE = mean(win_before.lower,'omitnan');
            MEAN_LLA_AFTER = mean(win_after.lower,'omitnan');
            
            MEAN_ULA_BEFORE = mean(win_before.upper,'omitnan');
            MEAN_ULA_AFTER = mean(win_after.upper,'omitnan');
            
            if NIRS_signal_exists
                
                MEAN_COX_BEFORE = mean(win_before.cox,'omitnan');
                MEAN_COX_AFTER = mean(win_after.cox,'omitnan');
                
                MEAN_NIRS_BEFORE = mean(win_before.nirs,'omitnan');
                MEAN_NIRS_AFTER = mean(win_after.nirs,'omitnan');
            
            else
                
                    
                MEAN_COX_BEFORE = NaN;
                MEAN_COX_AFTER = NaN;
                
                MEAN_NIRS_BEFORE = NaN;
                MEAN_NIRS_AFTER = NaN;
            
                
                
            end
            

              

            if any(contains(data.Properties.VariableNames,'PRx')) && any(~isnan(data.PRx)) % if PRx is a variable (isnt for hemosphere data)
                  MEAN_PRX_BEFORE = mean(win_before.PRx(~isnan(win_before.PRx)));
                  MEAN_PRX_AFTER = mean(win_after.PRx(~isnan(win_after.PRx)));
            else
                 MEAN_PRX_BEFORE = NaN;
                 MEAN_PRX_AFTER = NaN;
            end
            
            
     
           MEAN_MAPopt_BEFORE = mean(win_before.mapopt,'omitnan');
           MEAN_MAPopt_AFTER = mean(win_after.mapopt,'omitnan');
           
            
           delta_MAP_AFTER =  MEAN_ABP_AFTER -   MEAN_MAPopt_AFTER; % postive values indicate ABP is too high. 
           delta_MAP_BEFORE = MEAN_ABP_BEFORE -   MEAN_MAPopt_BEFORE;
          
           delta_COX =   MEAN_COX_AFTER -  MEAN_COX_BEFORE; % positive values indicate worsening in autoregulatory function
           delta_PRX = MEAN_PRX_AFTER - MEAN_PRX_BEFORE;
           
         
           if selected_analysis % entered for every med admin
                
             
               if epoch == 3 
                 
                   if PERCENT_BELOW_AFTER >= below_threshold % if they spent more than 50% of their time after below LLA....
                       MR_index = find(strcmp(temp_cohort,MR));
                       num_responses(MR_index) =  num_responses(MR_index) + 1; % then classify this admin as a response

                   end
                   
                  

                   MR_index = find(strcmp(temp_cohort,MR)); % index of MR in temp_cohort
                   num_events(MR_index) = num_events(MR_index) + 1; % iteratively add 
                   times_below{MR_index} = [times_below{MR_index} TIME_BELOW]; % iteratively add 
                   map_drops{MR_index} = [map_drops{MR_index} MEAN_ABP_AFTER - MEAN_ABP_BEFORE];
                   
                   if nirs_data_used_indicator % NIRS-derived autoregulatory data was used for this patient
                       nirs_data_used(MR_index) = 1;
                   end
             
                   % covariates for recording as a whole (not with respect
                   % to epoch)
                   viable_recording_time(MR_index) = sum(~any(isnan([data.lower,data.upper,data.(ABP)]),2)) / fs;
                   [~,~,time_below] = outside_limits(data,epochs(epoch)*60); overall_time_below(MR_index) = time_below; % apply outside_limits to the entire data (not just a window)
                   [~,~,~,time_above] = outside_limits(data,epochs(epoch)*60); overall_time_above(MR_index) = time_above; % apply outside_limits to the entire data (not just a window)
                   mean_LA_range(MR_index) = mean(data.upper - data.lower,'omitnan');
                   mean_LLA(MR_index) = mean(data.lower,'omitnan');
                   mean_ULA(MR_index) = mean(data.upper,'omitnan');
                   mean_MAPopt(MR_index) = mean(data.mapopt,'omitnan');
                   mean_ABP(MR_index) = mean(data.(ABP),'omitnan');
                   
                   
                   
                   
                   if NIRS_signal_exists

                        mean_TOx(MR_index) = mean(data.cox,'omitnan');
                  else


                        mean_TOx(MR_index) = NaN;


                   end
                   
                
                   if any(contains(data.Properties.VariableNames,'PRx')) && any(~isnan(data.PRx)) % if PRx is a variable (isnt for hemosphere data)
                   
                       mean_PRx(MR_index) = mean(data.PRx,'omitnan');
                  
                   else
                       
                       mean_PRx(MR_index) = NaN;
                       
                   end

                      
          
                   
                   
               end
                
               
%                if tertile == 1 Temporal NIRS series. (code here is not
%                changed to reflect averaging of L&R NIRS signals)
%                    
%                    
%                     NIRS_off_series_1 = vertcat(NIRS_off_series_1,[win_before.(nirs_off)' win_after.(nirs_off)']);
%                     NIRS_on_series_1 = vertcat(NIRS_on_series_1,[win_before.(nirs_on)' win_after.(nirs_on)']);
%                     ABP_series_1 = vertcat(ABP_series_1,[win_before.(ABP)' win_after.(ABP)']);
%                    
%                elseif tertile == 2
%                    
%                     NIRS_off_series_2 = vertcat(NIRS_off_series_2,[win_before.(nirs_off)' win_after.(nirs_off)']);
%                     NIRS_on_series_2 = vertcat(NIRS_on_series_2,[win_before.(nirs_on)' win_after.(nirs_on)']);
%                     ABP_series_2 = vertcat(ABP_series_2,[win_before.(ABP)' win_after.(ABP)']);
%                    
%                    
%                else
%                    
%                     NIRS_off_series_3 = vertcat(NIRS_off_series_3,[win_before.(nirs_off)' win_after.(nirs_off)']);
%                     NIRS_on_series_3 = vertcat(NIRS_on_series_3,[win_before.(nirs_on)' win_after.(nirs_on)']);
%                     ABP_series_3 = vertcat(ABP_series_3,[win_before.(ABP)' win_after.(ABP)']);
%                    
%                    
%                end
               cd(newdir)
               
            end
           
           % PACKAGE DATA
           

           
           %% DR. KIM CODE, ADDED 9/28
           if epoch == 3 && kim_analysis 
              
               % ismember(event_id,overlap.Var1)
               
               
               % package requested covariatses
                covar_kim_data = {event_id,MR,time_taken,PERCENT_BELOW_AFTER*100, MEAN_ABP_BEFORE, MEAN_ABP_AFTER, MEDIAN_ABP_BEFORE, MEDIAN_ABP_AFTER,...
                                  MIN_ABP_BEFORE,MIN_ABP_AFTER,MAX_ABP_BEFORE,MAX_ABP_AFTER};
                
                epoch_file_kim = [epoch_file_kim ; covar_kim_data];
                
                % create and save requested plots
                concat_data = [win_before ; win_after]; % vert concat before and after data windows
                time_mod = concat_data.DateTime - time_taken;    time_mod.Format = 'm';
                x_ticks = time_mod(1):seconds(5*60):time_mod(end);
                
                time_mod.Format = 'm'; x_ticks = round(x_ticks,'minutes');
               
            
                
                fig = figure('Renderer', 'painters', 'Position', [10 10 1000 500]);
                plot(time_mod,concat_data.(ABP),'black','LineWidth',1.5); hold on;
                plot(time_mod,concat_data.lower,'Color',[0 0.4470 0.7410]);
                plot(time_mod,concat_data.upper,'Color',[0 0.4470 0.7410]);
                plot(time_mod,concat_data.mapopt);
                xticks(x_ticks)

                
                xlabel('Time (min)'); ylabel('MAP (mmHg)'); title(['Event: ', num2str(event_id),', Time taken: ',datestr(time_taken)])
                saveas(fig,[num2str(event_id),'.jpg'])
                close all
     
                    
                time_mod_table = table(time_mod);  time_mod_table.Properties.VariableNames = {['time']};
                abp = table(concat_data.(ABP)); abp.Properties.VariableNames = {['abp']};
                lla = table(concat_data.lower); lla.Properties.VariableNames = {['lla']};
                ula = table(concat_data.upper); ula.Properties.VariableNames = {['ula']};
                mapopt = table(concat_data.mapopt); mapopt.Properties.VariableNames = {['mapopt']};
                    
                master_table = [time_mod_table,abp,lla,ula,mapopt];
                writetable(master_table,[num2str(event_id),'.csv'])
                
           end
          
            covar_data = {event_id, MR, med, dose, PERCENT_BELOW_BEFORE*100, PERCENT_ABOVE_BEFORE*100, PERCENT_BELOW_AFTER*100, ...
            PERCENT_ABOVE_AFTER*100, MEAN_ABP_BEFORE,  MEAN_ABP_AFTER,  MEAN_COX_BEFORE, ...
            MEAN_COX_AFTER,  MEAN_MAPopt_BEFORE, MEAN_MAPopt_AFTER, delta_MAP_AFTER,delta_MAP_BEFORE,MEAN_PRX_BEFORE, MEAN_PRX_AFTER,...
            MEAN_NIRS_BEFORE, MEAN_NIRS_AFTER,...
            MEAN_LLA_BEFORE,   MEAN_LLA_AFTER,...
            MEAN_ULA_BEFORE,   MEAN_ULA_AFTER,...
            time_taken};

           
            epoch_file = [epoch_file ; covar_data]; % vertcat this admin onto the epoch file
           
        end
        
  
    end
    
    
    cd(analysis_folder)
    writetable(epoch_file,['Intervention_Analysis_Epoch_', num2str(epochs(epoch)),'.csv'])
    
    if kim_analysis
        
        writetable(epoch_file_kim,['Intervention_Analysis_KIM_Epoch_', num2str(epochs(epoch)),'.csv']) % DR. KIM CODE, ADDED 9/28
        
    end
    
  
    display(['Just packaged Epoch ', num2str(epochs(epoch)), ' with ', num2str(height(epoch_file)), ' events.'])
    cd(analysis_folder)
    
end


cd(code_folder)


%%

%%
% Temporal NIRS series. (code here is not
%                changed to reflect averaging of L&R NIRS signals)
times_below_mean = cellfun(@mean,times_below);
map_drops_mean = cellfun(@mean, map_drops);


times_below_std = cellfun(@std,times_below);
% 
% NIRS_on_series_1_mean = mean(NIRS_on_series_1,1,'omitnan')';
% NIRS_on_series_2_mean = mean(NIRS_on_series_2,1,'omitnan')';
% NIRS_on_series_3_mean = mean(NIRS_on_series_3,1,'omitnan')';
% 
% NIRS_off_series_1_mean = mean(NIRS_off_series_1,1,'omitnan')';
% NIRS_off_series_2_mean = mean(NIRS_off_series_2,1,'omitnan')';
% NIRS_off_series_3_mean = mean(NIRS_off_series_3,1,'omitnan')';
% 
% NIRS_on_series_1_std = std(NIRS_on_series_1,1,'omitnan')';
% NIRS_on_series_2_std = std(NIRS_on_series_2,1,'omitnan')';
% NIRS_on_series_3_std = std(NIRS_on_series_3,1,'omitnan')';
% 
% NIRS_off_series_1_std =  std(NIRS_off_series_1,1,'omitnan')';
% NIRS_off_series_2_std =  std(NIRS_off_series_2,1,'omitnan')';
% NIRS_off_series_3_std =  std(NIRS_off_series_3,1,'omitnan')';
% 
% ABP_1 = std(ABP_series_1,1,'omitnan')';
% ABP_2 = std(ABP_series_2,1,'omitnan')';
% ABP_3 = std(ABP_series_3,1,'omitnan')';
% 



export = [num_events num_responses recording_time viable_recording_time overall_time_below ...
          overall_time_above mean_LA_range mean_LLA mean_ULA mean_TOx mean_MAPopt mean_ABP ...
          times_below_mean times_below_std map_drops_mean mean_PRx nirs_data_used];


