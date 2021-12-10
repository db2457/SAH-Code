%% Revision History
% 11/23: Added code to search for MAP or ABP. MAP exists in hemosphere
% files. ABP exists in ICM+ files.

% MAKE SURE IN MAIN DIRECTORY BEFORE RUNNING

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
cohort_filename = "SAH_cohort_ALL.csv";
eeg_overlap_filename = "EEG_overlap.csv";

epochs = [15, 30, 60, 90, 120]; % +/- mins

jdat = readtable(jdat_filename); jdat.MAR_TAKEN_TIME = datetime(jdat.MAR_TAKEN_TIME,'ConvertFrom','excel');
cohort = readtable(cohort_filename);
overlap = readtable(eeg_overlap_filename);
old_folder = cd;
fs = 0.1;
off_injury = 0; % 0 = analyze injured sides. 1 = analyze uninjured sides

covariates = {'event_id', 'MR', 'med', 'dose', 'PERCENT_BELOW_BEFORE', 'PERCENT_ABOVE_BEFORE', 'PERCENT_BELOW_AFTER', ...
                    'PERCENT_ABOVE_AFTER', 'MEAN_ABP_BEFORE',  'MEAN_ABP_AFTER',  'MEAN_COX_BEFORE', ...
                    'MEAN_COX_AFTER',  'MEAN_MAPopt_BEFORE', 'MEAN_MAPopt_AFTER', 'delta_MAP_AFTER','delta_MAP_BEFORE','MEAN_PRX_BEFORE', 'MEAN_PRX_AFTER',...
                    'MEAN_NIRS_ON_BEFORE', 'MEAN_NIRS_ON_AFTER','NIRS_ON_EFFECT',...
                    'MEAN_NIRS_OFF_BEFORE', 'MEAN_NIRS_OFF_AFTER','NIRS_OFF_EFFECT',...
                    'MEAN_NIRS_ON_BELOW','MEAN_NIRS_ON_NOTBELOW','MEAN_NIRS_ON_LA_EFFECT',...
                    'MEAN_NIRS_OFF_BELOW','MEAN_NIRS_OFF_NOTBELOW','MEAN_NIRS_OFF_LA_EFFECT',...
                    'time_taken'};
                    
                    
         
covariates_kim = {'event_id','mrn','date','PERCENT_BELOW_AFTER','mean_before','mean_after','median_before','median_after','min_before',...
                    'min_after','max_before','max_after'};  % DR. KIM CODE, ADDED 9/28
              
% event_interest = overlap.event_id; % DR. KIM CODE, ADDED 9/28

temp_cohort = readtable('SAH_cohort_SELECTED.csv'); temp_cohort_assignments = temp_cohort.tertile_assignment; temp_cohort = temp_cohort.MR;
num_events = zeros(1,length(temp_cohort))';
num_responses = zeros(1,length(temp_cohort))';
recording_time = zeros(1,length(temp_cohort))'; % total recroding time
viable_recording_time = zeros(1,length(temp_cohort))'; % viable recroding time (ALL signals exist)
times_below = cell(1,length(temp_cohort))'; 
overall_time_below = zeros(1,length(temp_cohort))';

%% 
selected_analysis = 0;


%% IMPORT ICM+ DATA

cd .. 
cd('CSV_Files\Cleaned_SAH')

ALL_DATA = {}; % cell array containing data of interest. each entry is another patient.
files = dir( fullfile(cd,'*.csv') );   %# list all .csv files in directory
files = {files.name}';                      %'# get file names
warning('off','MATLAB:table:ModifiedAndSavedVarnames') % turn annoying import warnings off

for file = 1:length(files) % iterate through patients
  
    filename = files{file}; 
    table = readtable(filename); % read CSV file
    
    % var check
    
%     if ~any(contains(table.Properties.VariableNames,'LLA_R')) % quick and dirty fix to exclude files without nirs data
%         display([filename, 'has no LA data'])
%         continue
%     end
%     
    if isempty(table)
        disp([filename, 'is empty'])
        continue
    end
    
    ALL_DATA{file} = table; % append to cell array
           
end

cd(old_folder)

%% ANALYSIS

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

for epoch = 3
    
    newdir_name = ['Epoch ', num2str(epochs(epoch)),' min'];
    mkdir(newdir_name)
    cd(newdir_name)    
    newdir = cd;
   
    epoch_file = cell2table(cell(0,length(covariates)));  epoch_file.Properties.VariableNames = covariates; % initialize empty epoch file
    epoch_file_kim  = cell2table(cell(0,length(covariates_kim)));  epoch_file_kim.Properties.VariableNames = covariates_kim; % initialize empty epoch file
%    
    
    for patient = 1:length(ALL_DATA)
        
        % collect information
        data = ALL_DATA{patient};

        
        if isempty(data)
     
           continue %skip patient 
           
        end
        
        time = data.DateTime;
        split_filename = split(files{patient},'_');
        MR = split_filename{1};% assumes format MR_ ...
        
  
        
        cohort_ind =  find(strcmp(cohort.mrn, MR)); % index of MR in cohort file
        
        if isempty(cohort_ind) % if the ICM+ CSV file does not exist in cohort
            display([MR, ' does not exist in the cohort'])
            continue % skip this patient
        end
        
        if ~isnan(cohort.Exclude(cohort_ind)) % if this patient is excluded per cohort file... (good files have NaNs)
            reason = cohort.Status(cohort_ind);
            reason = reason{1};
            display([MR, ' excluded: ',reason])
            continue % skip this patient
        end
        
        
        if selected_analysis 
            
            MR_index = find(strcmp(temp_cohort,MR)); % index of MR in temp_cohort
            recording_time(MR_index) = height(data) / fs; % total monitoring time (s)
            tertile = temp_cohort_assignments(MR_index); % get tertile assignment from selected patients
            
 
        
            if ~any( strcmp(temp_cohort,MR) ) % if MR does not exist in temp cohort...
                continue %skip patinet
            end
            
        end

        
        side = cohort.side(  cohort_ind ); % side of injury (L=0, R=1)

    
            
        % define side
         if side == 0
                lower = "LLA_L";
                upper = "ULA_L";
                cox = "COxl";
                nirs_on = "rSO2L"; % stroke side
                nirs_off = "rSO2R";
                
         elseif side == 1
                lower = "LLA_R";
                upper = "ULA_R";
                cox = "COxr";
                nirs_on = "rSO2R"; % stroke side
                nirs_off = "rSO2L";
         else
             
             display(['Injury side not defined for ', MR])
             continue % skip patient
         end
         
        if off_injury % if performing an off-injury analysis, all sides are switched to uninjured side.
            side = ~side;
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

%             win_before_half = win_before(round((height(win_before))/2):end,:);
%             viable_before = (sum(~any(isnan([win_before_half.(lower),win_before_half.(upper),win_before_half.ABP,win_before_half.(cox)]),2)) / fs) / (height(win_before_half)/fs);
%             
             % check data viability (%)
            viable_before = (sum(~any(isnan([win_before.(lower),win_before.(upper),win_before.ABP,win_before.(cox)]),2)) / fs) / (height(win_before)/fs);
            viable_after = (sum(~any(isnan([win_after.(lower),win_after.(upper),win_after.ABP,win_after.(cox)]),2)) / fs) / (height(win_after)/fs); 
            
           
            if viable_before < viable_threshold || viable_after < viable_threshold % if >90% of data is availabe on both sides

%                 fprintf(['\n',MR, ', Event ID: ', num2str(event_id),', Side: ',num2str(side),', Epoch: ', num2str(epoch),' is not viable for analysis (',num2str(viable_before*100),',',num2str(viable_after*100),')']);

                continue % skip this med admin
            end
            
          % -CALCULATE OVARIATES -----------------------------
      
            
         
            
            cd(old_folder) % to access needed .m files
            
    
            
        
            
            [PERCENT_BELOW_BEFORE, PERCENT_ABOVE_BEFORE] = outside_limits(win_before,side,epochs(epoch)*60);
            [PERCENT_BELOW_AFTER, PERCENT_ABOVE_AFTER,TIME_BELOW, MEAN_NIRS_ON_BELOW , MEAN_NIRS_ON_NOTBELOW,...
            MEAN_NIRS_OFF_BELOW,MEAN_NIRS_OFF_NOTBELOW] = outside_limits(win_after,side,epochs(epoch)*60);
       
            MEAN_NIRS_ON_LA_EFFECT = MEAN_NIRS_ON_BELOW - MEAN_NIRS_ON_NOTBELOW;
            MEAN_NIRS_OFF_LA_EFFECT = MEAN_NIRS_OFF_BELOW - MEAN_NIRS_OFF_NOTBELOW;
            
            if selected_analysis
                
               if epoch == 3 
                 
                   if PERCENT_BELOW_AFTER >= below_threshold % if they spent more than 50% of their time after below LLA....
                       MR_index = find(strcmp(temp_cohort,MR));
                       num_responses(MR_index) =  num_responses(MR_index) + 1; % then classify this admin as a response

                   end

                   MR_index = find(strcmp(temp_cohort,MR)); % index of MR in temp_cohort
                   num_events(MR_index) = num_events(MR_index) + 1; % iteratively add 
                   times_below{MR_index} = [times_below{MR_index} TIME_BELOW]; % iteratively add 
                   viable_recording_time(MR_index) = sum(~any(isnan([data.(lower),data.(upper),data.ABP,data.(cox)]),2)) / fs;
                   [~,~,time_below] = outside_limits(data,side,epochs(epoch)*60); overall_time_below(MR_index) = time_below; % apply outside_limits to the entire data (not just a window)
               end
                
               
               if tertile == 1
                   
                   
                    NIRS_off_series_1 = vertcat(NIRS_off_series_1,[win_before.(nirs_off)' win_after.(nirs_off)']);
                    NIRS_on_series_1 = vertcat(NIRS_on_series_1,[win_before.(nirs_on)' win_after.(nirs_on)']);
                    ABP_series_1 = vertcat(ABP_series_1,[win_before.ABP' win_after.ABP']);
                   
               elseif tertile == 2
                   
                    NIRS_off_series_2 = vertcat(NIRS_off_series_2,[win_before.(nirs_off)' win_after.(nirs_off)']);
                    NIRS_on_series_2 = vertcat(NIRS_on_series_2,[win_before.(nirs_on)' win_after.(nirs_on)']);
                    ABP_series_2 = vertcat(ABP_series_2,[win_before.ABP' win_after.ABP']);
                   
                   
               else
                   
                    NIRS_off_series_3 = vertcat(NIRS_off_series_3,[win_before.(nirs_off)' win_after.(nirs_off)']);
                    NIRS_on_series_3 = vertcat(NIRS_on_series_3,[win_before.(nirs_on)' win_after.(nirs_on)']);
                    ABP_series_3 = vertcat(ABP_series_3,[win_before.ABP' win_after.ABP']);
                   
                   
               end
               
               
            end
           
         
           
           
            cd(newdir)

            MEAN_ABP_BEFORE = mean(win_before.ABP(~isnan(win_before.ABP)));
            MEDIAN_ABP_BEFORE = median((win_before.ABP(~isnan(win_before.ABP))));
            MAX_ABP_BEFORE = max((win_before.ABP(~isnan(win_before.ABP))));
            MIN_ABP_BEFORE = min((win_before.ABP(~isnan(win_before.ABP))));
            
            MEAN_ABP_AFTER = mean(win_after.ABP(~isnan(win_after.ABP)));
            MEDIAN_ABP_AFTER = median(win_after.ABP(~isnan(win_after.ABP)));
            MAX_ABP_AFTER = max(win_after.ABP(~isnan(win_after.ABP)));
            MIN_ABP_AFTER = min(win_after.ABP(~isnan(win_after.ABP)));

            MEAN_COX_BEFORE = mean(win_before.(cox)(~isnan(win_before.(cox))));
            MEAN_COX_AFTER = mean(win_after.(cox)(~isnan(win_after.(cox))));
            
            MEAN_PRX_BEFORE = mean(win_before.PRx(~isnan(win_before.PRx)));
            MEAN_PRX_AFTER = mean(win_after.PRx(~isnan(win_after.PRx)));
            
            MEAN_NIRS_ON_BEFORE = mean(win_before.(nirs_on),'omitnan');
            MEAN_NIRS_ON_AFTER = mean(win_after.(nirs_on),'omitnan');
            NIRS_ON_EFFECT = MEAN_NIRS_ON_AFTER - MEAN_NIRS_ON_BEFORE;
            
            MEAN_NIRS_OFF_BEFORE = mean(win_before.(nirs_off),'omitnan');
            MEAN_NIRS_OFF_AFTER = mean(win_after.(nirs_off),'omitnan');
            NIRS_OFF_EFFECT = MEAN_NIRS_OFF_AFTER - MEAN_NIRS_OFF_BEFORE;
              
            if any(contains(data.Properties.VariableNames,'MAPopt_PRx')) && any(~isnan(data.MAPopt_PRx)) %  if non-NaN PRx values exist..
                
                MEAN_MAPopt_BEFORE = mean(win_before.('MAPopt_PRx')(~isnan(win_before.('MAPopt_PRx'))));
                MEAN_MAPopt_AFTER = mean(win_after.('MAPopt_PRx')(~isnan(win_after.('MAPopt_PRx'))));
                
            elseif side == 0 % use MAPopt_flex on the left side
                   
                MEAN_MAPopt_BEFORE = mean(win_before.('MAPopt_flex_L')(~isnan(win_before.('MAPopt_flex_L'))));
                MEAN_MAPopt_AFTER = mean(win_after.('MAPopt_flex_L')(~isnan(win_after.('MAPopt_flex_L'))));
           
            else  % use MAPopt_flex on the right side
                
                if any(contains(data.Properties.VariableNames,'MAPopt_flex_R_mmHg_'))
                    
                    MEAN_MAPopt_BEFORE = mean(win_before.('MAPopt_flex_R_mmHg_')(~isnan(win_before.('MAPopt_flex_R_mmHg_'))));
                    MEAN_MAPopt_AFTER = mean(win_after.('MAPopt_flex_R_mmHg_')(~isnan(win_after.('MAPopt_flex_R_mmHg_'))));
                   
                    
                    
                else
                    
                    MEAN_MAPopt_BEFORE = mean(win_before.('MAPopt_flex_R')(~isnan(win_before.('MAPopt_flex_R'))));
                    MEAN_MAPopt_AFTER = mean(win_after.('MAPopt_flex_R')(~isnan(win_after.('MAPopt_flex_R'))));
                
                    
                end
                

            end
            
           delta_MAP_AFTER =  MEAN_ABP_AFTER -   MEAN_MAPopt_AFTER; % postive values indicate ABP is too high. 
           delta_MAP_BEFORE = MEAN_ABP_BEFORE -   MEAN_MAPopt_BEFORE;
          
           delta_COX =   MEAN_COX_AFTER -  MEAN_COX_BEFORE; % positive values indicate worsening in autoregulatory function
           delta_PRX = MEAN_PRX_AFTER - MEAN_PRX_BEFORE;
           
           % PACKAGE DATA
           
           covar_data = {event_id, MR, med, dose, PERCENT_BELOW_BEFORE*100, PERCENT_ABOVE_BEFORE*100, PERCENT_BELOW_AFTER*100, ...
                    PERCENT_ABOVE_AFTER*100, MEAN_ABP_BEFORE,  MEAN_ABP_AFTER,  MEAN_COX_BEFORE, ...
                    MEAN_COX_AFTER,  MEAN_MAPopt_BEFORE, MEAN_MAPopt_AFTER, delta_MAP_AFTER,delta_MAP_BEFORE,MEAN_PRX_BEFORE, MEAN_PRX_AFTER,...
                    MEAN_NIRS_ON_BEFORE, MEAN_NIRS_ON_AFTER,NIRS_ON_EFFECT,...
                    MEAN_NIRS_OFF_BEFORE, MEAN_NIRS_OFF_AFTER,NIRS_OFF_EFFECT,...
                    MEAN_NIRS_ON_BELOW,MEAN_NIRS_ON_NOTBELOW,MEAN_NIRS_ON_LA_EFFECT,...
                    MEAN_NIRS_OFF_BELOW,MEAN_NIRS_OFF_NOTBELOW,MEAN_NIRS_OFF_LA_EFFECT,...
                    time_taken};
                    
                    
         
                
                
         
                
           
           epoch_file = [epoch_file ; covar_data]; % vertcat this admin onto the epoch file
           
           %% DR. KIM CODE, ADDED 9/28
           if epoch == 3 && PERCENT_BELOW_AFTER == 0 && (MEAN_ABP_AFTER - MEAN_ABP_BEFORE) > 0
               
               
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
                plot(time_mod,concat_data.ABP,'black','LineWidth',1.5); hold on;
                plot(time_mod,concat_data.(lower),'Color',[0 0.4470 0.7410]);
                plot(time_mod,concat_data.(upper),'Color',[0 0.4470 0.7410]);
                xticks(x_ticks)
                
                xlabel('Time (min)'); ylabel('MAP (mmHg)'); title(['Event: ', num2str(event_id),', Time taken: ',datestr(time_taken)])
               % saveas(fig,[num2str(event_id),'.jpg'])
     
                    
                time_mod_table = table(time_mod);  time_mod_table.Properties.VariableNames = {['time']};
                abp = table(concat_data.ABP); abp.Properties.VariableNames = {['abp']};
                lla = table(concat_data.(lower)); lla.Properties.VariableNames = {['lla']};
                ula = table(concat_data.(upper)); ula.Properties.VariableNames = {['ula']};
                    
                master_table = [time_mod_table,abp,lla,ula];
               % writetable(master_table,[num2str(event_id),'.csv'])
                
           end
          
           
        end
        
  
    end
    
    
    writetable(epoch_file,['Intervention_Analysis_Epoch_', num2str(epochs(epoch)),'.csv'])
    writetable(epoch_file_kim,['Intervention_Analysis_KIM_Epoch_', num2str(epochs(epoch)),'.csv']) % DR. KIM CODE, ADDED 9/28
    display(["Just packaged Epoch ", num2str(epochs(epoch))])
    cd(analysis_folder)
    
end


cd(old_folder)



%%

%%

times_below_mean = cellfun(@mean,times_below);

times_below_std = cellfun(@std,times_below);

NIRS_on_series_1_mean = mean(NIRS_on_series_1,1,'omitnan')';
NIRS_on_series_2_mean = mean(NIRS_on_series_2,1,'omitnan')';
NIRS_on_series_3_mean = mean(NIRS_on_series_3,1,'omitnan')';

NIRS_off_series_1_mean = mean(NIRS_off_series_1,1,'omitnan')';
NIRS_off_series_2_mean = mean(NIRS_off_series_2,1,'omitnan')';
NIRS_off_series_3_mean = mean(NIRS_off_series_3,1,'omitnan')';

NIRS_on_series_1_std = std(NIRS_on_series_1,1,'omitnan')';
NIRS_on_series_2_std = std(NIRS_on_series_2,1,'omitnan')';
NIRS_on_series_3_std = std(NIRS_on_series_3,1,'omitnan')';

NIRS_off_series_1_std =  std(NIRS_off_series_1,1,'omitnan')';
NIRS_off_series_2_std =  std(NIRS_off_series_2,1,'omitnan')';
NIRS_off_series_3_std =  std(NIRS_off_series_3,1,'omitnan')';

ABP_1 = std(ABP_series_1,1,'omitnan')';
ABP_2 = std(ABP_series_2,1,'omitnan')';
ABP_3 = std(ABP_series_3,1,'omitnan')';


