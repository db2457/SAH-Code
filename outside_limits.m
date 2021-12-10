function [percent_below,percent_above,time_below,...
      MEAN_COx_BELOW , MEAN_COx_NOTBELOW,...
      MEAN_PRx_BELOW,MEAN_PRx_NOTBELOW] = outside_limits(window,side,epoch)
%GET_EVENTS Summary of this function goes here
%   Return NaN if the window is COMPLETELY NaN
%   epoch: duration of epoch in seconds
% 11/9/2021: added code to compute absolute and percent differences in NIRS
% between times when ABP < LLA, within LA, and > ULA

    fs = 0.1;
    
    if side == 0
        lower = "LLA_L";
        upper = "ULA_L";
        nirs_on = "rSO2L"; % stroke side
        nirs_off = "rSO2R";
        COx = "COxl";
  
    else
        lower = "LLA_R";
        upper = "ULA_R";
        nirs_on = "rSO2R"; % stroke side
        nirs_off = "rSO2L";
        COx = "COxr";
    end
    
    PRx = window.PRx;
    COx = window.(COx);
    
    time = window.DateTime;
    ABP = window.ABP;
    LLA = window.(lower);
    ULA = window.(upper);
    NIRS_on = window.(nirs_on);
    NIRS_off = window.(nirs_off);
%     COx
%     PRx = window.
    
    
    hyper_index = ABP > ULA; % 1 where true
    within_index = ABP <= ULA | ABP >= LLA;
    hypo_index = ABP < LLA; % 1 where true
    
    
   % time_notNaN = sum(~any(isnan([LLA,ULA,ABP]),2)) / fs;  % time of
   % viable recording (s). Not using this as a the denominator anymore! We
   % are using epoch to be more accurate :) the data quality filtering step already took
   % place when we screened for viability
    
    percent_below = ( sum(hypo_index)/fs ) / epoch; time_below = ( sum(hypo_index)/fs );
    
    percent_above = ( sum (hyper_index)/fs ) / epoch;
    
    
    
    % NIRS analysis added 11/9/2021 per Nils request
    
    MEAN_NIRS_ON_BELOW = mean(NIRS_on(hypo_index),'omitnan');
    MEAN_NIRS_ON_WITHIN = mean(NIRS_on(within_index),'omitnan');
    MEAN_NIRS_ON_ABOVE = mean(NIRS_on(hyper_index),'omitnan');
    MEAN_NIRS_ON_NOTBELOW = mean(NIRS_on(~hypo_index),'omitnan');
    
    MEAN_NIRS_OFF_BELOW = mean(NIRS_off(hypo_index),'omitnan');
    MEAN_NIRS_OFF_WITHIN = mean(NIRS_off(within_index),'omitnan');
    MEAN_NIRS_OFF_ABOVE = mean(NIRS_off(hyper_index),'omitnan');
    MEAN_NIRS_OFF_NOTBELOW = mean(NIRS_off(~hypo_index),'omitnan');
    

    % COx and PRx analysis added 12/9/2021 per Nils request
    
    MEAN_COx_BELOW = mean(COx(hypo_index),'omitnan');
    MEAN_COx_WITHIN = mean(COx(within_index),'omitnan');
    MEAN_COx_ABOVE = mean(COx(hyper_index),'omitnan');
    MEAN_COx_NOTBELOW = mean(COx(~hypo_index),'omitnan');
    
     
    MEAN_PRx_BELOW = mean(PRx(hypo_index),'omitnan');
    MEAN_PRx_WITHIN = mean(PRx(within_index),'omitnan');
    MEAN_PRx_ABOVE = mean(PRx(hyper_index),'omitnan');
    MEAN_PRx_NOTBELOW = mean(PRx(~hypo_index),'omitnan');
    
 
  
% %DEBUGGING
%     if isnan( NIRS_OFF_PERCENT_DIFFERENCE)
%        display("Stop") 
%     end
%     
%     if isnan( NIRS_OFF_ABS_DIFFERENCE ) 
%        display("Stop") 
%     end
%                 
%   
end

%DEBUGGING
% plot(window.DateTime,ABP); hold on;
% 
% plot(window.DateTime,LLA)
% 
% plot(window.DateTime,ULA)
% 
% legend('ABP','LLA','ULA')

