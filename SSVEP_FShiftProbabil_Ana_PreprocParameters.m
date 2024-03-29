%% script to analyse preprocessing parameters

clearvars

%% parameters
% general parameters
F.PathIn            = 'E:\work\data\SSVEP_FShift_Probabil\eeg\SCADS\';
F.PathIn            = 'N:\AllgPsy\experimental_data\2023_FShift_Probabil\eeg\SCADS\';
F.PathIn            = 'N:\AllgPsy\experimental_data\2023_FShift_Probabil\eeg\SCADS_2stfa\';
F.subjects          = cellfun(@(x) sprintf('%02.0f',x),num2cell(1:40),'UniformOutput', false)';
% F.subs2use          = [9 10 11 12];%
F.subs2use          = [1 3:6 7 9 10 11 12 13 14 15 18 20:36];
F.con2an            = {[10] [20] [30] [40] [50] [60]};
F.con2an_label      = {'valid | target RDK1'; 'valid | target RDK2'; ...
    'invalid | target RDK2'; 'invalid | target RDK2';...
    'neutral | target RDK1'; 'neutral | target RDK2'};
% F.con2an_label      = {'valid | attend RDK1'; 'valid | attend RDK2'; ...
%     'invalid | attend RDK1'; 'invalid | attend RDK2';...
%     'neutral | attend both'; 'neutral | attend both'};
F.con2an            = {[10 20] [30 40] [50 60]};
F.con2an_label      = {'valid'; 'invalid'; 'neutral'};

%% loop across subjects to extract
for i_sub = 1:numel(F.subs2use)
    % load mat
    input = load(sprintf('%sVP%s_Preprocess_summary.mat',F.PathIn,F.subjects{F.subs2use(i_sub)}));
    % extract relevant data
    params.SCADS_InterpChansPerTrial(i_sub)=input.SumData.interpolated_channels_avgPERtrial;
    params.trialnr_rem_blinks(i_sub) = numel(find(input.PreProc.trial_blink(ismember(input.PreProc.trial_con,unique(cell2mat(F.con2an))))==0));
    params.trialnr_rem_eyemov(i_sub) = numel(find(input.PreProc.trial_eyemov(ismember(input.PreProc.trial_con,unique(cell2mat(F.con2an))))==0));
    params.trialnr_rem_SCADS(i_sub) = numel(find(input.PreProc.trial_SCADS(ismember(input.PreProc.trial_con,unique(cell2mat(F.con2an))))==0));
    params.trialnr_in_anaylsis(:,i_sub) = cellfun(...
        @(x) numel(find(input.PreProc.trial_blink & input.PreProc.trial_eyemov & input.PreProc.trial_SCADS & ...
        ismember(input.PreProc.trial_con, x))),...
        F.con2an)';
end

%% output
fprintf('interpolated channels per trial: M=%1.3f; Std=%1.3f\n',mean(params.SCADS_InterpChansPerTrial),std(params.SCADS_InterpChansPerTrial))
fprintf('removed trials with blinks: M=%1.3f; Std=%1.3f\n',mean(params.trialnr_rem_blinks),std(params.trialnr_rem_blinks))
fprintf('removed trials with eye movements: M=%1.3f; Std=%1.3f\n',mean(params.trialnr_rem_eyemov),std(params.trialnr_rem_eyemov))
fprintf('removed trials by SCADS: M=%1.3f; Std=%1.3f\n',mean(params.trialnr_rem_SCADS),std(params.trialnr_rem_SCADS))
for i_con = 1:numel(F.con2an_label)
    fprintf('number of trials that entered analysis for condition %s: M=%1.3f; Std=%1.3f\n',...
        F.con2an_label{i_con}, mean(params.trialnr_in_anaylsis(i_con,:),2), std(params.trialnr_in_anaylsis(i_con,:),1,2))
end
fprintf('average trial number per condition M = %1.3f; SD = %1.3f\n',...
     mean(mean(params.trialnr_in_anaylsis,2)), mean(std(params.trialnr_in_anaylsis,1,2)))
fprintf('overall rejection rate of trials M = %1.3f; SD = %1.3f\n',...
    mean(mean(100-params.trialnr_in_anaylsis./120.*100,1)),std(mean(100-params.trialnr_in_anaylsis./120.*100,1)))