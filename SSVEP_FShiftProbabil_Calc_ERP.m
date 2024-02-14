%% parameters
clearvars
F.Pathlocal             = 'E:\work\data\SSVEP_FShift_Probabil\';
F.Pathlocal             = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2023_FShift_Probabil\';

F.PathInEEG             = fullfile(F.Pathlocal, 'eeg\epoch_erp\');
F.PathInBehavior        = fullfile(F.Pathlocal, 'behavior\raw\');
F.PathInSCADS           = fullfile(F.Pathlocal, 'eeg\SCADS_erp\');
F.PathOut               = fullfile(F.Pathlocal, 'eeg\erp\'); 
F.subjects              = arrayfun(@(x) sprintf('%02.0f',x),1:30,'UniformOutput',false)';
F.sub2use               = [1 3 4 5 6 7 9 10 11 12 13 14 15 18 20 21 22 23 24 25];%:53;
F.sub2use               = [1 3 4 5 6 7 9 10 11 13 15 18 20 21 22 23 24 25];%:53; % for subject 12, 14: eeg and behavior data don't match
% F.sub2use               = [15 18 20 21 22 23 24 25];%:53;

F.trigger               = {[111] [112] [121] [122] [131] [132] [211] [212] [221] [222] [231] [232]}; % regular
% F.trigger               = {[100];[200]};

F.EEGChans              = 64;
F.ERPepoch              = [-0.5 1];
F.ERPbase               = [-0.1 0]; % erp baseline in s
F.CSD_flag              = 1; % 0 = no; 1 = yes
F.ERP_FiltFreq          = [0 15];

%F.TFAfreqs              = [5:(1/6):40];
F.con1name              = 'validity';
F.con1label             = {'valid';'valid';'invalid;';'invalid';'neutral';'neutral'; 'valid';'valid';'invalid;';'invalid';'neutral';'neutral'};
F.con2name              = 'target position';
F.con2label             = {'RDK1';'RDK1';'RDK1';'RDK1';'RDK1';'RDK1';'RDK2';'RDK2';'RDK2';'RDK2';'RDK2';'RDK2'};
F.con3name              = 'event type';
F.con2label             = {'chroma-';'chroma+';'chroma-';'chroma+';'chroma-';'chroma+';'chroma-';'chroma+';'chroma-';'chroma+';'chroma-';'chroma+'};
F.targtype              = {'chroma-';'chroma+'};

%% start processing
%% loop across subjects
for i_sub = 1:numel(F.sub2use)
    %% load files
    % EEG
    EEG = pop_loadset('filename',sprintf('VP%s_e.set',F.subjects{F.sub2use(i_sub)}),'filepath',F.PathInEEG);
    prep_input = load(fullfile(F.PathInSCADS,sprintf('VP%s_Preprocess_summary.mat',F.subjects{F.sub2use(i_sub)})));
    % pop_eegplot(EEG,1,1,1)
    % behavior (loads latest file)
    t.files = dir(fullfile(F.PathInBehavior,sprintf('VP%s_timing*.mat',F.subjects{F.sub2use(i_sub)})));
    [t.val t.idx ]=max([t.files.datenum]);
    behavior = load(fullfile(F.PathInBehavior,t.files(t.idx).name));
    
         
    %% do csd transform
    if F.CSD_flag==1
        if  i_sub == 1 % calculate CSD matrix
%             CSD.chanmat=ExtractMontage('C:\Users\psy05cvd\Dropbox\work\matlab\software\toolboxes\CSD\resource\10-5-System_Mastoids_EGI129.csd',{EEG.chanlocs.labels}');
            CSD.chanmat=ExtractMontage('C:\Users\EEG\Documents\MATLAB\christopher\general_functions\CSD\resource\10-5-System_Mastoids_EGI129.csd',{EEG.chanlocs.labels}');
            [CSD.G,CSD.H] = GetGH(CSD.chanmat);
        end
        fprintf(1,'\n###\ncalculating CSD transform\n###\n')
        for i_tr = 1:EEG.trials
            % csd of raw data
            EEG.data(:,:,i_tr)= CSDTransform(EEG.data(:,:,i_tr), CSD.G, CSD.H);
        end
    end
    % pop_eegplot(EEG,1,1,1)

    %% additional calculation and conditiona allocation
    % raw (for FFT)
    [EEG_ep, t.indices] = pop_epoch(EEG, num2cell(unique(cell2mat(F.trigger))), F.ERPepoch, 'epochinfo', 'yes');
    EEG_ep = pop_rmbase(EEG_ep,F.ERPbase.*1000);
    % filtered for ERP
    EEG_f = pop_eegfiltnew(EEG, F.ERP_FiltFreq(1), F.ERP_FiltFreq(2), 8*EEG.srate, 0, [], 0);
    [EEG_fep, t.indices] = pop_epoch(EEG_f, num2cell(unique(cell2mat(F.trigger))), F.ERPepoch, 'epochinfo', 'yes');
    EEG_fep = pop_rmbase(EEG_fep,F.ERPbase.*1000);

    t.behavior = horzcat(behavior.resp.experiment{:});
    t.prep_idx = prep_input.PreProc.trial_blink & prep_input.PreProc.trial_eyemov & prep_input.PreProc.trial_SCADS;
    t.behavior = t.behavior(t.prep_idx);


    % add some information
    t.ur_epoch = num2cell(find(t.prep_idx));
    [t.behavior.urepoch] = deal(t.ur_epoch{:});


    for i_event = 1:numel(t.behavior)
        % onset times after cue
        t.behavior(i_event).postcue_onset = ...
            (t.behavior(i_event).event_onset_times-t.behavior(i_event).pre_cue_times)*1000;
        % event color label
        t.behavior(i_event).eventRDK_col_label = behavior.RDK.RDK(t.behavior(i_event).eventpos).col_label;

        t.behavior(i_event).evnt_type_label = F.targtype{t.behavior(i_event).event_type};
    end

    %% some bookkeeping
    EP.behavior = t.behavior;
    EP.EEG_ep = EEG_ep;
    EP.EEG_fep = EEG_fep;
    EP.parameters = F;
    EP.RDK = behavior.RDK.RDK;
    
    %% save
    EP.savetime=datestr(now);
    if ~exist(F.PathOut); mkdir(F.PathOut); end
    fprintf(1,'|| saving file ||  %s\\VP%s_target_erp.mat ||\n', ...
        F.PathOut,F.subjects{F.sub2use(i_sub)})
    save(fullfile(F.PathOut, sprintf('VP%s_target_erp.mat',F.subjects{F.sub2use(i_sub)})), 'EP')
    
end