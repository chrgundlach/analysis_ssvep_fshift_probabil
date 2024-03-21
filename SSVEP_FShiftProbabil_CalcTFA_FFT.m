%% parameters
clearvars
F.Pathlocal             = 'E:\work\data\SSVEP_FShift_Probabil\';
F.Pathlocal             = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2023_FShift_Probabil\';

% F.PathInEEG             = fullfile(F.Pathlocal, 'eeg\epoch\');
F.PathInEEG             = fullfile(F.Pathlocal, 'eeg\epoch_2stfa\');
F.PathInBehavior        = fullfile(F.Pathlocal, 'behavior\raw\');
% F.PathInSCADS           = fullfile(F.Pathlocal, 'eeg\SCADS\');
F.PathInSCADS           = fullfile(F.Pathlocal, 'eeg\SCADS_2stfa\');
F.PathOut               = fullfile(F.Pathlocal, 'eeg\tfa_2s\'); % with FWHM 0.5
F.subjects              = arrayfun(@(x) sprintf('%02.0f',x),1:50,'UniformOutput',false)';
% F.sub2use               = [1 3 4 5 6 7 9 10 11 12 13 14 15 18 20:36];%:53;
F.sub2use               = [35:36];%:53;

F.trigger               = {[10] [20] [30] [40] [50] [60]}; % regular
% F.trigger               = {[100];[200]};

F.EEGChans              = 64;

%F.TFAfreqs              = [5:(1/6):40];
F.TFAfreqs              = [4:0.25:40];
F.TFAFlag               = [2]; % 1=wavelet; 2=gabor
F.TFAepoch              = [-1.5 2.5];
F.CSD_flag              = 1; % 0 = no; 1 = yes

F.FFT_timewins          = {[-1 0]; [0 1]; [0.5 1.5]; [1 2]};
F.FFT_freqres           = 16384;

F.conds=            {'valid | target RDK1'; 'valid | target RDK2'; ...
    'invalid | target RDK2'; 'invalid | target RDK2';...
    'neutral | target RDK1'; 'neutral | target RDK2'};
F.conds2=            {'valid | attend RDK1'; 'valid | attend RDK2'; ...
    'invalid | attend RDK2'; 'invalid | attend RDK1';...
    'neutral | attend RDK1'; 'neutral | attend RDK2'};
%% start processing
%% loop across subjects
for i_sub = 1:numel(F.sub2use)
    %% load files
    % EEG
    EEG = pop_loadset('filename',sprintf('VP%s_e.set',F.subjects{F.sub2use(i_sub)}),'filepath',F.PathInEEG);
    Preprocessing = load(fullfile(F.PathInSCADS,sprintf('VP%s_Preprocess_summary.mat',F.subjects{F.sub2use(i_sub)})));
    % pop_eegplot(EEG,1,1,1)
    % behavior (loads latest file)
    t.files = dir(fullfile(F.PathInBehavior,sprintf('VP%s_timing*.mat',F.subjects{F.sub2use(i_sub)})));
    [t.val t.idx ]=max([t.files.datenum]);
    behavior = load(fullfile(F.PathInBehavior,t.files(t.idx).name));
    
     % select certain epochs
    t.idx = [EEG.event(ismember([EEG.event.type],unique(cell2mat(F.trigger)))).epoch];
    EEG = pop_select(EEG,'trial',t.idx);
    
    % select certain time window
    EEG = pop_select(EEG, 'time', F.TFAepoch);
    
    %% do csd transform
    if F.CSD_flag==1
        if  i_sub == 1 % calculate CSD matrix
            CSD.chanmat=ExtractMontage('C:\Users\psy05cvd\Dropbox\work\matlab\software\toolboxes\CSD\resource\10-5-System_Mastoids_EGI129.csd',{EEG.chanlocs.labels}');
%             CSD.chanmat=ExtractMontage('C:\Users\EEG\Documents\MATLAB\christopher\general_functions\CSD\resource\10-5-System_Mastoids_EGI129.csd',{EEG.chanlocs.labels}');
            [CSD.G,CSD.H] = GetGH(CSD.chanmat);
        end
        fprintf(1,'\n###\ncalculating CSD transform\n###\n')
        for i_tr = 1:EEG.trials
            % csd of raw data
            EEG.data(:,:,i_tr)= CSDTransform(EEG.data(:,:,i_tr), CSD.G, CSD.H);
        end
    end
    % pop_eegplot(EEG,1,1,1)
    
    %% calculate induced TFA (average in frequency domain) and evoked TFA
    TFA.electrodes = EEG.chanlocs;
    TFA.RDK = behavior.RDK;
    
    TFA.params.srate = EEG.srate;
    
    %filter
    TFA.params.filter = {[] []};
    
    TFA.data_evo=nan(numel(F.TFAfreqs),ceil(EEG.pnts/2), F.EEGChans,numel(F.trigger),'single');
    TFA.data_ind=nan(numel(F.TFAfreqs),ceil(EEG.pnts/2), F.EEGChans,numel(F.trigger),'single');
    TFA.all_trialnum = EEG.trials;
    TFA.con_trigger = F.trigger;
    TFA.con_trialnum = cellfun(@(x) sum(ismember([EEG.event.type],x)), F.trigger);
    
    fprintf(1,'calculating TFA with %1.0f frequency steps:\n',numel(F.TFAfreqs))
    TFA.type = 'gabor';
    TFA.f = F.TFAfreqs;
    %         TFA.params.gabor_FWHM_freq = 0.5;
    TFA.params.gabor_FWHM_freq = 1;
    
    % index trials
    t.trialindex = cellfun(@(x) [EEG.event([ismember([EEG.event.type],x)]).epoch], F.trigger,'UniformOutput',false);
    
    %% loop across frequencies
    for i_freq = 1:numel(TFA.f)
        %% calculate gabor
        fprintf(1,'%1.0f of %1.0f\n',i_freq, numel(F.TFAfreqs))
        % induced TFA
        EEG_Gabor = eegF_Gabor(EEG, TFA.f(i_freq),TFA.params.gabor_FWHM_freq);
        EEG_m = pop_select(EEG,'trial',1);
        for i_con = 1:numel(TFA.con_trigger)
            % induced TFA
            TFA.data_ind(i_freq,:,:,i_con)=single(imresize(mean(EEG_Gabor.data(:,:,t.trialindex{i_con}),3),...
                [EEG_Gabor.nbchan,ceil(EEG_Gabor.pnts/2)])');
            
            % evoked TFA
            EEG_m.data = mean(EEG.data(:,:,t.trialindex{i_con}),3);
            EEG_m_Gabor = eegF_Gabor(EEG_m, TFA.f(i_freq),TFA.params.gabor_FWHM_freq);
            TFA.data_evo(i_freq,:,:,i_con) = single(imresize(EEG_m_Gabor.data(:,:),...
                [EEG_Gabor.nbchan,ceil(EEG_Gabor.pnts/2)])');
        end
    end
    TFA.t=imresize(EEG_Gabor.times,[1,ceil(EEG_Gabor.pnts/2)]);
    fprintf(1,'...done\n')
    
%     % plotting for checking
%     pl.elec2plot = {'Oz'};
%     %             pl.elec2plot = {'C3';'CP3'};
%     pl.base=[-500 0]; pl.base_i=dsearchn(TFA.t',pl.base');
%     pl.t2pl = [-1000 2000]; pl.t2pl_i = dsearchn(TFA.t',pl.t2pl');
%     pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
%     
%     pl.data = mean(TFA.data_ind(:,pl.t2pl_i(1):pl.t2pl_i(2),pl.elec2plot_i,:),3); pl.clims = [0 max(pl.data(:))];
%     pl.data = mean(TFA.data_evo(:,pl.t2pl_i(1):pl.t2pl_i(2),pl.elec2plot_i,:),3); pl.clims = [0 max(pl.data(:))];
%     pl.ti1 = pl.elec2plot; pl.ti2 = 'raw';  pl.clims = [0 1]*max(abs(pl.data(:)));
%     %         pl.data = mean(bsxfun(@minus, TFA.data(:,pl.t2pl_i(1):pl.t2pl_i(2),pl.elec2plot_i,:), mean(TFA.data(:,pl.base_i(1):pl.base_i(2),pl.elec2plot_i,:),2)),3);
%     %         pl.ti1 = pl.elec2plot; pl.ti2 = 'baseline corrected';  pl.clims = [-1 1]*max(abs(pl.data(:)));
%     
%     
%     figure;
%     for i_con = 1:size(pl.data,4)
%         subplot(size(pl.data,4)+1,1,i_con)
%         imagesc(TFA.t(pl.t2pl_i(1):pl.t2pl_i(2)),TFA.f,pl.data(:,:,:,i_con),pl.clims)
%         colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
%         set(gca,'YDir','normal')
%         title(sprintf('%s tfa %s attend %1.0f Hz', pl.ti2, vararg2str(pl.ti1), TFA.RDK.RDK(i_con).freq), 'FontSize',8,'Interpreter','none')
%         colorbar
%         xlabel('time in ms')
%         ylabel('frequency in Hz')
%         set(gca,'FontSize',8)
%     end
%     subplot(size(pl.data,4)+1,1,i_con+1)
%     imagesc(TFA.t(pl.t2pl_i(1):pl.t2pl_i(2)),TFA.f,pl.data(:,:,:,1)-pl.data(:,:,:,2),[-1 1]*max(max(abs(pl.data(:,:,:,1)-pl.data(:,:,:,2)))))
%     colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
%     set(gca,'YDir','normal')
%     title(sprintf('%s tfa %s attend %1.0f - attend %1.0f Hz', pl.ti2, vararg2str(pl.ti1), TFA.RDK.RDK(1:2).freq), 'FontSize',8,'Interpreter','none')
%     colorbar
%     xlabel('time in ms')
%     ylabel('frequency in Hz')
%     set(gca,'FontSize',8)
    
    
    %% save FFT transform
    TFA.FFT.res = F.FFT_freqres;
    TFA.FFT.timewin = F.FFT_timewins;
    TFA.FFT.data_ind=nan(TFA.FFT.res, F.EEGChans,numel(F.trigger),numel(TFA.FFT.timewin));
    TFA.FFT.data_evo= TFA.FFT.data_ind;
    fprintf('FFT: %1.0f windows for %1.0f electrodes:\n', numel(TFA.FFT.timewin), EEG.nbchan)
    for i_win = 1:numel(TFA.FFT.timewin)
        % select data
        EEGt = pop_select(EEG,'time',TFA.FFT.timewin{i_win});
        % surrogate mean data for evoked spectra
        EEGt_m = pop_select(EEGt,'trial',1);
        % detrend
        EEGt = eegF_Detrend(EEGt,[]);
        % loop across channels
        fprintf('channel:  ')
        t.FFTdata = nan(TFA.FFT.res, F.EEGChans,EEGt.trials);
        for i_el = 1:EEGt.nbchan
            fprintf('\b\b%02.0f',i_el)
            t.FFTdata(:,i_el,:) = squeeze(abs(fft(EEGt.data(i_el,:,:),TFA.FFT.res,2))*2/size(EEGt.data,2));
        end
        
        % loop across conditions
        for i_con = 1:numel(TFA.con_trigger)
            % induced
            TFA.FFT.data_ind(:,:,i_con,i_win)=mean(t.FFTdata(:,:,t.trialindex{i_con}),3);
            EEGt_m.data = detrend(mean(EEGt.data(:,:,t.trialindex{i_con}),3)')';
            TFA.FFT.data_evo(:,:,i_con,i_win) = squeeze(abs(fft(EEGt_m.data,TFA.FFT.res,2))*2/size(EEGt.data,2))';
        end
        TFA.FFT.freqs = ((0:size(TFA.FFT.data_ind,1)-1)/size(TFA.FFT.data_ind,1)) * EEGt.srate;
        fprintf('\n')
    end
    
    
%     %     plotting for checking
% %     pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'};
%     pl.elec2plot = {'Oz'};
%     pl.elec2plot = {'POz';'O1';'Oz';'O2';'Iz'};
%     pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
%     
%     pl.data = mean(TFA.FFT.data_ind(:,pl.elec2plot_i,:,:),2); pl.ylims = [0 max(pl.data(:))];
%     pl.data = mean(TFA.FFT.data_evo(:,pl.elec2plot_i,:,:),2); pl.ylims = [0 max(pl.data(:))];
%     figure;
%     for i_con = 1:size(pl.data,3)
%         subplot(size(pl.data,3),1,i_con)
%         plot(TFA.FFT.freqs,squeeze(pl.data(:,:,i_con,:)))
%         xlabel('frequency in Hz')
%         ylabel('amplitude in \muV')
%         xlim([0 50])
%         ylim(pl.ylims)
%         title(sprintf('%s | RDKs [1 2 3]: [%1.0f %1.0f %1.0f Hz]',F.conds{i_con},TFA.RDK.RDK.freq))
%         legend(cellfun(@(x) num2str(x),TFA.FFT.timewin,'UniformOutput',false))
%         set(gca,'FontSize',8)
%     end
%     figure;
%     for i_time = 1:size(pl.data,4)
%         subplot(size(pl.data,4),1,i_time)
%         plot(TFA.FFT.freqs,squeeze(pl.data(:,:,:,i_time)))
%         xlabel('frequency in Hz')
%         ylabel('amplitude in \muV')
%         xlim([0 50])
%         ylim(pl.ylims)
%         legend(F.conds)
%         title(sprintf('time [%1.0f %1.0f]s | RDKs [1 2 3]: [%1.0f %1.0f %1.0f Hz]',TFA.FFT.timewin{i_time}, TFA.RDK.RDK.freq))
%         set(gca,'FontSize',8)
%     end
%     
    %% save
    TFA.savetime=datestr(now);
    if ~exist(F.PathOut); mkdir(F.PathOut); end
    fprintf(1,'|| saving file ||  %s\\VP%s_tfa.mat ||\n', ...
        F.PathOut,F.subjects{F.sub2use(i_sub)})
    save(fullfile(F.PathOut, sprintf('VP%s_tfa.mat',F.subjects{F.sub2use(i_sub)})), 'TFA')
    
end