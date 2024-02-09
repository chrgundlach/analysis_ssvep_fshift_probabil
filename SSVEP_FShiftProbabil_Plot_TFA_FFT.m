%% plot TFA images
clearvars
% F.PathInEEG             = 'D:\work\data\SSVEP_FShiftAlpha\eeg\TFA'; % with FWHM 1
F.PathInEEG             = 'N:\AllgPsy\experimental_data\2023_FShift_Probabil\eeg\tfa'; % with FWHM 1
F.PathInEEG             = 'C:\Users\psy05cvd\Dropbox\work\projects\SSVEP_FShift_Probabil\EEG_data'; % with FWHM 1


F.Subs                  = arrayfun(@(x) sprintf('%02.0f',x),1:40,'UniformOutput',false)';
F.Subs2use              = [1 3 4 5 6 7 9 10 11 12 13 14 15 18 20 21 22 23 24 25]; 
                        % 2 and 8 are excluded as the didn't do the task properly, sub 11 has potentially low number of trials
F.TFA.baseline          = [-500 -250];



F.conds=            {'valid | target RDK1'; 'valid | target RDK2'; ...
    'invalid | target RDK2'; 'invalid | target RDK2';...
    'neutral | target RDK1'; 'neutral | target RDK2'};
F.conds2=            {'valid | attend RDK1'; 'valid | attend RDK2'; ...
    'invalid | attend RDK2'; 'invalid | attend RDK1';...
    'neutral | attend RDK1'; 'neutral | attend RDK2'};
F.stim_frequency = [18 21 24];
F.stim_colnames = {'green';'redish';'magenta'};
%% load data
for i_sub = 1:numel(F.Subs2use)
    fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%s_tfa.mat ||\n',...
        i_sub,numel(F.Subs2use),F.PathInEEG,F.Subs{F.Subs2use(i_sub)})
    
    temp.tfa = open(fullfile(F.PathInEEG,sprintf('VP%s_tfa.mat',F.Subs{F.Subs2use(i_sub)})));
    
    % convert to single
    temp.tfa.TFA.data_evo = single(temp.tfa.TFA.data_evo);
    temp.tfa.TFA.data_ind = single(temp.tfa.TFA.data_ind);
    temp.tfa.TFA.FFT.data_evo = single(temp.tfa.TFA.FFT.data_evo);
    temp.tfa.TFA.FFT.data_ind = single(temp.tfa.TFA.FFT.data_ind);
    
    
    % preallocate memory
    if i_sub == 1
        TFA.data_evo = single(nan([size(temp.tfa.TFA.data_evo),numel(F.Subs2use)]));
        TFA.data_ind = single(nan([size(temp.tfa.TFA.data_ind),numel(F.Subs2use)]));
        TFA.time = temp.tfa.TFA.t;
        TFA.frequency = temp.tfa.TFA.f;
        TFA.electrodes = temp.tfa.TFA.electrodes;
        TFA.con_trialnum = temp.tfa.TFA.con_trialnum;
        TFA.srate = temp.tfa.TFA.params.srate/2;
        TFA.fftdata_ind = nan([size(temp.tfa.TFA.FFT.data_ind),numel(F.Subs2use)]);
        TFA.fftdata_evo = nan([size(temp.tfa.TFA.FFT.data_evo),numel(F.Subs2use)]);
        TFA.ffttimewin = temp.tfa.TFA.FFT.timewin;
        TFA.fftfreqs = temp.tfa.TFA.FFT.freqs;
    end
    
    % assign data
    TFA.data_evo(:,:,:,:,i_sub) = temp.tfa.TFA.data_evo; % evoked data
    TFA.data_ind(:,:,:,:,i_sub) = temp.tfa.TFA.data_ind; % induced data
    %     TFA(i_exp).data_bc(:,:,:,i_sub) = bsxfun(@minus, temp.tfa.TFA.data_induced, ...
    %         mean(temp.tfa.TFA.data(:,eeg_time2points(TFA.baseline(1),TFA.time):eeg_time2points(TFA.baseline(2),TFA.time),:,:,:),2));
%     TFA(i_exp).data_bc(:,:,:,:,i_sub) = 100*((bsxfun(@rdivide, temp.tfa.TFA.data, ...
%         mean(temp.tfa.TFA.data(:,eeg_time2points(F.TFA.baseline(1),TFA(i_exp).time):eeg_time2points(F.TFA.baseline(2),TFA(i_exp).time),:,:,:),2)))-1);
    TFA.fftdata_evo(:,:,:,:,i_sub) = temp.tfa.TFA.FFT.data_evo;
    TFA.fftdata_ind(:,:,:,:,i_sub) = temp.tfa.TFA.FFT.data_ind;
    TFA.RDK(i_sub) = temp.tfa.TFA.RDK;
    TFA.con_trialnum(i_sub,:) = temp.tfa.TFA.con_trialnum;
    
    clear temp
    
end

pl.plegend = {'1';'0.5';'0.25';'0.1';'0.05';'0.01';'0.001';'0.0001';'0.00001'};
pl.pcorrect = [0 abs(log10([0.5 0.25 0.1 0.05 0.01 0.001 0.0001 0.00001]))];



%% check for descriptives
% frequencies
diag.freqs = cell2mat(cellfun(@(x) [x.freq],{TFA.RDK.RDK},'UniformOutput',false)');
diag.freqs_frequency = [F.stim_frequency' cell2mat(arrayfun(@(x) sum(diag.freqs==x), F.stim_frequency ,'UniformOutput',false)')];
diag.cols = cellfun(@(x) {x.col_label},{TFA.RDK.RDK},'UniformOutput',false)';
diag.cols = vertcat(diag.cols{:});
diag.cols_frequency = [F.stim_colnames(1:3) num2cell(cell2mat(cellfun(@(x) sum(strcmp(diag.cols,x)), F.stim_colnames(1:3) ,'UniformOutput',false)))];


%% plot grand mean FFT data | spectra | with focus on specific frequency
% plotting parameters
pl.elec2plot = {'O1';'Oz';'O2';'I1';'Iz';'I2'};
% pl.elec2plot_cluster = {{TFA.electrodes(33:64).labels}';{TFA.electrodes(1:32).labels}';{TFA.electrodes(1:32).labels}'};
pl.elec2plot_label = {'central'}; % which stimulus is best captured?

% pl.elec2plot = {'O1';'Oz';'O2';'Iz'};
% % pl.elec2plot_cluster = {{TFA.electrodes(33:64).labels}';{TFA.electrodes(1:32).labels}';{TFA.electrodes(1:32).labels}'};
% pl.elec2plot_label = {'central_small'}; % which stimulus is best captured?
% 
% pl.elec2plot = {'P5';'P7';'PO7';'PO3'; 'O1';'Oz';'O2';'I1';'Iz';'I2'; 'PO4'; 'P8';'P6';'PO8'};
% pl.elec2plot_label = {'large'}; % which stimulus is best captured?

% plotting parameters
pl.elec2plot = {'POz';'O1';'Oz';'O2';'Iz'};
% pl.elec2plot_cluster = {{TFA.electrodes(33:64).labels}';{TFA.electrodes(1:32).labels}';{TFA.electrodes(1:32).labels}'};
pl.elec2plot_label = {'central'}; % which stimulus is best captured?

% cluster analysis
% pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) startsWith({TFA(1).electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

% plot channels with labels
% findchannellabel( 1:64 , 1)

pl.time2plot = [2]; % [1 2 4]
pl.sub2plot = [1:numel(F.Subs2use)]; % sub 15 ith no SSVEP default = 1:numel(F.Subs2use)

% extract data
pl.data_ind = squeeze(mean(TFA.fftdata_ind(:,pl.elec2plot_i,:,pl.time2plot,pl.sub2plot),[2,3,4]));
pl.data_evo = squeeze(mean(TFA.fftdata_evo(:,pl.elec2plot_i,:,pl.time2plot,pl.sub2plot),[2,3,4]));

pl.foi = [TFA.RDK(1).RDK.freq];



% plotting
figure
set(gcf,'Position',[100 100 700 400],'PaperPositionMode','auto')
subplot(2,1,1)
plot(TFA.fftfreqs,squeeze(pl.data_ind),'Color',[0.5 0.5 0.5],'LineWidth',1)
hold on;
plot(TFA.fftfreqs,squeeze(mean(pl.data_ind,2)),'Color','k','LineWidth',2)

xlim([0 50])
xlabel('frequency in Hz')
ylabel('amplitude in \muV')
title(sprintf('induced GrandMean FFT spectra | N = %1.0f ',numel(pl.sub2plot)),'Interpreter','none')
vline(F.stim_frequency,'k:')

% draw topography with electrode positions
h.a1 = axes('position',[0.72 0.75 0.15 0.15],'Visible','off');
topoplot(find(pl.elec2plot_i),TFA(1).electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o',[255 133 4]./255,4,1}); % [204 107 36; 15 71 101; 131 208 173]./255


subplot(2,1,2)
plot(TFA.fftfreqs,squeeze(pl.data_evo),'Color',[0.5 0.5 0.5],'LineWidth',1)
hold on;
plot(TFA.fftfreqs,squeeze(mean(pl.data_evo,2)),'Color','k','LineWidth',2)

xlim([0 50])
xlabel('frequency in Hz')
ylabel('amplitude in \muV')
title(sprintf('evoked GrandMean FFT spectra | N = %1.0f ',numel(pl.sub2plot)),'Interpreter','none')
% vline(pl.foi,'k:')

% draw topography with electrode positions
h.a1 = axes('position',[0.72 0.28 0.15 0.15],'Visible','off');
topoplot(find(pl.elec2plot_i),TFA(1).electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o',[255 133 4]./255,4,1}); % [204 107 36; 15 71 101; 131 208 173]./255
box on


figure
set(gcf,'Position',[100 100 500 250],'PaperPositionMode','auto')
plot(TFA.fftfreqs,squeeze(pl.data_evo),'Color',[0.5 0.5 0.5],'LineWidth',1)
hold on;
plot(TFA.fftfreqs,squeeze(mean(pl.data_evo,2)),'Color','k','LineWidth',2)

xlim([0 50])
xlabel('frequency in Hz')
ylabel('amplitude in \muV')
title(sprintf('evoked GrandMean FFT spectra | N = %1.0f',numel(pl.sub2plot)),'Interpreter','none')
vline(pl.foi,'k:')

% draw topography with electrode positions
h.a1 = axes('position',[0.65 0.60 0.3 0.3],'Visible','off');
topoplot(find(pl.elec2plot_i),TFA(1).electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o',[255 133 4]./255,4,1}); % [204 107 36; 15 71 101; 131 208 173]./255

box on




% figure('Position',[100 100 350 300]);
sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_FShiftAlpha\figures\';
sav.filenames = 'FFT_Spectra_GrandMean_pre_cue_18Hz';
% print(gcf, fullfile(sav.pathout,sav.filenames),'-dpng','-r600')
% saveas(gcf,fullfile(sav.pathout,sav.filenames),'fig')
% print(gcf,fullfile(sav.pathout,sav.filenames),'-depsc2', '-painters','-r300')

%% grand mean alpha spectrum
% plotting parameters
pl.elec2plot = {{'P6';'P8';'P10';'PO4';'PO8';'O2';'I2'}; ... % for left stimulus/left hemifield SSVEP
    {'P9';'P7';'P5';'PO7';'PO3';'I1';'O1'}};
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), [pl.elec2plot{1}; pl.elec2plot{2}], 'UniformOutput',false)),1));


% plot channels with labels
% findchannellabel( 1:64 , 1)

pl.time2plot = [1]; % [1 2 4]
pl.sub2plot = [1:numel(F.Subs2use)]; % sub 15 ith no SSVEP default = 1:numel(F.Subs2use)


% extract data
pl.data_ind = squeeze(mean(TFA.fftdata_ind(:,pl.elec2plot_i,:,pl.time2plot,pl.sub2plot),[2 3 4]));



% plotting
figure
set(gcf,'Position',[100 100 700 300],'PaperPositionMode','auto')
plot(TFA.fftfreqs,squeeze(pl.data_ind),'Color',[0.5 0.5 0.5],'LineWidth',1)
hold on;
plot(TFA.fftfreqs,squeeze(mean(pl.data_ind,2)),'Color','k','LineWidth',2)

xlim([0 50])
xlabel('frequency in Hz')
ylabel('amplitude in \muV')
title(sprintf('induced GrandMean FFT spectra | alpha electrode cluster | N = %1.0f',numel(pl.sub2plot)),'Interpreter','none')
% vline(pl.foi,'k:')


% draw topography with electrode positions
h.a1 = axes('position',[0.69 0.65 0.25 0.25],'Visible','off');
topoplot(find(pl.elec2plot_i),TFA(1).electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o',[255 133 4]./255,4,1}); % [204 107 36; 15 71 101; 131 208 173]./255

sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_FShiftAlpha\figures\';
sav.filenames = 'FFT_Spectra_GrandMean_pre_cue_alpha';
% print(gcf, fullfile(sav.pathout,sav.filenames),'-dpng','-r600')
% saveas(gcf,fullfile(sav.pathout,sav.filenames),'fig')
% print(gcf,fullfile(sav.pathout,sav.filenames),'-depsc2', '-painters','-r300')

%% plot Grand Mean FFT data | topoplot for different frequencies of stimuli
pl.time2plot = [1];
pl.freqrange=[-0.1 0.1];
pl.sub2plot = [1:numel(F.Subs2use)]; % sub 15 ith no SSVEP default = 1:numel(F.Subs2use)

% extract data
pl.data_ind = []; pl.data_evo = []; pl.title = {}; h.s =[];

for i_freq = 1:numel(F.stim_frequency) % loop across positions
    % index frequencies of interest 
    t.fidx = dsearchn(TFA.fftfreqs', (F.stim_frequency(i_freq) + pl.freqrange)');

    pl.data_ind(i_freq,:) = squeeze(mean(TFA.fftdata_ind(t.fidx(1):t.fidx(2),:,:,pl.time2plot,pl.sub2plot),[1,3,5]));
    pl.data_evo(i_freq,:) = squeeze(mean(TFA.fftdata_evo(t.fidx(1):t.fidx(2),:,:,pl.time2plot,pl.sub2plot),[1,3,5]));
end



figure;
set(gcf,'Position',[100 100 700 300],'PaperPositionMode','auto')
for i_freq = 1:size(pl.data_ind,1)
    h.s(i_freq,1)=subplot(2,size(pl.data_ind,1),i_freq);
    
    topoplot( pl.data_ind(i_freq,:), TFA.electrodes(1:64), ...
        'shading', 'interp', 'numcontour', 0, 'maplimits',[0 max(pl.data_ind,[],'all')],'conv','on','colormap',fake_parula,...
        'whitebk','on');
    title(sprintf('induced %1.0fHz',F.stim_frequency(i_freq)))
    colorbar
    
    h.s(i_freq,2)=subplot(2,size(pl.data_ind,1),i_freq+size(pl.data_ind,1));
    
%     topoplot( t.pldata(:,i_freq), TFA.electrodes(1:64), ...
%         'shading', 'interp', 'numcontour', 0, 'maplimits',[0 max(t.pldata(:))],'conv','on','colormap',fake_parula,...
%         'whitebk','on');
%     topoplot( pl.data_evo(i_freq,:), TFA.electrodes(1:64), ...
%         'shading', 'interp', 'numcontour', 0, 'maplimits',[0 max(pl.data_evo,[],'all')],'conv','on','colormap',fake_parula,...
%         'whitebk','on');
    topoplot( pl.data_evo(i_freq,:), TFA.electrodes(1:64), ...
        'shading', 'interp', 'numcontour', 0, 'maplimits',[0 max( pl.data_evo(i_freq,:),[],'all')],'conv','on','colormap',fake_parula,...
        'whitebk','on');
    title(sprintf('evoked %1.0fHz',F.stim_frequency(i_freq)))
    colorbar
end

figure;
set(gcf,'Position',[100 100 700 200],'PaperPositionMode','auto')
for i_freq = 1:size(pl.data_evo,1)
    h.s(i_freq,1)=subplot(1,size(pl.data_evo,1),i_freq);
    t.pldata = squeeze(mean( pl.data_evo,3));
%     topoplot( t.pldata(:,i_freq), TFA.electrodes(1:64), ...
%         'shading', 'interp', 'numcontour', 0, 'maplimits',[0 max(t.pldata(:))],'conv','on','colormap',fake_parula,...
%         'whitebk','on');
    topoplot( pl.data_evo(i_freq,:), TFA.electrodes(1:64), ...
        'shading', 'interp', 'numcontour', 0, 'maplimits',[0 max( pl.data_evo(i_freq,:),[],'all')],'conv','on','colormap',fake_parula,...
        'whitebk','on');
    title(sprintf('evoked %1.0fHz',F.stim_frequency(i_freq)))
    colorbar
end

% % % draw colorbar into new axis
% t.pos2 = get(h.s{i_pl,1},'Position');
% t.pos3 = get(h.s{i_pl,1},'OuterPosition');
% h.a1 = axes('position',[t.pos3(1) t.pos2(2) t.pos3(3) t.pos2(4)],'Visible','off');
% caxis(pl.clims(1,:));
% colormap(gca,'jet')
% h.cb{i_pl,1}=colorbar;
% t.pos4 = get(h.cb{i_pl,1},'Position');
% set(h.cb{i_pl,1},'Position',[t.pos4(1)+0.04 t.pos4(2)+(t.pos4(4)*(1/6)) t.pos4(3)/2 t.pos4(4)*(2/3)])

sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_FShiftAlpha\figures\';
sav.filenames = 'FFT_TOPO_GrandMean_pre_cue_SSVEP';
% print(gcf, fullfile(sav.pathout,sav.filenames),'-dpng','-r600')
% saveas(gcf,fullfile(sav.pathout,sav.filenames),'fig')
% print(gcf,fullfile(sav.pathout,sav.filenames),'-depsc2', '-painters','-r300')




%% plot alpha grand mean topographies
pl.time2plot_pre = [1];
pl.time2plot_post = [3];
pl.sub2plot = [1:numel(F.Subs2use)]; % sub 15 ith no SSVEP default = 1:numel(F.Subs2use)
pl.alpha_range = [8 12];
% pl.alpha_range = [25 30];
pl.alpha_range_i = dsearchn(TFA.fftfreqs', pl.alpha_range');
pl.rdk_att_index = [1 2 1];


% extract data
pl.data_ind = []; %[elec x time x sub]
pl.title = {}; h.s =[];

pl.data_ind(:,1,:) = ...
    squeeze(mean(TFA.fftdata_ind(...
    pl.alpha_range_i(1):pl.alpha_range_i(2),:,:,pl.time2plot_pre,pl.sub2plot),...
    [1,3,4]));
pl.data_ind(:,2,:) = ...
    squeeze(mean(TFA.fftdata_ind(...
    pl.alpha_range_i(1):pl.alpha_range_i(2),:,:,pl.time2plot_post,pl.sub2plot),...
    [1,3,4]));

% plot raw data
pl.data = pl.data_ind;
pl.clim1 = [0 max(mean(pl.data,3),[],'all')];
pl.lab_time = {'pre';'post'};

figure;
set(gcf,'Position',[100 100 1000 200],'PaperPositionMode','auto')
h.s = [];
% plot data

for i_time = 1:2
    h.s(i_time)=subplot(1,4,i_time);
    t.pldata = mean(pl.data(:,i_time,:),3);
    
    topoplot( t.pldata, TFA.electrodes(1:64), ...
        'numcontour', 0, 'maplimits',pl.clim1,'conv','on','colormap',fake_parula,...
        'whitebk','on');
    title(pl.lab_time{i_time})
    colorbar
    
end

% difference
h.s(i_time)=subplot(1,4,i_time+1);
t.pldata = mean(pl.data(:,2,:)-pl.data(:,1,:),3);

topoplot( t.pldata, TFA.electrodes(1:64), ...
    'numcontour', 0, 'maplimits','absmax','conv','on','colormap',flipud(cbrewer2('RdBu')),...
    'whitebk','on');
title('diff')
colorbar

% modulation in percent
h.s(i_time)=subplot(1,4,i_time+2);
t.pldata = mean(((pl.data(:,2,:)./pl.data(:,1,:)-1).*100),3);
topoplot( t.pldata, TFA.electrodes(1:64), ...
    'numcontour', 0, 'maplimits','absmax','conv','on','colormap',flipud(cbrewer2('RdBu')),...
    'whitebk','on');
title('pre-to-post modulation')
colorbar





sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_FShiftAlpha\figures\';
sav.filenames = 'FFT_TOPO_GrandMean_alpha';
% print(gcf, fullfile(sav.pathout,sav.filenames),'-dpng','-r600')
% saveas(gcf,fullfile(sav.pathout,sav.filenames),'fig')
% print(gcf,fullfile(sav.pathout,sav.filenames),'-depsc2', '-painters','-r300')



%% plot only electrode positions
% % large cluster
% pl.elec2plot_cluster = {{'O1';'Oz';'O2';'I1';'Iz';'I2'};...
%     {'P1';'Pz';'P2';'PO3';'POz';'PO4';'O1';'Oz';'O2';'P6';'P8';'PO8';'P10';'I2'};...
%     {'P1';'Pz';'P2';'PO3';'POz';'PO4';'O1';'Oz';'O2';'P5';'P7';'PO7';'P9';'I1'}};
% 
% % central
% pl.elec2plot_cluster = {{'Oz';'Iz'};...
%     {'POz';'PO3'};...
%     {'POz';'PO4'}};
% % 
% % peripheral
% pl.elec2plot_cluster = {{'Oz';'Iz'};...
%     {'PO8';'P8'};...
%     {'PO7';'P7'}};
% 
% pl.elec2plot_cluster_label = {'central';'left';'right'}; % which stimulus is best captured?
% pl.col = {[255 133 4]./255; [15 71 101]./255; [15 71 101]./255};
% figure;
% set(gcf,'Position',[100 100 100 500],'PaperPositionMode','auto')


pl.elec2plot_cluster = {{'P6';'P8';'P10';'PO4';'PO8';'O2';'I2'}; ... % for left stimulus/left hemifield SSVEP
    {'P9';'P7';'P5';'PO7';'PO3';'I1';'O1'}}; % for right stimulus/right hemifield SSVEP
pl.elec2plot_cluster_label = {'left';'right'}; % which stimulus is best captured?
pl.col = { [15 71 101]./255; [15 71 101]./255};
figure;
set(gcf,'Position',[100 100 100 300],'PaperPositionMode','auto')


for i_clust = 1:numel(pl.elec2plot_cluster_label)
    subplot(numel(pl.elec2plot_cluster_label),1,i_clust)
    pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x)strcmpi({TFA.electrodes.labels},x),pl.elec2plot_cluster{i_clust}, 'UniformOutput',false)),1));
    topoplot(find(pl.elec2plot_i),TFA(1).electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
        'emarker2',{find(pl.elec2plot_i),'o',pl.col{i_clust},4,1}); % [204 107 36; 15 71 101; 131 208 173]./255
end
sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_FShiftAlpha\figures\';
sav.filenames = 'ELECS_POS_alpha';
% print(gcf, fullfile(sav.pathout,sav.filenames),'-dpng','-r600')
% saveas(gcf,fullfile(sav.pathout,sav.filenames),'fig')
% print(gcf,fullfile(sav.pathout,sav.filenames),'-depsc2', '-painters','-r300')

%% extract amplitude values for FFT-derived SSVEPs
pl.time2plot = [1 2 3 4];
pl.freqrange=[-0.1 0.1];
pl.sub2plot = 1:numel(F.Subs2use);


% standard cluster:
pl.elec2plot_cluster = {{'O1';'Oz';'O2';'I1';'Iz';'I2'};...
    {'O1';'Oz';'O2';'I1';'Iz';'I2'};...
    {'O1';'Oz';'O2';'I1';'Iz';'I2'}};

% small cluster
pl.elec2plot_cluster = {{'POz';'O1';'Oz';'O2';'Iz'};...
    {'POz';'O1';'Oz';'O2';'Iz'};...
    {'POz';'O1';'Oz';'O2';'Iz'}};

% % broad cluster
% pl.elec2plot_cluster = { {'P5';'P7';'PO7';'PO3'; 'O1';'Oz';'O2';'I1';'Iz';'I2'; 'PO4'; 'P8';'P6';'PO8'};...
%     {'P5';'P7';'PO7';'PO3'; 'O1';'Oz';'O2';'I1';'Iz';'I2'; 'PO4'; 'P8';'P6';'PO8'};...
%     {'P5';'P7';'PO7';'PO3'; 'O1';'Oz';'O2';'I1';'Iz';'I2'; 'PO4'; 'P8';'P6';'PO8'}};


pl.elec2plot_cluster_label = {'central';'central';'central'}; % which stimulus is best captured?
pl.data_ind = []; pl.data_evo = [];
pl.con1label = {'valid';'valid';'invalid';'invalid';'neutral';'neutral'}; % cue validity
pl.con2label = {'RDK1';'RDK2';'RDK2';'RDK1';'RDK1';'RDK2'}; % which RDK is target?
pl.con3label = {'RDK1';'RDK2';'RDK1';'RDK2';'both';'both'}; % which RDK is attended
pl.RDKlabel = {'RDK1';'RDK2';'RDK3'};
pl.timelabel = cellfun(@(x) vararg2str(x),TFA.ffttimewin(pl.time2plot),'UniformOutput',false);

% extract data
pl.data_ind = nan(numel(pl.RDKlabel), size(TFA.fftdata_ind,3), numel(pl.time2plot), numel(pl.sub2plot));
pl.data_evo = pl.data_ind;
pl.label.color_label = repmat({''},numel(pl.RDKlabel),numel(pl.sub2plot));
pl.label.freq = repmat({[]},numel(pl.RDKlabel),numel(pl.sub2plot));
pl.label.elec = repmat({[]},numel(pl.RDKlabel),numel(pl.sub2plot));


% extract data
for i_sub = 1:numel(pl.sub2plot)
    for i_RDK = 1:numel(pl.RDKlabel) % loop across positions
        % extract all relevant information for RDK
        pl.label.color_label{i_RDK,i_sub} = TFA.RDK(pl.sub2plot(i_sub)).RDK(i_RDK).col_label;
        pl.label.freq{i_RDK,i_sub} = TFA.RDK(pl.sub2plot(i_sub)).RDK(i_RDK).freq;
        
        % index color for data extraction
        t.fidx1=cellfun(@(x) x+pl.freqrange, pl.label.freq(i_RDK,i_sub), 'UniformOutput', false);
        t.fidx2 = cell2mat(cellfun(@(x) dsearchn(TFA.fftfreqs',(x(1))):dsearchn(TFA.fftfreqs',(x(2))),t.fidx1,'UniformOutput',false));
        
        pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), pl.elec2plot_cluster{i_RDK}, 'UniformOutput',false)),1));
        pl.label.elec{i_RDK,i_sub} = vararg2str(pl.elec2plot_cluster{i_RDK});
        
        
        
        
        pl.data_ind(i_RDK,:,:,i_sub) = squeeze(mean(TFA.fftdata_ind(t.fidx2,pl.elec2plot_i,:,pl.time2plot,pl.sub2plot(i_sub)),[1,2]));
        pl.data_evo(i_RDK,:,:,i_sub) = squeeze(mean(TFA.fftdata_evo(t.fidx2,pl.elec2plot_i,:,pl.time2plot,pl.sub2plot(i_sub)),[1,2]));
    end
end



% into long format
pl.RDK.data_ind_l = pl.data_ind(:);
pl.RDK.data_evo_l = pl.data_evo(:);
% pl.data_ind_bc = 100*(bsxfun(@rdivide, pl.data_ind, mean(pl.data_ind(:,:,1,:),2))-1);
% pl.data_evo_bc = 100*(bsxfun(@rdivide, pl.data_evo, mean(pl.data_evo(:,:,1,:),2))-1);
pl.data_ind_bc = 100*(bsxfun(@rdivide, pl.data_ind, pl.data_ind(:,:,1,:))-1);
pl.data_evo_bc = 100*(bsxfun(@rdivide, pl.data_evo, pl.data_evo(:,:,1,:))-1);
% pl.data_ind_bc = bsxfun(@minus, pl.data_ind, pl.data_ind(:,:,1,:));
% pl.data_evo_bc = bsxfun(@minus, pl.data_evo, pl.data_evo(:,:,1,:));
pl.RDK.data_ind_bc_l = pl.data_ind_bc(:);
pl.RDK.data_evo_bc_l = pl.data_evo_bc(:);
% pl.RDK.subs = repmat(pl.sub2plot,prod(size(pl.data_ind,[1 2 3])),1); pl.RDK.subs = pl.RDK.subs(:);
pl.RDK.subs = repmat(F.Subs2use(pl.sub2plot),prod(size(pl.data_ind,[1 2 3])),1); pl.RDK.subs = pl.RDK.subs(:);
pl.RDK.RDK = repmat(pl.RDKlabel,1,prod(size(pl.data_ind,[2 3 4]))); pl.RDK.RDK=pl.RDK.RDK(:);
pl.RDK.col = permute(repmat(pl.label.color_label,[1, 1, prod(size(pl.data_ind,[2 3]))]), [1 3 2]); pl.RDK.col=pl.RDK.col(:);
pl.RDK.elec = permute(repmat(pl.label.elec,[1, 1, prod(size(pl.data_ind,[2 3]))]), [1 3 2]); pl.RDK.elec=pl.RDK.elec(:);
pl.RDK.freq = permute(repmat(pl.label.freq,[1, 1, prod(size(pl.data_ind,[2 3]))]), [1 3 2]); pl.RDK.freq=pl.RDK.freq(:);
pl.RDK.con1 = repmat(pl.con1label',size(pl.data_ind,1),prod(size(pl.data_ind,[3 4]))); pl.RDK.con1 = pl.RDK.con1(:);
pl.RDK.con2 = repmat(pl.con2label',size(pl.data_ind,1),prod(size(pl.data_ind,[3 4]))); pl.RDK.con2 = pl.RDK.con2(:);
pl.RDK.con3 = repmat(pl.con3label',size(pl.data_ind,1),prod(size(pl.data_ind,[3 4]))); pl.RDK.con3 = pl.RDK.con3(:);
pl.RDK.time = repmat(pl.timelabel',prod(size(pl.data_ind,[1 2])),prod(size(pl.data_ind,[4]))); pl.RDK.time = pl.RDK.time(:);

% remove grand mean induced data
pl.data_evo_bc = bsxfun(@minus, pl.data_evo_bc, mean(pl.data_ind_bc, [2,4]) );
pl.RDK.data_evo_bc_l = pl.data_evo_bc(:);

% plot 1
clear g
% plot with gramm
g(1,1)=gramm('x',pl.RDK.time,'y',pl.RDK.data_ind_l,'color',pl.RDK.con1);
g(1,1).facet_grid([],pl.RDK.RDK);
g(1,1).stat_summary('geom',{'bar','black_errorbar'});
g(1,1).set_names('x','time windows','y','amplitude in \muV','color','attention','column','SSVEP');
g(1,1).set_title('amplitude values of SSVEPs (induced analysis)');
% g.stat_boxplot();
% figure('Position',[100 100 800 350]);
% g.draw();

% % troubleshooting
% [pl.RDK.time pl.RDK.con pl.RDK.RDK num2cell(pl.RDK.data_ind_l)]

g(2,1)=gramm('x',pl.RDK.time,'y',pl.RDK.data_evo_l,'color',pl.RDK.con1);
g(2,1).facet_grid([],pl.RDK.RDK);
g(2,1).stat_summary('geom',{'bar','black_errorbar'});
g(2,1).set_names('x','time windows','y','amplitude in \muV','color','attention','column','SSVEP');
g(2,1).set_title('amplitude values of SSVEPs (evoked analysis)');
% g.stat_boxplot();
% figure('Position',[100 100 800 350]);
% g.draw();

% plot baseline corrected data (modulation in %)
% plot with gramm
g(1,2)=gramm('x',pl.RDK.time,'y',pl.RDK.data_ind_bc_l,'color',pl.RDK.con1);
g(1,2).facet_grid([],pl.RDK.RDK);
g(1,2).stat_summary('geom',{'bar','black_errorbar'});
g(1,2).set_names('x','time windows','y','modulation in %','color','attention','column','SSVEP');
g(1,2).set_title('amplitude modulations of SSVEPs (induced analysis)');
% g.stat_boxplot();
% figure('Position',[100 100 800 350]);
% g.draw();

g(2,2)=gramm('x',pl.RDK.time,'y',pl.RDK.data_evo_bc_l,'color',pl.RDK.con1);
g(2,2).facet_grid([],pl.RDK.RDK);
g(2,2).stat_summary('geom',{'bar','black_errorbar'});
g(2,2).set_names('x','time windows','y','modulation in %','color','attention','column','SSVEP');
g(2,2).set_title('amplitude modulations of SSVEPs (evoked analysis)');
% g.stat_boxplot();
g.set_text_options('base_size',8);
g.axe_property('YGrid','on','Box','on');
g.set_color_options('map',[255 133 4; 41 60 74; 25 138 131]./255);
figure('Position',[100 100 1600 700]);
g.draw();


% long data not collapsed
R_Mat.all = cell2struct([num2cell([pl.RDK.data_ind_l pl.RDK.data_evo_l pl.RDK.data_ind_bc_l pl.RDK.data_evo_bc_l pl.RDK.subs]) ...
    pl.RDK.RDK pl.RDK.col  pl.RDK.freq pl.RDK.elec pl.RDK.con1 pl.RDK.con2 pl.RDK.con3 pl.RDK.time]',...
    {'amplitude_induced','amplitude_evoked','modulation_induced','modulation_evoked','subjects',...
    'RDK_id', 'RDK_col', 'RDK_freq', 'RDK_elec', 'cue_validity', 'target_RDK', 'attended_RDK' 'time'});
R_Mat.all_t = struct2table(R_Mat.all);
t.path = 'C:\Users\psy05cvd\Dropbox\work\R-statistics\experiments\ssvep_fshiftprobabil\data';
t.datestr = datestr(now,'mm-dd-yyyy_HH-MM');
% t.filename = 'FFT_SSVEP_Amp_data_withoutBehav_sepRDK_peripheral_2elecs';
t.filename = 'FFT_SSVEP_largeclust';
% t.filename = 'FFT_SSVEP_smallcentralclust';
t.filename = 'FFT_SSVEP_smallcentralclust_fiveelecs';
% write to textfile
% writetable(R_Mat.all_t,fullfile(t.path,sprintf('%s.csv',t.filename)),'Delimiter',';')






%% plot SSVEP timecourse
pl = [];
pl.time2plot = [1 2 3 4];
pl.sub2plot = [1:numel(F.Subs2use)]; % sub 15 ith no SSVEP default = 1:numel(F.Subs2use)
% pl.sub2plot = [1:11 13:14 16:20 22:35]; % sub 15 with no SSVEP default = 1:numel(F.Subs2use); 12 21 bad in baseline?
pl.con1label = {'valid';'valid';'invalid';'invalid';'neutral';'neutral'}; % cue validity
pl.con2label = {'RDK1';'RDK2';'RDK2';'RDK1';'RDK1';'RDK2'}; % which RDK is target?
pl.con3label = {'RDK1';'RDK2';'RDK1';'RDK2';'both';'both'}; % which RDK is attended
pl.RDKlabel = {'RDK1';'RDK2';'RDK3'};
pl.base = F.TFA.baseline;
pl.base = [-250 0];
pl.base = [-500 -250];
pl.base_i = dsearchn( TFA.time', pl.base');
pl.xlims=[-500 1750]; % index time 2 plot
pl.xlims_i = dsearchn(TFA.time', pl.xlims');

pl.conlabel = {'cued';'uncued';'neutral';'irr_cue';'irr_neutr'};
% pl.concols = num2cell([255 133 4; 25 138 131; 255 194 64; 37 207 197]'./255,1);
pl.concols = num2cell([25 138 130; 241 131 26; 41 60 74; 240 63 63; 255 120 120]'./255,1);


% % standard cluster:
% pl.elec2plot_cluster = {{'O1';'Oz';'O2';'I1';'Iz';'I2'};...
%     {'O1';'Oz';'O2';'I1';'Iz';'I2'};...
%     {'O1';'Oz';'O2';'I1';'Iz';'I2'}};

pl.elec2plot_cluster = {{'O1';'Oz';'O2';'Iz'};...
    {'O1';'Oz';'O2';'Iz'};...
    {'O1';'Oz';'O2';'Iz'}};

pl.elec2plot_cluster = {{'O1';'Oz';'O2';'I1';'Iz';'I2'};...
    {'O1';'Oz';'O2';'Iz'};...
    {'O1';'Oz';'O2';'Iz'}};

% pl.elec2plot_cluster = { {'P5';'P7';'PO7';'PO3'; 'O1';'Oz';'O2';'I1';'Iz';'I2'; 'PO4'; 'P8';'P6';'PO8'};...
%     {'P5';'P7';'PO7';'PO3'; 'O1';'Oz';'O2';'I1';'Iz';'I2'; 'PO4'; 'P8';'P6';'PO8'};...
%     {'P5';'P7';'PO7';'PO3'; 'O1';'Oz';'O2';'I1';'Iz';'I2'; 'PO4'; 'P8';'P6';'PO8'}};

% exploratory cluster:
% % central
% pl.elec2plot_cluster = {{'Oz';'Iz'};...
%     {'POz';'PO3'};...
%     {'POz';'PO4'}};
% 
% % peripheral
% pl.elec2plot_cluster = {{'Oz';'Iz'};...
%     {'PO8';'P8'};...
%     {'PO7';'P7'}};


pl.elec2plot_cluster_label = {'central';'central';'central'}; % which stimulus is best captured?
% extract data
pl.data_RDK = []; pl.data_RDK_ind = [];pl.label.color_label = [];l.label.freq=[];
for i_sub = 1:numel(pl.sub2plot)
    for i_RDK = 1:numel(pl.RDKlabel) % loop across positions
        % extract all relevant information for RDK
        pl.label.color_label{i_RDK,i_sub} = TFA.RDK(pl.sub2plot(i_sub)).RDK(i_RDK).col_label;
        pl.label.freq{i_RDK,i_sub} = TFA.RDK(pl.sub2plot(i_sub)).RDK(i_RDK).freq;
        
        % index freq for data extraction
        t.fidx1 = dsearchn(TFA.frequency', pl.label.freq{i_RDK,i_sub});
        
        pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), pl.elec2plot_cluster{i_RDK}, 'UniformOutput',false)),1));
        pl.label.elec{i_RDK,i_sub} = vararg2str(pl.elec2plot_cluster{i_RDK});
        
        % extract data
        pl.data_RDK(i_RDK,:,:,i_sub) = ...
            squeeze(mean(TFA.data_evo(t.fidx1,:,pl.elec2plot_i,:,pl.sub2plot(i_sub)),3));
        pl.data_RDK_ind(i_RDK,:,:,i_sub) = ...
            squeeze(mean(TFA.data_ind(t.fidx1,:,pl.elec2plot_i,:,pl.sub2plot(i_sub)),3));
            
    end
end

% collapse across RDKs and conditions to represent data like [cued, uncued, neutral] x [task relevant, task irrelevant]
%#######
% induced
t.ridx = {'RDK1'; 'RDK2';'RDK1'};
pl.tdata = [];
for i_rdk = 1:numel(t.ridx)-1
    t.idx = strcmp(pl.con3label,t.ridx{i_rdk}); % index attended RDK
    pl.tdata(i_rdk,:,1,:) = mean(pl.data_RDK_ind( i_rdk,:,t.idx,:),3);
    t.idx = strcmp(pl.con3label,t.ridx{i_rdk+1});  % index unattended RDK
    pl.tdata(i_rdk,:,2,:) = mean(pl.data_RDK_ind( i_rdk,:,t.idx,:),3);
    t.idx = strcmp(pl.con3label,'both'); % index neutral trials
    pl.tdata(i_rdk,:,3,:) = mean(pl.data_RDK_ind( i_rdk,:,t.idx,:),3);
end
% collapse across RDKs
pl.tdata = mean(pl.tdata,1);

% get third RDK and collapse for single attention cue or both color cue
t.idx = ~strcmp(pl.con1label,'neutral');
pl.tdata(:,:,4,:)=mean(pl.data_RDK_ind(3,:,t.idx,:),3); % index single attended color
pl.tdata(:,:,5,:)=mean(pl.data_RDK_ind(3,:,~t.idx,:),3); % index neutral
pl.data_ind = squeeze(pl.tdata);

%#######
% evoked
t.ridx = {'RDK1'; 'RDK2';'RDK1'};
pl.tdata = [];
for i_rdk = 1:numel(t.ridx)-1
    t.idx = strcmp(pl.con3label,t.ridx{i_rdk}); % index attended RDK
    pl.tdata(i_rdk,:,1,:) = mean(pl.data_RDK( i_rdk,:,t.idx,:),3);
    t.idx = strcmp(pl.con3label,t.ridx{i_rdk+1});  % index unattended RDK
    pl.tdata(i_rdk,:,2,:) = mean(pl.data_RDK( i_rdk,:,t.idx,:),3);
    t.idx = strcmp(pl.con3label,'both'); % index neutral trials
    pl.tdata(i_rdk,:,3,:) = mean(pl.data_RDK( i_rdk,:,t.idx,:),3);
end
% collapse across RDKs
pl.tdata = mean(pl.tdata,1);

% get third RDK and collapse for single attention cue or both color cue
t.idx = ~strcmp(pl.con1label,'neutral');
pl.tdata(:,:,4,:)=mean(pl.data_RDK(3,:,t.idx,:),3);
pl.tdata(:,:,5,:)=mean(pl.data_RDK(3,:,~t.idx,:),3);
pl.data = squeeze(pl.tdata);


% baselince corrected
pl.data_bc = 100.*(bsxfun(@rdivide, pl.data, mean(pl.data(pl.base_i(1):pl.base_i(2),:,:),1))-1);
pl.data_ind_bc = 100.*(bsxfun(@rdivide, pl.data_ind, mean(pl.data_ind(pl.base_i(1):pl.base_i(2),:,:),1))-1);

% mean and sem data
pl.mdata = mean(pl.data,3);
pl.semdata = std(pl.data,1,3)./sqrt(numel(pl.sub2plot));
pl.mdata_bc = mean(pl.data_bc,3);
pl.semdata_bc = std(pl.data_bc,1,3)./sqrt(numel(pl.sub2plot));
pl.mdata_ind = mean(pl.data_ind,3);
pl.semdata_ind = std(pl.data_ind,1,3)./sqrt(numel(pl.sub2plot));
pl.mdata_ind_bc = mean(pl.data_ind_bc,3);
pl.semdata_ind_bc = std(pl.data_ind_bc,1,3)./sqrt(numel(pl.sub2plot));



% plot data with line plots
%  evoked raw
pl.col = pl.concols;
pl.col2 = [0.6 0.6 0.6];
pl.line = {'-';'-';'-';'-';'-'};
figure;
set(gcf,'Position',[100 100 700 400],'PaperPositionMode','auto')
h.pl = {}; h.plsem=[];  h.plm = []; h.pls = []; h.plst = []; pl.con_label = [];
t.idx = 1;

for i_con = 1:size(pl.mdata,2)
    % data index
    pl.idx = pl.xlims_i(1):pl.xlims_i(2);

    % plot SEM as boundary
    % create data
    pl.xconf = [TFA.time(pl.idx) TFA.time(pl.idx(end:-1:1))] ;
    pl.yconf = [pl.mdata(pl.idx,i_con)'+pl.semdata(pl.idx,i_con)' ...
        pl.mdata(pl.idx(end:-1:1),i_con)'-pl.semdata(pl.idx(end:-1:1),i_con)'];
    % plot
    h.plsem{t.idx} = fill(pl.xconf,pl.yconf,pl.col{t.idx}','EdgeColor','none','FaceAlpha',0.3);
    hold on

    % plot mean lines
    h.plm{t.idx}=plot(TFA.time(pl.idx), pl.mdata(pl.idx,i_con),'Color',pl.col{t.idx},'LineStyle',pl.line{t.idx},'LineWidth',2);

    

    t.idx = t.idx+1;
end

title(sprintf('raw SSVEP amplitude timecourses | %s',vararg2str(pl.elec2plot_cluster{1})))
xlim(pl.xlims)
xlabel('time in ms')
ylabel('amplitude in \muV/cm²')
legend([h.plm{:}],pl.conlabel,'Location','SouthOutside','Orientation','horizontal','Interpreter','none')
legend('boxoff')
grid on

sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_FShiftAlpha\figures\';
sav.filenames = 'TFA_SSVEPmod_largeclust_inducedremoved';
% print(gcf, fullfile(sav.pathout,sav.filenames),'-dpng','-r600')
% saveas(gcf,fullfile(sav.pathout,sav.filenames),'fig')
% print(gcf,fullfile(sav.pathout,sav.filenames),'-depsc2', '-painters','-r300')


%  evoked bc
pl.col = pl.concols;
pl.col2 = [0.6 0.6 0.6];
pl.line = {'-';'-';'-';'-';'-'};
figure;
set(gcf,'Position',[100 100 700 400],'PaperPositionMode','auto')
h.pl = {}; h.plsem=[];  h.plm = []; h.pls = []; h.plst = []; pl.con_label = [];
t.idx = 1;

for i_con = 1:size(pl.mdata,2)
    % data index
    pl.idx = pl.xlims_i(1):pl.xlims_i(2);

    % plot SEM as boundary
    % create data
    pl.xconf = [TFA.time(pl.idx) TFA.time(pl.idx(end:-1:1))] ;
    pl.yconf = [pl.mdata_bc(pl.idx,i_con)'+pl.semdata_bc(pl.idx,i_con)' ...
        pl.mdata_bc(pl.idx(end:-1:1),i_con)'-pl.semdata_bc(pl.idx(end:-1:1),i_con)'];
    % plot
    h.plsem{t.idx} = fill(pl.xconf,pl.yconf,pl.col{t.idx}','EdgeColor','none','FaceAlpha',0.3);
    hold on

    % plot mean lines
    h.plm{t.idx}=plot(TFA.time(pl.idx), pl.mdata_bc(pl.idx,i_con),'Color',pl.col{t.idx},'LineStyle',pl.line{t.idx},'LineWidth',2);

    

    t.idx = t.idx+1;
end

title(sprintf('baseline-corrected SSVEP amplitude timecourses | %s | base: [%1.0f %1.0f]ms', ...
    vararg2str(pl.elec2plot_cluster{1}), pl.base))
xlim(pl.xlims)
xlabel('time in ms')
ylabel('modulation in %')
legend([h.plm{:}],pl.conlabel,'Location','SouthOutside','Orientation','horizontal','Interpreter','none')
legend('boxoff')
grid on

sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_FShiftAlpha\figures\';
sav.filenames = 'TFA_SSVEPmod_largeclust_inducedremoved';
% print(gcf, fullfile(sav.pathout,sav.filenames),'-dpng','-r600')
% saveas(gcf,fullfile(sav.pathout,sav.filenames),'fig')
% print(gcf,fullfile(sav.pathout,sav.filenames),'-depsc2', '-painters','-r300')






%% plot alpha timecourse
pl.sub2plot = 1:numel(F.Subs2use);
pl.sub2plot = [1:14 16:35]; % sub 15 ith no SSVEP default = 1:numel(F.Subs2use)
pl.alpha_range = [8 12];
pl.alpha_range_i = dsearchn(TFA.frequency', pl.alpha_range');

pl.lab_con = {'contra attended';'contra unattended'};
pl.cols_con = {[240 63 63];[25 58 61]};

pl.elec2plot_cluster = {{'P6';'P8';'P10';'PO4';'PO8';'O2';'I2'}; ... % for left stimulus/left hemifield SSVEP
    {'P9';'P7';'P5';'PO7';'PO3';'I1';'O1'}}; % for right stimulus/right hemifield SSVEP
pl.elec2plot_cluster_label = {'alpha_contra_left_stim';'alpha_contra_right_stim'}; % which stimulus is best captured?
pl.elec2plot_cluster_label2 = {'left';'right'}; % which stimulus is best captured?
pl.elec2plot_i=cellfun( @(y) ...
    logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), y, 'UniformOutput',false)),1)),...
    pl.elec2plot_cluster, 'UniformOutput',false);

pl.time2plot = [-750 2000];
pl.time2plot_i = dsearchn(TFA.time', pl.time2plot');

pl.freq2plot = [2 45];
pl.freq2plot_i = dsearchn(TFA.frequency', pl.freq2plot');

pl.base = [-500 -250];
pl.base_i = dsearchn(TFA.time', pl.base');


% extract data
pl.data_ind_prime = nan(numel(TFA.frequency), numel(TFA.time),2,numel(pl.sub2plot),2);

% extract data
for i_sub = 1:numel(pl.sub2plot)
    for i_pos = 1:2 % loop across positions
        t.pos_att = {TFA.RDK(pl.sub2plot(i_sub)).RDK(3:4).position};
        
        if strcmp(t.pos_att{i_pos},'left')
            t.att_idx = [1 2];
        else
            t.att_idx = [2 1];
        end
        
        % extract data
        pl.data_ind_prime(:,:,:,i_sub,i_pos) = ...
            squeeze(mean(TFA.data_ind(:,:,pl.elec2plot_i{i_pos},t.att_idx,i_sub),3));
        
        pl.data_ind_prime(:,:,:,i_sub,i_pos) = ...
            squeeze(mean(TFA.data_evo(:,:,pl.elec2plot_i{i_pos},t.att_idx,i_sub),3));

    end
end

% average across left and right positions

pl.data_ind = mean(pl.data_ind_prime,5);
pl.data_ind_bc = 100.*(bsxfun(@rdivide, ...
    pl.data_ind, ...
    mean(pl.data_ind(:,pl.base_i(1):pl.base_i(2),:,:),2))-1);
pl.data_ind = pl.data_ind(pl.freq2plot_i(1):pl.freq2plot_i(2),pl.time2plot_i(1):pl.time2plot_i(2),:,:);
pl.data_ind_bc = pl.data_ind_bc(pl.freq2plot_i(1):pl.freq2plot_i(2),pl.time2plot_i(1):pl.time2plot_i(2),:,:);

% plot TFA plots
figure;
set(gcf,'Position',[100 100 1400 600],'PaperPositionMode','auto')
pl.clims1 = [0 max(mean(pl.data_ind,4),[],'all')];
pl.clims2 = [-1 1].*max(abs(mean(pl.data_ind_bc,4)),[],'all');



% raw data
for i_con = 1:2
    subplot(2,3,i_con)
    pl.data = mean(pl.data_ind(:,:,i_con,:),4);
    imagesc(TFA.time(pl.time2plot_i(1):pl.time2plot_i(2)),TFA.frequency(pl.freq2plot_i(1):pl.freq2plot_i(2)),pl.data,pl.clims1)
    colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    set(gca,'YDir','normal')
    title(sprintf('%s | raw | induced | N=%1.0f',pl.lab_con{i_con}, numel(pl.sub2plot)), ...
        'FontSize',8,'Interpreter','none')
    xlabel('time in ms')
    ylabel('frequency in Hz')
    set(gca,'FontSize',8)
    % set(gca,'ColorScale','log')
    cb = colorbar();
    
    if i_con == 1
        % draw topography with electrode positions
        h.a1 = axes('position',[0 0.45 0.14 0.14],'Visible','off');
        topoplot(1:64,TFA.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
            'emarker2',{find(sum([pl.elec2plot_i{1};pl.elec2plot_i{2}])),'o',[255 133 4]./255,4,1});
    end   
end

subplot(2,3,i_con+1)
pl.data = mean(pl.data_ind(:,:,1,:)-pl.data_ind(:,:,2,:),4);
pl.climst = [-1 1].*max(abs(pl.data),[],'all');
imagesc(TFA.time(pl.time2plot_i(1):pl.time2plot_i(2)),TFA.frequency(pl.freq2plot_i(1):pl.freq2plot_i(2)),pl.data,pl.climst)
colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
set(gca,'YDir','normal')
title(sprintf('diff | raw | induced | N=%1.0f', numel(pl.sub2plot)), ...
    'FontSize',8,'Interpreter','none')
xlabel('time in ms')
ylabel('frequency in Hz')
set(gca,'FontSize',8)
% set(gca,'ColorScale','log')
cb = colorbar();


% bc data
for i_con = 1:2
    subplot(2,3,i_con+3)
    pl.data = mean(pl.data_ind_bc(:,:,i_con,:),4);
    imagesc(TFA.time(pl.time2plot_i(1):pl.time2plot_i(2)),TFA.frequency(pl.freq2plot_i(1):pl.freq2plot_i(2)),pl.data,pl.clims2)
    colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    set(gca,'YDir','normal')
    title(sprintf('%s | modulation | induced | N=%1.0f',pl.lab_con{i_con}, numel(pl.sub2plot)), ...
        'FontSize',8,'Interpreter','none')
    xlabel('time in ms')
    ylabel('frequency in Hz')
    set(gca,'FontSize',8)
    % set(gca,'ColorScale','log')
    cb = colorbar();
    
    if i_con == 1
        % draw topography with electrode positions
        h.a1 = axes('position',[0 0.45 0.14 0.14],'Visible','off');
        topoplot(1:64,TFA.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
            'emarker2',{find(sum([pl.elec2plot_i{1};pl.elec2plot_i{2}])),'o',[255 133 4]./255,4,1});
    end   
end

subplot(2,3,i_con+1+3)
pl.data = mean(pl.data_ind_bc(:,:,1,:)-pl.data_ind_bc(:,:,2,:),4);
pl.climst = [-1 1].*max(abs(pl.data),[],'all');
imagesc(TFA.time(pl.time2plot_i(1):pl.time2plot_i(2)),TFA.frequency(pl.freq2plot_i(1):pl.freq2plot_i(2)),pl.data,pl.climst)
colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
set(gca,'YDir','normal')
title(sprintf('diff | modulation | induced | N=%1.0f', numel(pl.sub2plot)), ...
    'FontSize',8,'Interpreter','none')
xlabel('time in ms')
ylabel('frequency in Hz')
set(gca,'FontSize',8)
% set(gca,'ColorScale','log')
cb = colorbar();


% line plots
pl.col = pl.cols_con;
pl.col2 = [0.6 0.6 0.6];
pl.line = {'-';'-'};

pl.data_alpha = mean(pl.data_ind_prime(pl.alpha_range_i(1):pl.alpha_range_i(2),pl.time2plot_i(1):pl.time2plot_i(2),:,:,:),[1,5]);
pl.data_alpha_bc = 100.*(bsxfun(@rdivide, ...
    mean(pl.data_ind_prime(pl.alpha_range_i(1):pl.alpha_range_i(2),pl.time2plot_i(1):pl.time2plot_i(2),:,:,:),[1,5]),...
    mean(pl.data_ind_prime(pl.alpha_range_i(1):pl.alpha_range_i(2),pl.base_i(1):pl.base_i(2),:,:,:),[1,2,5]))-1);

pl.mdata_alpha = squeeze(mean(pl.data_alpha,4));
pl.mdata_alpha_bc = squeeze(mean(pl.data_alpha_bc,4));
pl.semdata_alpha = squeeze(std(pl.data_alpha,1,4))./sqrt(numel(pl.sub2plot));
pl.semdata_alpha_bc = squeeze(std(pl.data_alpha_bc,1,4))./sqrt(numel(pl.sub2plot));

figure;
set(gcf,'Position',[100 100 500 300],'PaperPositionMode','auto')
h.pl = {}; h.plsem=[];  h.plm = []; h.pls = []; h.plst = []; pl.con_label = [];
for i_con = 1:2
    % data index
    pl.idx = pl.time2plot_i(1):pl.time2plot_i(2);
    
    % plot SEM as boundary
    % create data
    pl.xconf = [TFA.time(pl.idx) TFA.time(pl.idx(end:-1:1))] ;
    pl.yconf = [pl.mdata_alpha(:,i_con)+pl.semdata_alpha(:,i_con); ...
        pl.mdata_alpha(end:-1:1,i_con)-pl.semdata_alpha(end:-1:1,i_con)];
    % plot
    h.plsem{ i_con} = fill(pl.xconf,pl.yconf',pl.col{i_con}./255,'EdgeColor','none','FaceAlpha',0.3);
    hold on
    
    % plot mean lines
    h.plm{ i_con}=plot(TFA.time(pl.idx), pl.mdata_alpha(:,i_con),'Color',pl.col{i_con}./255,'LineStyle',pl.line{t.idx},'LineWidth',2);
end
xlim(pl.time2plot)
xlabel('time in ms')
ylabel('amplitude in \muV/cm^3')
legend([h.plm{:}],pl.lab_con,'Location','SouthOutside','Orientation','horizontal')
legend('boxoff')
grid on
title(sprintf('alpha band amplitude | [%1.1f %1.1f]Hz',pl.alpha_range))

figure;
set(gcf,'Position',[100 100 500 300],'PaperPositionMode','auto')
h.pl = {}; h.plsem=[];  h.plm = []; h.pls = []; h.plst = []; pl.con_label = [];
for i_con = 1:2
    % data index
    pl.idx = pl.time2plot_i(1):pl.time2plot_i(2);
    
    % plot SEM as boundary
    % create data
    pl.xconf = [TFA.time(pl.idx) TFA.time(pl.idx(end:-1:1))] ;
    pl.yconf = [pl.mdata_alpha_bc(:,i_con)+pl.semdata_alpha_bc(:,i_con); ...
        pl.mdata_alpha_bc(end:-1:1,i_con)-pl.semdata_alpha_bc(end:-1:1,i_con)];
    % plot
    h.plsem{ i_con} = fill(pl.xconf,pl.yconf',pl.col{i_con}./255,'EdgeColor','none','FaceAlpha',0.3);
    hold on
    
    % plot mean lines
    h.plm{ i_con}=plot(TFA.time(pl.idx), pl.mdata_alpha_bc(:,i_con),'Color',pl.col{i_con}./255,'LineStyle',pl.line{t.idx},'LineWidth',2);
end
xlim(pl.time2plot)
xlabel('time in ms')
ylabel('modulation in %')
legend([h.plm{:}],pl.lab_con,'Location','SouthOutside','Orientation','horizontal')
legend('boxoff')
grid on
title(sprintf('alpha band modulation | [%1.1f %1.1f]Hz',pl.alpha_range))

sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_FShiftAlpha\figures\';
sav.filenames = 'TFA_alpha_mod_lineplotr_latelectclust';
% print(gcf, fullfile(sav.pathout,sav.filenames),'-dpng','-r600')
% saveas(gcf,fullfile(sav.pathout,sav.filenames),'fig')
% print(gcf,fullfile(sav.pathout,sav.filenames),'-depsc2', '-painters','-r300')

%% actual plotting data | TFA Grand Mean timecourse
% plotting parameters
% pl.elec2plot = {'Oz';'Iz'};
% pl.elec2plot = {'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';% vis alpha II
pl.elec2plot = {'P5';'P7';'P9';'P6';'P8';'P10';'PO3';'PO7';'PO4';'PO8';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';
% cluster analysis
% pl.elec2plot = {TFA.electrodes(1:64).labels}';
% pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) startsWith({TFA(1).electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.flims= TFA.frequency([1 end]);
pl.flims_i=dsearchn(TFA.frequency', pl.flims');

pl.xlims=[-500 1800]; % index time 2 plot
pl.xlims_i = dsearchn(TFA.time', pl.xlims');

pl.base = F.TFA.baseline;
pl.base_i = dsearchn(TFA.time', pl.base');

pl.sub2plot = F.Subs2use;
% pl.sub2plot = [7];

pl.data_ind = []; pl.data_evo = []; pl.data_ind_bc = []; pl.data_evo_bc = [];
pl.conlabel = {'attend RDK1';'attend RDK2'};
pl.RDKlabel = {'RDK1';'RDK2';'RDK3'};
% pl.timelabel = cellfun(@(x) vararg2str(x),TFA.ffttimewin(pl.time2plot),'UniformOutput',false);


pl.data_ind = squeeze(mean(mean(TFA.data_ind(pl.flims_i(1):pl.flims_i(2),pl.xlims_i(1):pl.xlims_i(2),pl.elec2plot_i,:,pl.sub2plot),3),4));
pl.data_evo = squeeze(mean(mean(TFA.data_evo(pl.flims_i(1):pl.flims_i(2),pl.xlims_i(1):pl.xlims_i(2),pl.elec2plot_i,:,pl.sub2plot),3),4));
pl.data_ind_bc = 100*(bsxfun(@rdivide, pl.data_ind, ...
    mean(squeeze(mean(mean(TFA.data_ind(pl.flims_i(1):pl.flims_i(2),pl.base_i(1):pl.base_i(2),pl.elec2plot_i,:,pl.sub2plot),3),4)),2))-1);
pl.data_evo_bc = 100*(bsxfun(@rdivide, pl.data_evo, ...
    mean(squeeze(mean(mean(TFA.data_evo(pl.flims_i(1):pl.flims_i(2),pl.base_i(1):pl.base_i(2),pl.elec2plot_i,:,pl.sub2plot),3),4)),2))-1);

figure;
set(gcf,'Position',[100 100 800 600],'PaperPositionMode','auto')
subplot(2,1,1)
pl.data = mean(pl.data_ind,3); pl.clims1=[0 max(pl.data(:))];
imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),TFA.frequency(pl.flims_i(1):pl.flims_i(2)),pl.data,pl.clims1)
colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
set(gca,'YDir','normal')
title(sprintf('Grand Mean TFA-amplitude | raw | induced\n for channel [%s]', vararg2str(pl.elec2plot)), ...
    'FontSize',8,'Interpreter','none')
xlabel('time in ms')
ylabel('frequency in Hz')
set(gca,'FontSize',8)
% set(gca,'ColorScale','log')
cb = colorbar();

subplot(2,1,2)
pl.data = mean(pl.data_evo,3); pl.clims1=[0 max(pl.data(:))];
imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),TFA.frequency(pl.flims_i(1):pl.flims_i(2)),pl.data,pl.clims1)
colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
set(gca,'YDir','normal')
title(sprintf('Grand Mean TFA-amplitude | raw | evoked\n for channel [%s]',vararg2str(pl.elec2plot)), ...
    'FontSize',8,'Interpreter','none')
xlabel('time in ms')
ylabel('frequency in Hz')
set(gca,'FontSize',8)
% set(gca,'ColorScale','log')
cb = colorbar();

% draw topography with electrode positions
h.a1 = axes('position',[0.85 0.45 0.14 0.14],'Visible','off');
topoplot(find(pl.elec2plot_i),TFA.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',4,1});


figure;
set(gcf,'Position',[100 100 800 600],'PaperPositionMode','auto')
subplot(2,1,1)
pl.data = mean(pl.data_ind_bc,3); pl.clims1=[-1 1]* max(abs(pl.data(:)));
imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),TFA.frequency(pl.flims_i(1):pl.flims_i(2)),pl.data,pl.clims1)
colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
set(gca,'YDir','normal')
title(sprintf('Grand Mean TFA-amplitude | baseline corrected | induced\n for channel [%s]',vararg2str(pl.elec2plot)), ...
    'FontSize',8,'Interpreter','none')
xlabel('time in ms')
ylabel('frequency in Hz')
set(gca,'FontSize',8)
% set(gca,'ColorScale','log')
cb = colorbar();

subplot(2,1,2)
pl.data = mean(pl.data_evo_bc,3); pl.clims1=[-1 1]* max(abs(pl.data(:)));
imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),TFA.frequency(pl.flims_i(1):pl.flims_i(2)),pl.data,pl.clims1)
colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
set(gca,'YDir','normal')
title(sprintf('Grand Mean TFA-amplitude | raw | evoked\n for channel [%s]',vararg2str(pl.elec2plot)), ...
    'FontSize',8,'Interpreter','none')
xlabel('time in ms')
ylabel('frequency in Hz')
set(gca,'FontSize',8)
% set(gca,'ColorScale','log')
cb = colorbar();

% draw topography with electrode positions
h.a1 = axes('position',[0.85 0.45 0.14 0.14],'Visible','off');
topoplot(find(pl.elec2plot_i),TFA.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',4,1});


%% actual plotting data | TFA timecourse
% plotting parameters
% pl.elec2plot = {'Oz';'Iz'};
% pl.elec2plot = {'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';% vis alpha II
pl.elec2plot = {'P5';'P7';'P9';'P6';'P8';'P10';'PO3';'PO7';'PO4';'PO8';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';
% cluster analysis
% pl.elec2plot = {TFA.electrodes(1:64).labels}';
% pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) startsWith({TFA(1).electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.freqrange=[-5 5];

pl.xlims=[-500 1800]; % index time 2 plot
pl.xlims_i = dsearchn(TFA.time', pl.xlims');

pl.base = F.TFA.baseline;
pl.base_i = dsearchn(TFA.time', pl.base');

pl.sub2plot = F.Subs2use;
% pl.sub2plot = [7];

pl.data_ind = []; pl.data_evo = []; pl.data_ind_bc = []; pl.data_evo_bc = [];
pl.conlabel = {'attend RDK1';'attend RDK2'};
pl.RDKlabel = {'RDK1';'RDK2';'RDK3'};
% pl.timelabel = cellfun(@(x) vararg2str(x),TFA.ffttimewin(pl.time2plot),'UniformOutput',false);


% extract data
for i_sub = 1:numel(pl.sub2plot)
    t.fidx = cell2mat(arrayfun(@(x) dsearchn(TFA.frequency', x+pl.freqrange'),[TFA.RDK(i_sub).RDK.freq],'UniformOutput',false));
    pl.freqlabel = TFA.frequency(t.fidx(1,1):t.fidx(2,1))-TFA.RDK(i_sub).RDK(1).freq;
    for i_RDK = 1:size(t.fidx,2)
        % raw
        pl.data_ind(:,:,:,i_RDK,i_sub) = ...
            squeeze(mean(TFA.data_ind(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.xlims_i(1):pl.xlims_i(2),pl.elec2plot_i,:,i_sub),3));
        pl.data_evo(:,:,:,i_RDK,i_sub) = ...
            squeeze(mean(TFA.data_evo(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.xlims_i(1):pl.xlims_i(2),pl.elec2plot_i,:,i_sub),3));
        % baseline corrected
        pl.data_ind_bc(:,:,:,i_RDK,i_sub) = ...
            100*(bsxfun(@rdivide, pl.data_ind(:,:,:,i_RDK,i_sub), ...
            mean(squeeze(mean(TFA.data_ind(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.base_i(1):pl.base_i(2),pl.elec2plot_i,:,i_sub),3)),2))-1);
        pl.data_evo_bc(:,:,:,i_RDK,i_sub) = ...
            100*(bsxfun(@rdivide, pl.data_evo(:,:,:,i_RDK,i_sub), ...
            mean(squeeze(mean(TFA.data_evo(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.base_i(1):pl.base_i(2),pl.elec2plot_i,:,i_sub),3)),2))-1);
    end   
end

% collapse across RDKs
t.data = pl.data_ind(:,:,:,1,:); t.data(:,:,:,:,:,2)= pl.data_ind(:,:,[2 1],2,:);
pl.data_ind_coll = squeeze(mean(t.data,6));
pl.data_ind_coll(:,:,3,:)=mean(pl.data_ind(:,:,:,3,:),3);
t.data = pl.data_evo(:,:,:,1,:); t.data(:,:,:,:,:,2)= pl.data_evo(:,:,[2 1],2,:);
pl.data_evo_coll = squeeze(mean(t.data,6));
pl.data_evo_coll(:,:,3,:)=mean(pl.data_evo(:,:,:,3,:),3);
t.data = pl.data_ind_bc(:,:,:,1,:); t.data(:,:,:,:,:,2)= pl.data_ind_bc(:,:,[2 1],2,:);
pl.data_ind_bc_coll = squeeze(mean(t.data,6));
pl.data_ind_bc_coll(:,:,3,:)=mean(pl.data_ind_bc(:,:,:,3,:),3);
t.data = pl.data_evo_bc(:,:,:,1,:); t.data(:,:,:,:,:,2)= pl.data_evo_bc(:,:,[2 1],2,:);
pl.data_evo_bc_coll = squeeze(mean(t.data,6));
pl.data_evo_bc_coll(:,:,3,:)=mean(pl.data_evo_bc(:,:,:,3,:),3);
pl.conlabel2 = {'attended';'unattended';'irrelevant'};


% raw
pl.data = mean(pl.data_ind_coll,4);
pl.clims1=[0 1].*squeeze(max(max(pl.data)));
pl.clims1=[0 1].*repmat(max(max(max(pl.data))),3,1);
figure;
set(gcf,'Position',[100 100 1500 600],'PaperPositionMode','auto')
for i_con = 1:size(pl.data,3)
    subplot(size(pl.data,3),2,1+2*(i_con-1))
    imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),pl.freqlabel,pl.data(:,:,i_con),pl.clims1(i_con,:))
    colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    set(gca,'YDir','normal')
    if i_con == 1
        title(sprintf('TFA-amplitude | raw | induced | %s\n for channel [%s]',pl.conlabel2{i_con}, vararg2str(pl.elec2plot)), ...
            'FontSize',8,'Interpreter','none')
    end
    
    if i_con ~= 3
%         set(gca,'XTickLabel',[])
    else
        xlabel('time in ms')
    end
    ylabel({pl.conlabel2{i_con};'delta f_R_D_K in Hz'})
    xlim(pl.xlims)
    set(gca,'FontSize',8)
    % set(gca,'ColorScale','log')
    cb = colorbar();
    
end

pl.data = mean(pl.data_evo_coll,4);
pl.clims1=[0 1].*squeeze(max(max(pl.data)));
pl.clims1=[0 1].*repmat(max(max(max(pl.data))),3,1);
for i_con = 1:size(pl.data,3)
    subplot(size(pl.data,3),2,2+2*(i_con-1))
    imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),pl.freqlabel,pl.data(:,:,i_con),pl.clims1(i_con,:))
    colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    set(gca,'YDir','normal')
    if i_con == 1
        title(sprintf('TFA-amplitude | baseline corrected | evoked \n for channel [%s]', vararg2str(pl.elec2plot)), ...
            'FontSize',8,'Interpreter','none')
    end
    
    if i_con ~= 3
%         set(gca,'XTickLabel',[])
    else
        xlabel('time in ms')
    end
    ylabel({pl.conlabel2{i_con};'delta f_R_D_K in Hz'})
    xlim(pl.xlims)
    set(gca,'FontSize',8)
    % set(gca,'ColorScale','log')
    cb = colorbar();
    
end


% draw topography with electrode positions
h.a1 = axes('position',[0.425 0.75 0.15 0.15],'Visible','off');
% topoplot(find(pl.elec2plot_i),TFA(1).electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
%     'emarker2',{find(pl.elec2plot_i),'o','r',4,1});
topoplot(find(pl.elec2plot_i),TFA(1).electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o',[255 133 4]./255,4,1});


% baseline corrected
pl.data = mean(pl.data_ind_bc_coll,4);
pl.clims1=[-1 1].*squeeze(max(max(abs(pl.data))));
pl.clims1=[-1 1].*repmat(max(max(max(abs(pl.data)))),3,1);
figure;
set(gcf,'Position',[100 100 1500 600],'PaperPositionMode','auto')
for i_con = 1:size(pl.data,3)
    subplot(size(pl.data,3),2,1+2*(i_con-1))
    imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),pl.freqlabel,pl.data(:,:,i_con),pl.clims1(i_con,:))
    colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    set(gca,'YDir','normal')
    if i_con == 1
        title(sprintf('TFA-amplitude | baseline corrected | induced | %s\n for channel [%s]',pl.conlabel2{i_con}, vararg2str(pl.elec2plot)), ...
            'FontSize',8,'Interpreter','none')
    end
    
    if i_con ~= 3
%         set(gca,'XTickLabel',[])
    else
        xlabel('time in ms')
    end
    ylabel({pl.conlabel2{i_con};'delta f_R_D_K in Hz'})
    xlim(pl.xlims)
    set(gca,'FontSize',8)
    % set(gca,'ColorScale','log')
    cb = colorbar();
    
end

pl.data = mean(pl.data_evo_bc_coll,4);
pl.clims1=[-1 1].*squeeze(max(max(abs(pl.data))));
pl.clims1=[-1 1].*repmat(max(max(max(abs(pl.data)))),3,1);
for i_con = 1:size(pl.data,3)
    subplot(size(pl.data,3),2,2+2*(i_con-1))
    imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),pl.freqlabel,pl.data(:,:,i_con),pl.clims1(i_con,:))
    colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    set(gca,'YDir','normal')
    if i_con == 1
        title(sprintf('TFA-amplitude | raw | evoked \n for channel [%s]', vararg2str(pl.elec2plot)), ...
            'FontSize',8,'Interpreter','none')
    end
    
    if i_con ~= 3
%         set(gca,'XTickLabel',[])
    else
        xlabel('time in ms')
    end
    ylabel({pl.conlabel2{i_con};'frequency in Hz'})
    xlim(pl.xlims)
    set(gca,'FontSize',8)
    % set(gca,'ColorScale','log')
    cb = colorbar();
    
end


% draw topography with electrode positions
h.a1 = axes('position',[0.425 0.75 0.15 0.15],'Visible','off');
% topoplot(find(pl.elec2plot_i),TFA(1).electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
%     'emarker2',{find(pl.elec2plot_i),'o','r',4,1});
topoplot(find(pl.elec2plot_i),TFA(1).electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o',[255 133 4]./255,4,1});






