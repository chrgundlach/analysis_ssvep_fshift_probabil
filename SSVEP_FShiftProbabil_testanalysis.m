%% scrript to do first preliminary analysis
% includes FFT so far
clearvars


%% parameters
p.pathin=           'E:\work\data\SSVEP_FShift_Probabil\eeg\epoch\';
p.subs=             cellfun(@(x) sprintf('%02.0f',x),num2cell(1:40),'UniformOutput', false)';
% p.subs2use=         [9 10 11 12 13];
p.subs2use =        [1 3 4 5 6 7 9 10 11 12];
p.freqs=            [18 21 24];

p.conds=            {'valid | target RDK1'; 'valid | target RDK2'; ...
    'invalid | target RDK2'; 'invalid | target RDK2';...
    'neutral | target RDK1'; 'neutral | target RDK2'};
p.conds2=            {'valid | attend RDK1'; 'valid | attend RDK2'; ...
    'invalid | attend RDK2'; 'invalid | attend RDK1';...
    'neutral | attend RDK1'; 'neutral | attend RDK2'};

p.conds_trig=       {[10] [20] [30] [40] [50] [60]}; % regular
% p.evnt_trig=        [81 82 86 87];
% p.epoch_evnt=       [-2 1];

FFT.fftres=         10000; % frequency resolution
FFT.time=           {[-1000 0]; [-500 500]; [0 1000]; [500 1500]; [1000 2000]; [1500 2500]; [2000 3000]; [-1000 3000]};


%% loop across subjects [extracting data]
for i_sub = 1:numel(p.subs2use)
    % load file
    EEG = pop_loadset(sprintf('VP%s_e.set',p.subs{p.subs2use(i_sub)}), p.pathin);
    % pop_eegplot(EEG,1,1,1)
%     figure; plot(EEG.times, mean(EEG.data(:,:,:),3))
    
    if i_sub == 1
        FFT.data=nan(EEG.nbchan,FFT.fftres,numel(p.subs2use),numel(p.conds_trig),numel(FFT.time),2);
        FFT.data2 = nan(EEG.nbchan,FFT.fftres,numel(p.subs2use),numel(p.conds_trig),1,2);
    end
    
    %% loop for conditions
    for i_con = 1:numel(p.conds_trig)
        % select only portion of trials
        EEGtmp = pop_selectevent( EEG, 'type',p.conds_trig{i_con},'deleteevents','off','deleteepochs','on','invertepochs','off');
        % pop_eegplot(EEGtmp,1,1,1)
%         figure; plot(EEG.times, mean(EEGtmp.data(:,:,:),3))
        param.trialsused(i_con,i_sub)=EEGtmp.trials;
        
        %% analysis cue locked
        % preallocate memory
        t.fftdata = nan(EEG.nbchan,FFT.fftres,numel(FFT.time),2);
        for i_time = 1:numel(FFT.time)
            % index time points
            t.index = eeg_time2points(FFT.time{i_time}(1),EEGtmp.times):eeg_time2points(FFT.time{i_time}(2),EEGtmp.times)-1;
            % induced data
            % loop across trials
            t.fftdata_ind = nan(EEG.nbchan,FFT.fftres,EEGtmp.trials);
            for i_tr = 1:EEGtmp.trials
                t.data = detrend(EEGtmp.data(:,t.index,i_tr)')';
                t.fftdata_ind(:,:,i_tr) = abs(fft(t.data,FFT.fftres,2))*2/size(t.data,2);
            end
            % average across trials
            t.fftdata(:,:,i_time,1)= mean(t.fftdata_ind,3);
            
             % evoked data
            t.data_m = detrend(mean(EEGtmp.data(:,t.index,:),3)')';
            t.fftdata(:,:,i_time,2)= abs(fft(t.data_m,FFT.fftres,2))*2/size(t.data_m,2);           
        end
        
        % average across time windows...not
%         FFT.data(:,:,i_sub,i_con,:)=mean(t.fftdata,3);
        FFT.data(:,:,i_sub,i_con,:,:)=t.fftdata;
        
%        
        
    end
end
FFT.freqs = ((0:size(FFT.data,2)-1)/size(FFT.data,2)) * EEGtmp.srate;



%% plot grand mean data for SSVEPs/alpha
% pl.freq2plot=[17.9 18.1]; p.freqlabel='SSVEP1_left_stim';
% pl.freq2plot=[20.9 21.1]; p.freqlabel='SSVEP2_right_stim';
pl.freq2plot=[23.9 24.1]; p.freqlabel='SSVEP1_left_stim';
% pl.freq2plot=[8 12]; p.freqlabel='alpha';
pl.freqindex = eeg_time2points(pl.freq2plot,FFT.freqs);
% pl.dat2plot ='induced'; pl.dat2plot_i = 1;
pl.dat2plot ='evoked'; pl.dat2plot_i = 2;

pl.subs2use = 1:numel(p.subs2use);
pl.clims = []; % absmax if pl.clims = [];

pl.time2plot = [8];

figure;


pl.data = mean(FFT.data(:,pl.freqindex(1):pl.freqindex(2),pl.subs2use,:,pl.time2plot,pl.dat2plot_i),[2 3 4 5]);
if isempty(pl.clims)
    h.tp=topoplot( pl.data, EEGtmp.chanlocs(1:64), ...
        'numcontour', 0, 'conv', 'on', 'maplimits',[0 max(pl.data)],'whitebk','on','colormap',fake_parula);
else
    h.tp=topoplot( pl.data, EEGtmp.chanlocs(1:64), ...
        'numcontour', 0, 'conv', 'on', 'maplimits',pl.clims,'whitebk','on','colormap',fake_parula);
end
    
title(sprintf('Grand mean topo | all conditions and time windows\n[%1.2f %1.2f] Hz | %s | %s',...
    pl.freq2plot, pl.dat2plot, p.freqlabel),'FontSize',8,'Interpreter','none')
% colormap((plasma)) % flipud(bone)
colorbar

set(gcf,'Position',[100 100 400 330],'PaperPositionMode','auto')
% print(fullfile(p.pathout_fft,'figures',sprintf('TOPO_Grand_mean_%s_%s',pl.dat2plot,p.freqlabel)),'-djpeg','-r300')
% saveas(gcf,fullfile(p.pathout_fft,'figures',sprintf('TOPO_Grand_mean_%s_%s',pl.dat2plot,p.freqlabel)),'fig')


%% plot fft
% pl.elec2plot = {'POz';'Oz'}; % periphery (large cluster)
% pl.elec2plot = {'P04';'POz';'PO8';'O2';'Oz';'I2';'Iz'}; % for left stimuls 
% pl.elec2plot = {'P03';'POz';'PO7';'O1';'Oz';'I1';'Iz'}; % for right stimuls
pl.elec2plot = {'O1';'Oz';'O2';'I1';'Iz';'I2'}; % periphery (large cluster)
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({EEG.chanlocs.labels},x), lower(pl.elec2plot), 'UniformOutput',false)),1));
pl.xlims = [0 30];
pl.freqcolor_tr={[1 0 0];[0 0 1];[0 1 0]};
pl.freqcolor_n={'c';'b';'m'};
% pl.dat2plot ='induced'; pl.dat2plot_i = 1;
pl.dat2plot ='evoked'; pl.dat2plot_i = 2;

pl.time2plot = [8];

figure;
pl.data=squeeze(mean(FFT.data(pl.elec2plot_i,:,:,:,pl.time2plot,pl.dat2plot_i),[1, 3, 5]));
h.p=plot(FFT.freqs,pl.data);
xlim(pl.xlims)
xlabel('frequency in Hz')
ylabel('amplitude in \muV')
grid on
ylim([0 1.1*max(max(pl.data(eeg_time2points(pl.xlims(1),FFT.freqs):eeg_time2points(pl.xlims(2),FFT.freqs),:)))])

for i_l = 1:numel(p.freqs) 
    try vline(p.freqs(i_l),pl.freqcolor_tr{i_l}); 
    catch; vline(p.freqs(i_l),pl.freqcolor_n{i_l}); 
    end 
end
title(sprintf('amplitude spectra | at %s | [%1.0f %1.0f]ms', ...
    vararg2str(pl.elec2plot),FFT.time{pl.time2plot(1)}(1),FFT.time{pl.time2plot(end)}(2)))
legend(p.conds)

figure;
pl.data=squeeze(mean(FFT.data(pl.elec2plot_i,:,:,:,pl.time2plot,pl.dat2plot_i),[1,4,5]));
h.p=plot(FFT.freqs,pl.data,'Color',[0.4 0.4 0.4]);
hold on; h.p2=plot(FFT.freqs,mean(pl.data,2),'Color',[0.8 0.4 0.4],'LineWidth',1);
xlim(pl.xlims)
xlabel('frequency in Hz')
ylabel('amplitude in \muV')
grid on
ylim([0 1.1*max(max(pl.data(eeg_time2points(pl.xlims(1),FFT.freqs):eeg_time2points(pl.xlims(2),FFT.freqs))))])

for i_l = 1:3 
    try vline(p.freqs(i_l),pl.freqcolor_tr{i_l}); 
    catch; vline(p.freqs(i_l),pl.freqcolor_n{i_l}); 
    end 
end
title(sprintf('grand mean amplitude spectrum | at %s | [%1.0f %1.0f]ms',vararg2str(pl.elec2plot),FFT.time{1}(1),FFT.time{end}(2)))



