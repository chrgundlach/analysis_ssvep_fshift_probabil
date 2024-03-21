%% plot TFA images
clearvars
F.PathInEEG             = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2023_FShift_Probabil\eeg\erp'; 


F.Subs                  = arrayfun(@(x) sprintf('%02.0f',x),1:40,'UniformOutput',false)';
F.Subs2use              = [1 3 4 5 6 7 9 10 11 13 15 18 20:28 30:34]; 
                        % 2 and 8 are excluded as the didn't do the task properly, sub 11 has potentially low number of trials
                        % for subject 12, 14, 39: eeg and behavior data don't match



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
    fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%s_target_erp.mat ||\n',...
        i_sub,numel(F.Subs2use),F.PathInEEG,F.Subs{F.Subs2use(i_sub)})
    
    temp.erp = open(fullfile(F.PathInEEG,sprintf('VP%s_target_erp.mat',F.Subs{F.Subs2use(i_sub)})));
    
    % preallocate memory
    if i_sub == 1
        EP.filtdata = repmat({[]},1,numel(F.Subs2use));
        EP.behavior = repmat({[]},1,numel(F.Subs2use));
        EP.RDK = repmat({[]},1,numel(F.Subs2use));
        EP.params = temp.erp.EP.parameters;
        EP.electrodes = temp.erp.EP.EEG_fep.chanlocs;
        EP.time = temp.erp.EP.EEG_fep.times;
        EP.srate = temp.erp.EP.EEG_fep.srate;
    end
    
    % assign data
    EP.filtdata{i_sub} = temp.erp.EP.EEG_fep;
    EP.behavior{i_sub} = temp.erp.EP.behavior;
    EP.RDK{i_sub} = temp.erp.EP.RDK;
    
    clear temp
    
end

pl.plegend = {'1';'0.5';'0.25';'0.1';'0.05';'0.01';'0.001';'0.0001';'0.00001'};
pl.pcorrect = [0 abs(log10([0.5 0.25 0.1 0.05 0.01 0.001 0.0001 0.00001]))];



%% plot ERPs exploratively in topo array
    
% loop through conditions defined in contrasts
pl.con_contrast = {... % contrasts by 1st dim; averaged across second dim
    'trial_timing_type', {{'regular'}}; ...
    'cue_validity_label', {{'valid'};{'neutral'};{'invalid'}}};
pl.sub2plot = 1:numel(F.Subs2use);


% preallocate data
pl.dat2plot = nan([size(EP.electrodes,2), ...           % 1st dim channels
    size(EP.time,2), ...                                % 2nd dim time
    cellfun(@(x) size(x,1), pl.con_contrast(:,2))', ... % 3rd to n dim = levels defind in pl.con_contrast
    numel(pl.sub2plot)]);                               % n+1 dim participants

% extract data across contrasts
% define how to loop through contrasts
t.cont = cellfun(@(x) 1:size(x,1), pl.con_contrast(:,2), 'UniformOutput', false);
t.contidx = combvec(t.cont{:});
% loop across participants
for i_sub = 1:numel(pl.sub2plot)
    t.behavior = EP.behavior{pl.sub2plot(i_sub)};
    t.epdata = EP.filtdata{pl.sub2plot(i_sub)}.data;

    % loop through contrasts
    for i_cont = 1:size(t.contidx,2)
        % define evaluation syntax
                
        % preallocate logical idx for trials
        t.idx = false(size(t.contidx,1),numel(t.behavior ));
        for i_contstep = 1:size(t.contidx,1)
            % index step
            t.idx(i_contstep,:) = ...
                any( ... % if it is more than one level, combine
                cell2mat( ...
                cellfun(@(x) strcmp({t.behavior.(pl.con_contrast{i_contstep,1})},x), ...
                pl.con_contrast{i_contstep,2}{t.contidx(i_contstep,i_cont)}, 'UniformOutput', false)' ...
                ) ...
                ,1);
        end
        % combine all contrast conditions
        t.idx_all = all(t.idx,1);
        
        % average across trials
        t.text = sprintf("%1.0f,",t.contidx(:,i_cont));
        t.evaltext = sprintf("pl.dat2plot(:,:,%si_sub) = mean(t.epdata(:,:,t.idx_all),3);",t.text);
        eval(t.evaltext)
    end
end

% reshape pl.data
pl.dat2plot_rs = mean(squeeze(pl.dat2plot),4);

pl.con_label = [pl.con_contrast{2,2}{:}];

plot_data_topoarray(EP.electrodes, pl.dat2plot_rs,'ERP','times',EP.time,'conds',pl.con_label)


%% plot ERPs exploratively for specified electrodes | cue validity

% loop through conditions defined in contrasts
pl.con_contrast = {... % contrasts by 1st dim; averaged across second dim
    'trial_timing_type', {{'regular'}}; ...
    'event_response_type', {{'hit'}}; ...
%     'evnt_type_label', {{'chroma+'}}; ...
%     'event_response_type', {{'hit','FA','error','miss'}}; ...
    'cue_validity_label', {{'valid'};{'neutral'};{'invalid'}}};
pl.sub2plot = 1:numel(F.Subs2use);
% pl.elec2plot = {'P8';'PO8';'P10';'P7';'PO7';'P9'}; % for P1 component lateral !
% pl.elec2plot = {'P7';'PO7';'P9';'O1';'I1';'Oz'; 'Iz';'O2';'I2';'P8';'PO8';'P10';}; % for N2 component posterior!
pl.elec2plot = {'P3';'P1';'Pz';'P4';'P2';'POz';'PO3';'PO4'}; % for P300 component centro-parietal!
% pl.elec2plot = {'POz';'Oz';'O1';'O2';'Iz'}; % for N1 component centro-parietal
% pl.elec2plot = {'P6';'P8';'PO8';'P10';'P5';'P7';'PO7';'P9'}; % for N2 component lateral
% pl.elec2plot = {'POz'}; % early N2 SN?
% pl.elec2plot = {'CPz';'Cz'}; % early N2 SN?
pl.elec2plot = {'FCz';'Fz'}; % early frontal?

pl.con_label = {'valid';'neutral';'invalid'};

pl.elec2plot_i = ...
    any(cell2mat(cellfun(@(x) strcmp({EP.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1);

pl.time2plot = [-100 500]; % time in ms
pl.time2plot = [-100 700]; % time in ms
pl.time2plot_i = dsearchn(EP.time', pl.time2plot');

pl.concols = num2cell([25 138 130; 41 60 74; 241 131 26]'./255,1);


% preallocate data
pl.dat2plot = nan([sum(pl.elec2plot_i), ...             % 1st dim channels
    size(EP.time,2), ...                                % 2nd dim time
    cellfun(@(x) size(x,1), pl.con_contrast(:,2))', ... % 3rd to n dim = levels defind in pl.con_contrast
    numel(pl.sub2plot)]);                               % n+1 dim participants

% extract data across contrasts
% define how to loop through contrasts
t.cont = cellfun(@(x) 1:size(x,1), pl.con_contrast(:,2), 'UniformOutput', false);
t.contidx = combvec(t.cont{:});
% loop across participants
for i_sub = 1:numel(pl.sub2plot)
    t.behavior = EP.behavior{pl.sub2plot(i_sub)};
    t.epdata = EP.filtdata{pl.sub2plot(i_sub)}.data;

    % loop through contrasts
    for i_cont = 1:size(t.contidx,2)
        % define evaluation syntax
                
        % preallocate logical idx for trials
        t.idx = false(size(t.contidx,1),numel(t.behavior ));
        for i_contstep = 1:size(t.contidx,1)
            % index step
            t.idx(i_contstep,:) = ...
                any( ... % if it is more than one level, combine
                cell2mat( ...
                cellfun(@(x) strcmp({t.behavior.(pl.con_contrast{i_contstep,1})},x), ...
                pl.con_contrast{i_contstep,2}{t.contidx(i_contstep,i_cont)}, 'UniformOutput', false)' ...
                ) ...
                ,1);
        end
        % combine all contrast conditions
        t.idx_all = all(t.idx,1);
        
        % average across trials
        t.text = sprintf("%1.0f,",t.contidx(:,i_cont));
        t.evaltext = sprintf("pl.dat2plot(:,:,%si_sub) = mean(t.epdata(pl.elec2plot_i,:,t.idx_all),3);",t.text);
        eval(t.evaltext)
    end
end

% reshape pl.data
pl.dat2plot_rs = squeeze(mean(pl.dat2plot,1));

pl.data_m = mean(pl.dat2plot_rs,3);
pl.data_sem = std(pl.dat2plot_rs,1,3)./sqrt(numel(pl.sub2plot));



figure;
set(gcf,'Position',[100 100 600 300],'PaperPositionMode','auto')
pl.idx = pl.time2plot_i(1):pl.time2plot_i(2);
h.plsem=[];  h.plm = []; h.pls = []; h.plst = [];

for i_con = 1:size(pl.data_m,2)
    % plot SEM as boundary
    % create data
    pl.xconf = [EP.time(pl.idx) EP.time(pl.idx(end:-1:1))] ;
    pl.yconf = [pl.data_m(pl.idx,i_con)' + pl.data_sem(pl.idx,i_con)' ...
        pl.data_m(pl.idx(end:-1:1),i_con)' - pl.data_sem(pl.idx(end:-1:1),i_con)'];
    % plot
    h.plsem{i_con} = fill(pl.xconf,pl.yconf,pl.concols{i_con}','EdgeColor','none','FaceAlpha',0.3);
    hold on
    
    % plot mean lines
    h.plm{i_con}=...
        plot(EP.time(pl.idx), pl.data_m(pl.idx,i_con),'Color',pl.concols{i_con},'LineWidth',2);
end
xlim(pl.time2plot)
xlabel('time in ms')
ylabel('amplitude in \muV/cm²')
legend([h.plm{:}],pl.con_label,'Location','EastOutside','box','off')
grid on
set(gca, 'Ydir','reverse')


h.ax2 = axes('Position', [0.725, .63, .25, .25],'Visible','off');
topoplot(find(pl.elec2plot_i),EP.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',3,1});

set(gcf, 'Color', [1 1 1]);

%% plot ERPs exploratively for specified electrodes | hit vs miss

% loop through conditions defined in contrasts
pl.con_contrast = {... % contrasts by 1st dim; averaged across second dim
    'trial_timing_type', {{'regular'}}; ...
    'evnt_type_label', {{'chroma+'}}; ...
    'event_response_type', {{'hit'};{'FA','error','miss'}}};
pl.sub2plot = 1:numel(F.Subs2use);
% pl.elec2plot = {'P8';'PO8';'P10';'P7';'PO7';'P9'}; % for P1 component lateral !
% pl.elec2plot = {'P7';'PO7';'P9';'O1';'I1';'Oz'; 'Iz';'O2';'I2';'P8';'PO8';'P10';}; % for N2 component posterior!
pl.elec2plot = {'P3';'P1';'Pz';'P4';'P2';'POz';'PO3';'PO4'}; % for P300 component centro-parietal!
% pl.elec2plot = {'POz';'Oz';'O1';'O2';'Iz'}; % for N1 component centro-parietal
% pl.elec2plot = {'P6';'P8';'PO8';'P10';'P5';'P7';'PO7';'P9'}; % for N2 component lateral

% pl.elec2plot = {'POz'}; % early N2 SN?
% pl.elec2plot = {'CPz';'Cz'}; % early N2 SN?
pl.elec2plot = {'FCz';'Fz'}; % early frontal?


pl.elec2plot_i = ...
    any(cell2mat(cellfun(@(x) strcmp({EP.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1);

pl.time2plot = [-100 500]; % time in ms
pl.time2plot = [-100 700]; % time in ms
pl.time2plot_i = dsearchn(EP.time', pl.time2plot');

pl.concols = num2cell([63 63 240; 240 63 63]'./255,1);
pl.con_label = {'hit';'error+FA+miss'};


% preallocate data
pl.dat2plot = nan([sum(pl.elec2plot_i), ...             % 1st dim channels
    size(EP.time,2), ...                                % 2nd dim time
    cellfun(@(x) size(x,1), pl.con_contrast(:,2))', ... % 3rd to n dim = levels defind in pl.con_contrast
    numel(pl.sub2plot)]);                               % n+1 dim participants

% extract data across contrasts
% define how to loop through contrasts
t.cont = cellfun(@(x) 1:size(x,1), pl.con_contrast(:,2), 'UniformOutput', false);
t.contidx = combvec(t.cont{:});
% loop across participants
for i_sub = 1:numel(pl.sub2plot)
    t.behavior = EP.behavior{pl.sub2plot(i_sub)};
    t.epdata = EP.filtdata{pl.sub2plot(i_sub)}.data;

    % loop through contrasts
    for i_cont = 1:size(t.contidx,2)
        % define evaluation syntax
                
        % preallocate logical idx for trials
        t.idx = false(size(t.contidx,1),numel(t.behavior ));
        for i_contstep = 1:size(t.contidx,1)
            % index step
            t.idx(i_contstep,:) = ...
                any( ... % if it is more than one level, combine
                cell2mat( ...
                cellfun(@(x) strcmp({t.behavior.(pl.con_contrast{i_contstep,1})},x), ...
                pl.con_contrast{i_contstep,2}{t.contidx(i_contstep,i_cont)}, 'UniformOutput', false)' ...
                ) ...
                ,1);
        end
        % combine all contrast conditions
        t.idx_all = all(t.idx,1);
        
        % average across trials
        t.text = sprintf("%1.0f,",t.contidx(:,i_cont));
        t.evaltext = sprintf("pl.dat2plot(:,:,%si_sub) = mean(t.epdata(pl.elec2plot_i,:,t.idx_all),3);",t.text);
        eval(t.evaltext)
    end
end

% reshape pl.data
pl.dat2plot_rs = squeeze(mean(pl.dat2plot,1));

pl.data_m = mean(pl.dat2plot_rs,3);
pl.data_sem = std(pl.dat2plot_rs,1,3)./sqrt(numel(pl.sub2plot));



figure;
set(gcf,'Position',[100 100 600 300],'PaperPositionMode','auto')
pl.idx = pl.time2plot_i(1):pl.time2plot_i(2);
h.plsem=[];  h.plm = []; h.pls = []; h.plst = [];

for i_con = 1:size(pl.data_m,2)
    % plot SEM as boundary
    % create data
    pl.xconf = [EP.time(pl.idx) EP.time(pl.idx(end:-1:1))] ;
    pl.yconf = [pl.data_m(pl.idx,i_con)' + pl.data_sem(pl.idx,i_con)' ...
        pl.data_m(pl.idx(end:-1:1),i_con)' - pl.data_sem(pl.idx(end:-1:1),i_con)'];
    % plot
    h.plsem{i_con} = fill(pl.xconf,pl.yconf,pl.concols{i_con}','EdgeColor','none','FaceAlpha',0.3);
    hold on
    
    % plot mean lines
    h.plm{i_con}=...
        plot(EP.time(pl.idx), pl.data_m(pl.idx,i_con),'Color',pl.concols{i_con},'LineWidth',2);
end
xlim(pl.time2plot)
xlabel('time in ms')
ylabel('amplitude in \muV/cm²')
legend([h.plm{:}],pl.con_label,'Location','EastOutside','box','off')
grid on
set(gca, 'Ydir','reverse')


h.ax2 = axes('Position', [0.725, .63, .25, .25],'Visible','off');
topoplot(find(pl.elec2plot_i),EP.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',3,1});

set(gcf, 'Color', [1 1 1]);

%% plot ERPs exploratively for specified electrodes | chroma+ or chroma-
% loop through conditions defined in contrasts
pl.con_contrast = {... % contrasts by 1st dim; averaged across second dim
    'trial_timing_type', {{'regular'}}; ...
    'event_response_type', {{'hit'}}; ...
%     'event_response_type', {{'hit','FA','error','miss'}}; ...
    'evnt_type_label', {{'chroma+'};{'chroma-'}}};
pl.sub2plot = 1:numel(F.Subs2use);
% pl.elec2plot = {'P8';'PO8';'P10';'P7';'PO7';'P9'}; % for P1 component lateral !
% pl.elec2plot = {'P7';'PO7';'P9';'O1';'I1';'Oz'; 'Iz';'O2';'I2';'P8';'PO8';'P10';}; % for N2 component posterior!
pl.elec2plot = {'P3';'P1';'Pz';'P4';'P2';'POz';'PO3';'PO4'}; % for P300 component centro-parietal!
% pl.elec2plot = {'POz';'Oz';'O1';'O2';'Iz'}; % for N1 component centro-parietal
% pl.elec2plot = {'P6';'P8';'PO8';'P10';'P5';'P7';'PO7';'P9'}; % for N2 component lateral

% pl.elec2plot = {'POz'}; % early N2 SN?
% pl.elec2plot = {'CPz';'Cz'}; % early N2 SN?
pl.elec2plot = {'FCz';'Fz'}; % early frontal?


pl.elec2plot_i = ...
    any(cell2mat(cellfun(@(x) strcmp({EP.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1);

pl.time2plot = [-100 500]; % time in ms
pl.time2plot = [-100 700]; % time in ms
pl.time2plot_i = dsearchn(EP.time', pl.time2plot');

pl.concols = num2cell([240 63 240; 100 63 100]'./255,1);
pl.con_label = {'chroma+';'chroma-'};


% preallocate data
pl.dat2plot = nan([sum(pl.elec2plot_i), ...             % 1st dim channels
    size(EP.time,2), ...                                % 2nd dim time
    cellfun(@(x) size(x,1), pl.con_contrast(:,2))', ... % 3rd to n dim = levels defind in pl.con_contrast
    numel(pl.sub2plot)]);                               % n+1 dim participants

% extract data across contrasts
% define how to loop through contrasts
t.cont = cellfun(@(x) 1:size(x,1), pl.con_contrast(:,2), 'UniformOutput', false);
t.contidx = combvec(t.cont{:});
% loop across participants
for i_sub = 1:numel(pl.sub2plot)
    t.behavior = EP.behavior{pl.sub2plot(i_sub)};
    t.epdata = EP.filtdata{pl.sub2plot(i_sub)}.data;

    % loop through contrasts
    for i_cont = 1:size(t.contidx,2)
        % define evaluation syntax
                
        % preallocate logical idx for trials
        t.idx = false(size(t.contidx,1),numel(t.behavior ));
        for i_contstep = 1:size(t.contidx,1)
            % index step
            t.idx(i_contstep,:) = ...
                any( ... % if it is more than one level, combine
                cell2mat( ...
                cellfun(@(x) strcmp({t.behavior.(pl.con_contrast{i_contstep,1})},x), ...
                pl.con_contrast{i_contstep,2}{t.contidx(i_contstep,i_cont)}, 'UniformOutput', false)' ...
                ) ...
                ,1);
        end
        % combine all contrast conditions
        t.idx_all = all(t.idx,1);
        
        % average across trials
        t.text = sprintf("%1.0f,",t.contidx(:,i_cont));
        t.evaltext = sprintf("pl.dat2plot(:,:,%si_sub) = mean(t.epdata(pl.elec2plot_i,:,t.idx_all),3);",t.text);
        eval(t.evaltext)
    end
end

% reshape pl.data
pl.dat2plot_rs = squeeze(mean(pl.dat2plot,1));

pl.data_m = mean(pl.dat2plot_rs,3);
pl.data_sem = std(pl.dat2plot_rs,1,3)./sqrt(numel(pl.sub2plot));



figure;
set(gcf,'Position',[100 100 600 300],'PaperPositionMode','auto')
pl.idx = pl.time2plot_i(1):pl.time2plot_i(2);
h.plsem=[];  h.plm = []; h.pls = []; h.plst = [];

for i_con = 1:size(pl.data_m,2)
    % plot SEM as boundary
    % create data
    pl.xconf = [EP.time(pl.idx) EP.time(pl.idx(end:-1:1))] ;
    pl.yconf = [pl.data_m(pl.idx,i_con)' + pl.data_sem(pl.idx,i_con)' ...
        pl.data_m(pl.idx(end:-1:1),i_con)' - pl.data_sem(pl.idx(end:-1:1),i_con)'];
    % plot
    h.plsem{i_con} = fill(pl.xconf,pl.yconf,pl.concols{i_con}','EdgeColor','none','FaceAlpha',0.3);
    hold on
    
    % plot mean lines
    h.plm{i_con}=...
        plot(EP.time(pl.idx), pl.data_m(pl.idx,i_con),'Color',pl.concols{i_con},'LineWidth',2);
end
xlim(pl.time2plot)
xlabel('time in ms')
ylabel('amplitude in \muV/cm²')
legend([h.plm{:}],pl.con_label,'Location','EastOutside','box','off')
grid on
set(gca, 'Ydir','reverse')


h.ax2 = axes('Position', [0.725, .63, .25, .25],'Visible','off');
topoplot(find(pl.elec2plot_i),EP.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',3,1});

set(gcf, 'Color', [1 1 1]);

%% plot topagraphies | cue validity


% loop through conditions defined in contrasts
pl.con_contrast = {... % contrasts by 1st dim; averaged across second dim
    'trial_timing_type', {{'regular'}}; ...
    'event_response_type', {{'hit'}}; ...
    'cue_validity_label', {{'valid'};{'neutral'};{'invalid'}}};
pl.sub2plot = 1:numel(F.Subs2use);

pl.time2plot = [160 190]; % time in ms
pl.time2plot = [160 250]; % time in ms
pl.time2plot = [100 150]; % time in ms

% pl.time2plot = [100 190]; % time in ms
pl.time2plot = [250 350]; % time in ms
% pl.time2plot = [400 550]; % time in ms
pl.time2plot_i = dsearchn(EP.time', pl.time2plot');


% preallocate data
pl.dat2plot = nan([size(EP.electrodes,2), ...           % 1st dim channels
    size(EP.time,2), ...                                % 2nd dim time
    cellfun(@(x) size(x,1), pl.con_contrast(:,2))', ... % 3rd to n dim = levels defind in pl.con_contrast
    numel(pl.sub2plot)]);                               % n+1 dim participants

% extract data across contrasts
% define how to loop through contrasts
t.cont = cellfun(@(x) 1:size(x,1), pl.con_contrast(:,2), 'UniformOutput', false);
t.contidx = combvec(t.cont{:});
% loop across participants
for i_sub = 1:numel(pl.sub2plot)
    t.behavior = EP.behavior{pl.sub2plot(i_sub)};
    t.epdata = EP.filtdata{pl.sub2plot(i_sub)}.data;

    % loop through contrasts
    for i_cont = 1:size(t.contidx,2)
        % define evaluation syntax
                
        % preallocate logical idx for trials
        t.idx = false(size(t.contidx,1),numel(t.behavior ));
        for i_contstep = 1:size(t.contidx,1)
            % index step
            t.idx(i_contstep,:) = ...
                any( ... % if it is more than one level, combine
                cell2mat( ...
                cellfun(@(x) strcmp({t.behavior.(pl.con_contrast{i_contstep,1})},x), ...
                pl.con_contrast{i_contstep,2}{t.contidx(i_contstep,i_cont)}, 'UniformOutput', false)' ...
                ) ...
                ,1);
        end
        % combine all contrast conditions
        t.idx_all = all(t.idx,1);
        
        % average across trials
        t.text = sprintf("%1.0f,",t.contidx(:,i_cont));
        t.evaltext = sprintf("pl.dat2plot(:,:,%si_sub) = mean(t.epdata(:,:,t.idx_all),3);",t.text);
        eval(t.evaltext)
    end
end

% reshape pl.data
pl.idx = pl.time2plot_i(1):pl.time2plot_i(2);
pl.dat2plot_rs = squeeze(mean(pl.dat2plot(:,pl.idx,:,:,:,:),2));

pl.data_m = mean(pl.dat2plot_rs,3);

pl.con_label = [pl.con_contrast{end,2}{:}];

pl.clims = [-1 1]*max(abs(pl.data_m),[],"all");

figure;
set(gcf,'Position',[100 100 1100 200],'PaperPositionMode','auto')

for i_con = 1:size(pl.data_m,2)
    subplot(1,size(pl.data_m,2)+1,i_con)
    
    topoplot(pl.data_m(:,i_con), EP.electrodes(1:64), ...
            'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clims,'conv','on','colormap',flipud(cbrewer2('RdBu')),...
            'whitebk','on'); % 'colormap',fake_parula; 'colormap',flipud(cbrewer2('RdBu'))
    
    title(sprintf("%s\n[%1.0f %1.0f]ms",pl.con_label{i_con},pl.time2plot))
    colorbar
end
% mean across conditions
subplot(1,size(pl.data_m,2)+1,i_con+1)

topoplot(mean(pl.data_m,2), EP.electrodes(1:64), ...
    'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clims,'conv','on','colormap',flipud(cbrewer2('RdBu')),...
    'whitebk','on'); % 'colormap',fake_parula; 'colormap',flipud(cbrewer2('RdBu'))

title(sprintf("mean con\n[%1.0f %1.0f]ms",pl.time2plot))
colorbar

%% export single trial erp data for R analysis
pl.sub2plot = 1:numel(F.Subs2use);
pl.erpparams = { ...
    'PD130', [100 190], {'P8';'PO8';'P10';'P7';'PO7';'P9'}; ...
    'N2', [250 330], {'P7';'PO7';'P9';'O1';'I1';'Oz'; 'Iz';'O2';'I2';'P8';'PO8';'P10'}; ...
    'P300', [400 550], {'P3';'P1';'Pz';'P4';'P2';'POz';'PO3';'PO4'} ...
    };

% append behavioral data
out.data = horzcat(EP.behavior{pl.sub2plot});

% add participant number
t.subs = [];
for i_sub = 1:numel(pl.sub2plot)
    t.subs = [t.subs; ...
        repmat( ...
        num2strcell(F.Subs2use(pl.sub2plot(i_sub))), ...
        numel(EP.behavior{pl.sub2plot(i_sub)}), ...
        1)];
end
[out.data.participants] = deal(t.subs{:});


% extract data
for i_erp = 1:size(pl.erpparams,1)
    t.erpdata = [];
    for i_sub = 1:numel(pl.sub2plot)
        % index electrodes
        pl.elec2plot_i = ...
            any(cell2mat(cellfun(@(x) strcmp({EP.electrodes.labels},x), pl.erpparams{i_erp,3}, 'UniformOutput',false)),1);
        pl.time2plot_i = dsearchn(EP.time', pl.erpparams{i_erp,2}');
        % extract erpdata
        t.erpdata = [t.erpdata; ...
            squeeze(mean(EP.filtdata{pl.sub2plot(i_sub)}.data(pl.elec2plot_i,pl.time2plot_i(1):pl.time2plot_i(2),:),[1,2]))];
    end
    % append values
    t.erpdata_cell = num2cell(t.erpdata);
    [out.data.(pl.erpparams{i_erp,1})] = deal(t.erpdata_cell{:});
    % add electrode values
    t.string = sprintf('%s_electrodes',pl.erpparams{i_erp,1});
    [out.data.(t.string)] = deal(vararg2str(pl.erpparams{i_erp,3}));
    % add time values
    t.string = sprintf('%s_timewindow',pl.erpparams{i_erp,1});
    [out.data.(t.string)] = deal(sprintf('[%1.0f %1.0f]ms',pl.erpparams{i_erp,2}));
end

t.path = 'C:\Users\psy05cvd\Dropbox\work\R-statistics\experiments\ssvep_fshiftprobabil\data';
t.datestr = datestr(now,'mm-dd-yyyy_HH-MM');
% t.filename = 'FFT_SSVEP_Amp_data_withoutBehav_sepRDK_peripheral_2elecs';
t.filename = 'ERP_data';
out.data_t = struct2table(out.data);
% remove some data
out.data_t = removevars(out.data_t,{'color_attended','eventRDK_color','event_color','button_presses_fr','button_presses_t'});

% write to textfile
% writetable(out.data_t,fullfile(t.path,sprintf('%s.csv',t.filename)),'Delimiter',';')


%% plot some erp images across time [cue validity]
% loop through conditions defined in contrasts
pl.con_contrast = {... % contrasts by 1st dim; averaged across second dim
    'trial_timing_type', {{'regular'}}; ...
    'event_response_type', {{'hit'}}; ...
%     'evnt_type_label', {{'chroma+'}}; ...
%     'event_response_type', {{'hit','FA','error','miss'}}; ...
    'cue_validity_label', {{'valid'};{'neutral'};{'invalid'}}};
pl.sub2plot = 1:numel(F.Subs2use);

pl.p_time = [25 25 -100 800]; % width step min max
pl.posScale = 1.1;

pl.pl.p_time_i = dsearchn(EP.time', pl.p_time(3:4)');

pl.concols = num2cell([240 63 240; 100 63 100]'./255,1);
pl.con_label = {'valid';'neutral';'invalid'};

pl.plegend = {'1';'0.5';'0.25';'0.1';'0.05';'0.01';'0.001';'0.0001';'0.00001'};
pl.pcorrect = [0 abs(log10([0.5 0.25 0.1 0.05 0.01 0.001 0.0001 0.00001]))];


% preallocate data
pl.dat2plot = nan([size(EP.electrodes,2), ...           % 1st dim channels
    size(EP.time,2), ...                                % 2nd dim time
    cellfun(@(x) size(x,1), pl.con_contrast(:,2))', ... % 3rd to n dim = levels defind in pl.con_contrast
    numel(pl.sub2plot)]);                               % n+1 dim participants

% extract data across contrasts
% define how to loop through contrasts
t.cont = cellfun(@(x) 1:size(x,1), pl.con_contrast(:,2), 'UniformOutput', false);
t.contidx = combvec(t.cont{:});
% loop across participants
for i_sub = 1:numel(pl.sub2plot)
    t.behavior = EP.behavior{pl.sub2plot(i_sub)};
    t.epdata = EP.filtdata{pl.sub2plot(i_sub)}.data;

    % loop through contrasts
    for i_cont = 1:size(t.contidx,2)
        % define evaluation syntax
                
        % preallocate logical idx for trials
        t.idx = false(size(t.contidx,1),numel(t.behavior ));
        for i_contstep = 1:size(t.contidx,1)
            % index step
            t.idx(i_contstep,:) = ...
                any( ... % if it is more than one level, combine
                cell2mat( ...
                cellfun(@(x) strcmp({t.behavior.(pl.con_contrast{i_contstep,1})},x), ...
                pl.con_contrast{i_contstep,2}{t.contidx(i_contstep,i_cont)}, 'UniformOutput', false)' ...
                ) ...
                ,1);
        end
        % combine all contrast conditions
        t.idx_all = all(t.idx,1);
        
        % average across trials
        t.text = sprintf("%1.0f,",t.contidx(:,i_cont));
        t.evaltext = sprintf("pl.dat2plot(:,:,%si_sub) = mean(t.epdata(:,:,t.idx_all),3);",t.text);
        eval(t.evaltext)
    end
end

% squeeze data
pl.dat2plot = squeeze(pl.dat2plot);

% calcualte contrast differences
t.contridx = nchoosek(1:size(pl.dat2plot,3),2);
pl.typeidx = [ones(1,size(pl.dat2plot,3)) ones(1,size(t.contridx,1))+1 ones(1,size(t.contridx,1))+2];
pl.dat2plot2 = pl.dat2plot;
for i_contr = 1:size(t.contridx,1)
    pl.dat2plot2(:,:,end+1,:) = pl.dat2plot(:,:,t.contridx(i_contr,1),:) - pl.dat2plot(:,:,t.contridx(i_contr,2),:);
    pl.con_label{end+1} = sprintf('%s-%s',pl.con_label{t.contridx(i_contr,1)},pl.con_label{t.contridx(i_contr,2)});
end
pl.con_label(pl.typeidx==3) = pl.con_label(pl.typeidx==2); % append for ttests

% plot amplitudes of all values
t.time=[];
t.timedot=[];
for i_st = 1:floor((pl.p_time(4)-pl.p_time(1)-pl.p_time(3))/pl.p_time(2))+1
    t.time(i_st,:)=pl.p_time(3)+pl.p_time(2)*(i_st-1)+[0 pl.p_time(1)];
    t.timedot(i_st,1:2)=dsearchn(EP.time',t.time(end,1:2)');
end
% t.t = get(0,'MonitorPositions');
% t.row = round(sqrt((size(t.time,1)+2)/(1/(t.t(1,4)/t.t(1,3)))));
% t.col = ceil(sqrt((size(t.time,1)+2)/(t.t(1,4)/t.t(1,3))));
t.t = [1920 1080];
t.col = ceil(sqrt((size(t.time,1)+2)/(1/(t.t(1)/t.t(2)))));
t.row = ceil(sqrt((size(t.time,1)+2)/(t.t(1)/t.t(2))));
% create plotdata

plotdata=[];
for i_pl = 1:size(t.time,1)
    plotdata(:,:,i_pl)=squeeze(mean(pl.dat2plot2(:,t.timedot(i_pl,1):t.timedot(i_pl,2),:,:),[2,4]));
end

% add ttest data
t.idx = size(pl.dat2plot2,3); % index where data should be appended
t.idx2 = size(pl.dat2plot,3); % index where contrasts start
for i_pl = 1:size(t.time,1)
    for i_tt = 1:size(t.contridx,1)
        t.data = squeeze(mean(pl.dat2plot2(:,t.timedot(i_pl,1):t.timedot(i_pl,2),t.idx2+i_tt,:),[2]));
        [tt.h tt.p tt.ci tt.stats] = ttest(t.data,0,'Dim',2);
        plotdata(:,t.idx+i_tt,i_pl) = tt.p;
    end
end


h.fig = []; h.sp = [];

for i_fig = 1:size(plotdata,2)
    h.fig(i_fig)=figure;
    set(gcf,'Position',[100 100 1200 700],'PaperPositionMode','auto')
    if pl.typeidx(i_fig) == 1
        t.lims = [-1 1]*max(abs(plotdata(:,pl.typeidx==1,:)),[],'all');
    elseif pl.typeidx(i_fig) == 2
        t.lims = [-1 1]*max(abs(plotdata(:,pl.typeidx==2,:)),[],'all');
        % t.lims = [-1 1]*max(abs(plotdata(:,i_fig,:)),[],'all');
    elseif pl.typeidx(i_fig) == 3
        % t.lims = [0 max(abs(log10(plotdata(:,pl.typeidx==3,:))),[],'all')];
        t.lims = [0 max(abs(log10(plotdata(:,i_fig,:))),[],'all')];
        
        % define colormaps
        t.pcriterion = abs(log10(0.05));
        if t.lims(2)<t.pcriterion
            t.colormap = repmat(linspace(1,0.3,1000)',1,3);
        else
            t.border = ceil((t.pcriterion / t.lims(2))*1000);
            % temp.colormap = [repmat([0.5 0.5 0.5],temp.border,1); [linspace(0,1,100-temp.border)' repmat(0,100-temp.border,1) linspace(1,0,100-temp.border)']];
            t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
        end
            
    end
    for i_spl = 1:size(plotdata,3)
        h.sp(i_spl)=subplot(t.row,t.col,i_spl);
        if pl.typeidx(i_fig) == 1 || pl.typeidx(i_fig) == 2
%             topoplot( plotdata(:,i_fig,i_spl), EP.electrodes(1:64), ...
%                 'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',t.lims, ...
%                 'electrodes','off','colormap',fake_parula,'whitebk','on');
            topoplot( plotdata(:,i_fig,i_spl), EP.electrodes(1:64), ...
                'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',t.lims, ...
                'electrodes','off','colormap',flipud(cbrewer2('RdBu')),'whitebk','on');
        elseif pl.typeidx(i_fig) == 3
            topoplot( abs(log10(plotdata(:,i_fig,i_spl))), EP.electrodes(1:64), ...
                'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',t.lims, ...
                'electrodes','off','colormap',t.colormap,'whitebk','on');
        end
        title(sprintf('[%1.0f %1.0f]',t.time(i_spl,1),t.time(i_spl,2)),'FontSize',6)
        t.pos = get(h.sp(i_spl),'Position');
        set(h.sp(i_spl),'Position',[t.pos(1:2)-(t.pos(3:4).*((pl.posScale-1)/2)) t.pos(3:4).*pl.posScale])
    end
    h.sp(i_spl+1)=subplot(t.row,t.col,i_spl+2);
    topoplot( [], EP.electrodes(1:64),  ...
        'style','blank','whitebk','on');
    title(sprintf('%s',pl.con_label{i_fig}),'FontSize',8)
    t.pos = get(h.sp(i_spl+1),'Position');
    set(h.sp(i_spl+1),'Position',[t.pos(1:2)-(t.pos(3:4).*((pl.posScale-1)/2)) t.pos(3:4).*pl.posScale])
    
    t.pos2 = get(h.sp(i_spl+1),'Position');
    t.pos3 = get(h.sp(i_spl+1),'OuterPosition');
    h.a1 = axes('position',[t.pos3(1) t.pos2(2) t.pos3(3) t.pos2(4)],'Visible','off');
    if pl.typeidx(i_fig) == 3
        colormap(gca, t.colormap)
        clim(t.lims);
        h.c3 = colorbar();
        t.pos4 = get(h.c3,'Position');
        set(h.c3,'Position',[t.pos4(1)+0.065 t.pos4(2)+(t.pos4(4)*(1/6)) t.pos4(3) t.pos4(4)*(2/3)])
        t.yticks = get(h.c3,'YTick');
        set(h.c3,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.lims(end),1,'last')), ...
            'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.lims(end),1,'last')))
    else
        colormap(gca,flipud(cbrewer2('RdBu')))
        clim(t.lims);
        h.c3 = colorbar();
        t.pos4 = get(h.c3,'Position');
        set(h.c3,'Position',[t.pos4(1)+0.065 t.pos4(2)+(t.pos4(4)*(1/6)) t.pos4(3) t.pos4(4)*(2/3)])
    end
    axcopy(h.fig(i_fig))
end



%% plot some erp images across time [hit or no hit]
% loop through conditions defined in contrasts
pl.con_contrast = {... % contrasts by 1st dim; averaged across second dim
    'trial_timing_type', {{'regular'}}; ...
    'evnt_type_label', {{'chroma+'}}; ...
    'event_response_type', {{'hit'};{'FA','error','miss'}}};
pl.sub2plot = 1:numel(F.Subs2use);

pl.p_time = [25 25 -100 800]; % width step min max
pl.posScale = 1.1;

pl.pl.p_time_i = dsearchn(EP.time', pl.p_time(3:4)');

pl.concols = num2cell([240 63 240; 100 63 100]'./255,1);
pl.con_label = {'hit';'error+FA+miss'};

pl.plegend = {'1';'0.5';'0.25';'0.1';'0.05';'0.01';'0.001';'0.0001';'0.00001'};
pl.pcorrect = [0 abs(log10([0.5 0.25 0.1 0.05 0.01 0.001 0.0001 0.00001]))];


% preallocate data
pl.dat2plot = nan([size(EP.electrodes,2), ...           % 1st dim channels
    size(EP.time,2), ...                                % 2nd dim time
    cellfun(@(x) size(x,1), pl.con_contrast(:,2))', ... % 3rd to n dim = levels defind in pl.con_contrast
    numel(pl.sub2plot)]);                               % n+1 dim participants

% extract data across contrasts
% define how to loop through contrasts
t.cont = cellfun(@(x) 1:size(x,1), pl.con_contrast(:,2), 'UniformOutput', false);
t.contidx = combvec(t.cont{:});
% loop across participants
for i_sub = 1:numel(pl.sub2plot)
    t.behavior = EP.behavior{pl.sub2plot(i_sub)};
    t.epdata = EP.filtdata{pl.sub2plot(i_sub)}.data;

    % loop through contrasts
    for i_cont = 1:size(t.contidx,2)
        % define evaluation syntax
                
        % preallocate logical idx for trials
        t.idx = false(size(t.contidx,1),numel(t.behavior ));
        for i_contstep = 1:size(t.contidx,1)
            % index step
            t.idx(i_contstep,:) = ...
                any( ... % if it is more than one level, combine
                cell2mat( ...
                cellfun(@(x) strcmp({t.behavior.(pl.con_contrast{i_contstep,1})},x), ...
                pl.con_contrast{i_contstep,2}{t.contidx(i_contstep,i_cont)}, 'UniformOutput', false)' ...
                ) ...
                ,1);
        end
        % combine all contrast conditions
        t.idx_all = all(t.idx,1);
        
        % average across trials
        t.text = sprintf("%1.0f,",t.contidx(:,i_cont));
        t.evaltext = sprintf("pl.dat2plot(:,:,%si_sub) = mean(t.epdata(:,:,t.idx_all),3);",t.text);
        eval(t.evaltext)
    end
end

% squeeze data
pl.dat2plot = squeeze(pl.dat2plot);

% calcualte contrast differences
t.contridx = nchoosek(1:size(pl.dat2plot,3),2);
pl.typeidx = [ones(1,size(pl.dat2plot,3)) ones(1,size(t.contridx,1))+1 ones(1,size(t.contridx,1))+2];
pl.dat2plot2 = pl.dat2plot;
for i_contr = 1:size(t.contridx,1)
    pl.dat2plot2(:,:,end+1,:) = pl.dat2plot(:,:,t.contridx(i_contr,1),:) - pl.dat2plot(:,:,t.contridx(i_contr,2),:);
    pl.con_label{end+1} = sprintf('%s-%s',pl.con_label{t.contridx(i_contr,1)},pl.con_label{t.contridx(i_contr,2)});
end
pl.con_label(pl.typeidx==3) = pl.con_label(pl.typeidx==2); % append for ttests

% plot amplitudes of all values
t.time=[];
t.timedot=[];
for i_st = 1:floor((pl.p_time(4)-pl.p_time(1)-pl.p_time(3))/pl.p_time(2))+1
    t.time(i_st,:)=pl.p_time(3)+pl.p_time(2)*(i_st-1)+[0 pl.p_time(1)];
    t.timedot(i_st,1:2)=dsearchn(EP.time',t.time(end,1:2)');
end
% t.t = get(0,'MonitorPositions');
% t.row = round(sqrt((size(t.time,1)+2)/(1/(t.t(1,4)/t.t(1,3)))));
% t.col = ceil(sqrt((size(t.time,1)+2)/(t.t(1,4)/t.t(1,3))));
t.t = [1920 1080];
t.col = ceil(sqrt((size(t.time,1)+2)/(1/(t.t(1)/t.t(2)))));
t.row = ceil(sqrt((size(t.time,1)+2)/(t.t(1)/t.t(2))));
% create plotdata

plotdata=[];
for i_pl = 1:size(t.time,1)
    plotdata(:,:,i_pl)=squeeze(mean(pl.dat2plot2(:,t.timedot(i_pl,1):t.timedot(i_pl,2),:,:),[2,4]));
end

% add ttest data
t.idx = size(pl.dat2plot2,3); % index where data should be appended
t.idx2 = size(pl.dat2plot,3); % index where contrasts start
for i_pl = 1:size(t.time,1)
    for i_tt = 1:size(t.contridx,1)
        t.data = squeeze(mean(pl.dat2plot2(:,t.timedot(i_pl,1):t.timedot(i_pl,2),t.idx2+i_tt,:),[2]));
        [tt.h tt.p tt.ci tt.stats] = ttest(t.data,0,'Dim',2);
        plotdata(:,t.idx+i_tt,i_pl) = tt.p;
    end
end

h.fig = []; h.sp = [];

for i_fig = 1:size(plotdata,2)
    h.fig(i_fig)=figure;
    set(gcf,'Position',[100 100 1200 700],'PaperPositionMode','auto')
    if pl.typeidx(i_fig) == 1
        t.lims = [-1 1]*max(abs(plotdata(:,pl.typeidx==1,:)),[],'all');
    elseif pl.typeidx(i_fig) == 2
        t.lims = [-1 1]*max(abs(plotdata(:,pl.typeidx==2,:)),[],'all');
        % t.lims = [-1 1]*max(abs(plotdata(:,i_fig,:)),[],'all');
    elseif pl.typeidx(i_fig) == 3
        % t.lims = [0 max(abs(log10(plotdata(:,pl.typeidx==3,:))),[],'all')];
        t.lims = [0 max(abs(log10(plotdata(:,i_fig,:))),[],'all')];

        % define colormaps
        t.pcriterion = abs(log10(0.05));
        if t.lims(2)<t.pcriterion
            t.colormap = repmat(linspace(1,0.3,1000)',1,3);
        else
            t.border = ceil((t.pcriterion / t.lims(2))*1000);
            % temp.colormap = [repmat([0.5 0.5 0.5],temp.border,1); [linspace(0,1,100-temp.border)' repmat(0,100-temp.border,1) linspace(1,0,100-temp.border)']];
            t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
        end
    end
    for i_spl = 1:size(plotdata,3)
        h.sp(i_spl)=subplot(t.row,t.col,i_spl);
        if pl.typeidx(i_fig) == 1 || pl.typeidx(i_fig) == 2
            topoplot(plotdata(:,i_fig,i_spl), EP.electrodes(1:64), ...
                'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',t.lims, ...
                'electrodes','off','colormap',flipud(cbrewer2('RdBu')),'whitebk','on');
        elseif pl.typeidx(i_fig) == 3
            topoplot( abs(log10(plotdata(:,i_fig,i_spl))), EP.electrodes(1:64), ...
                'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',t.lims, ...
                'electrodes','off','colormap',t.colormap,'whitebk','on');
        end
       
        title(sprintf('[%1.0f %1.0f]',t.time(i_spl,1),t.time(i_spl,2)),'FontSize',6)
        t.pos = get(h.sp(i_spl),'Position');
        set(h.sp(i_spl),'Position',[t.pos(1:2)-(t.pos(3:4).*((pl.posScale-1)/2)) t.pos(3:4).*pl.posScale])
    end
    h.sp(i_spl+1)=subplot(t.row,t.col,i_spl+2);
    topoplot( [], EP.electrodes(1:64),  ...
        'style','blank','whitebk','on');
    title(sprintf('%s',pl.con_label{i_fig}),'FontSize',8)
    t.pos = get(h.sp(i_spl+1),'Position');
    set(h.sp(i_spl+1),'Position',[t.pos(1:2)-(t.pos(3:4).*((pl.posScale-1)/2)) t.pos(3:4).*pl.posScale])
    
    t.pos2 = get(h.sp(i_spl+1),'Position');
    t.pos3 = get(h.sp(i_spl+1),'OuterPosition');
    h.a1 = axes('position',[t.pos3(1) t.pos2(2) t.pos3(3) t.pos2(4)],'Visible','off');
    if pl.typeidx(i_fig) == 3
        colormap(gca, t.colormap)
        clim(t.lims);
        h.c3 = colorbar();
        t.pos4 = get(h.c3,'Position');
        set(h.c3,'Position',[t.pos4(1)+0.065 t.pos4(2)+(t.pos4(4)*(1/6)) t.pos4(3) t.pos4(4)*(2/3)])
        t.yticks = get(h.c3,'YTick');
        set(h.c3,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.lims(end),1,'last')), ...
            'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.lims(end),1,'last')))
    else
        colormap(gca,flipud(cbrewer2('RdBu')))
        clim(t.lims);
        h.c3 = colorbar();
        t.pos4 = get(h.c3,'Position');
        set(h.c3,'Position',[t.pos4(1)+0.065 t.pos4(2)+(t.pos4(4)*(1/6)) t.pos4(3) t.pos4(4)*(2/3)])
    end
    axcopy(h.fig(i_fig))
end


%% plot some erp images across time [chroma]
% loop through conditions defined in contrasts
pl.con_contrast = {... % contrasts by 1st dim; averaged across second dim
    'trial_timing_type', {{'regular'}}; ...
    'event_response_type', {{'hit'}}; ...
%     'event_response_type', {{'hit','FA','error','miss'}}; ...
    'evnt_type_label', {{'chroma+'};{'chroma-'}}};
pl.sub2plot = 1:numel(F.Subs2use);

pl.p_time = [25 25 -100 800]; % width step min max
pl.posScale = 1.1;

pl.pl.p_time_i = dsearchn(EP.time', pl.p_time(3:4)');

pl.concols = num2cell([240 63 240; 100 63 100]'./255,1);
pl.con_label = {'chroma+';'chroma-'};

pl.plegend = {'1';'0.5';'0.25';'0.1';'0.05';'0.01';'0.001';'0.0001';'0.00001'};
pl.pcorrect = [0 abs(log10([0.5 0.25 0.1 0.05 0.01 0.001 0.0001 0.00001]))];


% preallocate data
pl.dat2plot = nan([size(EP.electrodes,2), ...           % 1st dim channels
    size(EP.time,2), ...                                % 2nd dim time
    cellfun(@(x) size(x,1), pl.con_contrast(:,2))', ... % 3rd to n dim = levels defind in pl.con_contrast
    numel(pl.sub2plot)]);                               % n+1 dim participants

% extract data across contrasts
% define how to loop through contrasts
t.cont = cellfun(@(x) 1:size(x,1), pl.con_contrast(:,2), 'UniformOutput', false);
t.contidx = combvec(t.cont{:});
% loop across participants
for i_sub = 1:numel(pl.sub2plot)
    t.behavior = EP.behavior{pl.sub2plot(i_sub)};
    t.epdata = EP.filtdata{pl.sub2plot(i_sub)}.data;

    % loop through contrasts
    for i_cont = 1:size(t.contidx,2)
        % define evaluation syntax
                
        % preallocate logical idx for trials
        t.idx = false(size(t.contidx,1),numel(t.behavior ));
        for i_contstep = 1:size(t.contidx,1)
            % index step
            t.idx(i_contstep,:) = ...
                any( ... % if it is more than one level, combine
                cell2mat( ...
                cellfun(@(x) strcmp({t.behavior.(pl.con_contrast{i_contstep,1})},x), ...
                pl.con_contrast{i_contstep,2}{t.contidx(i_contstep,i_cont)}, 'UniformOutput', false)' ...
                ) ...
                ,1);
        end
        % combine all contrast conditions
        t.idx_all = all(t.idx,1);
        
        % average across trials
        t.text = sprintf("%1.0f,",t.contidx(:,i_cont));
        t.evaltext = sprintf("pl.dat2plot(:,:,%si_sub) = mean(t.epdata(:,:,t.idx_all),3);",t.text);
        eval(t.evaltext)
    end
end

% squeeze data
pl.dat2plot = squeeze(pl.dat2plot);

% calcualte contrast differences
t.contridx = nchoosek(1:size(pl.dat2plot,3),2);
pl.typeidx = [ones(1,size(pl.dat2plot,3)) ones(1,size(t.contridx,1))+1 ones(1,size(t.contridx,1))+2];
pl.dat2plot2 = pl.dat2plot;
for i_contr = 1:size(t.contridx,1)
    pl.dat2plot2(:,:,end+1,:) = pl.dat2plot(:,:,t.contridx(i_contr,1),:) - pl.dat2plot(:,:,t.contridx(i_contr,2),:);
    pl.con_label{end+1} = sprintf('%s-%s',pl.con_label{t.contridx(i_contr,1)},pl.con_label{t.contridx(i_contr,2)});
end
pl.con_label(pl.typeidx==3) = pl.con_label(pl.typeidx==2); % append for ttests

% plot amplitudes of all values
t.time=[];
t.timedot=[];
for i_st = 1:floor((pl.p_time(4)-pl.p_time(1)-pl.p_time(3))/pl.p_time(2))+1
    t.time(i_st,:)=pl.p_time(3)+pl.p_time(2)*(i_st-1)+[0 pl.p_time(1)];
    t.timedot(i_st,1:2)=dsearchn(EP.time',t.time(end,1:2)');
end
% t.t = get(0,'MonitorPositions');
% t.row = round(sqrt((size(t.time,1)+2)/(1/(t.t(1,4)/t.t(1,3)))));
% t.col = ceil(sqrt((size(t.time,1)+2)/(t.t(1,4)/t.t(1,3))));
t.t = [1920 1080];
t.col = ceil(sqrt((size(t.time,1)+2)/(1/(t.t(1)/t.t(2)))));
t.row = ceil(sqrt((size(t.time,1)+2)/(t.t(1)/t.t(2))));
% create plotdata

plotdata=[];
for i_pl = 1:size(t.time,1)
    plotdata(:,:,i_pl)=squeeze(mean(pl.dat2plot2(:,t.timedot(i_pl,1):t.timedot(i_pl,2),:,:),[2,4]));
end

% add ttest data
t.idx = size(pl.dat2plot2,3); % index where data should be appended
t.idx2 = size(pl.dat2plot,3); % index where contrasts start
for i_pl = 1:size(t.time,1)
    for i_tt = 1:size(t.contridx,1)
        t.data = squeeze(mean(pl.dat2plot2(:,t.timedot(i_pl,1):t.timedot(i_pl,2),t.idx2+i_tt,:),[2]));
        [tt.h tt.p tt.ci tt.stats] = ttest(t.data,0,'Dim',2);
        plotdata(:,t.idx+i_tt,i_pl) = tt.p;
    end
end


h.fig = []; h.sp = [];

for i_fig = 1:size(plotdata,2)
    h.fig(i_fig)=figure;
    set(gcf,'Position',[100 100 1200 700],'PaperPositionMode','auto')
    if pl.typeidx(i_fig) == 1
        t.lims = [-1 1]*max(abs(plotdata(:,pl.typeidx==1,:)),[],'all');
    elseif pl.typeidx(i_fig) == 2
        t.lims = [-1 1]*max(abs(plotdata(:,pl.typeidx==2,:)),[],'all');
        % t.lims = [-1 1]*max(abs(plotdata(:,i_fig,:)),[],'all');
    elseif pl.typeidx(i_fig) == 3
        % t.lims = [0 max(abs(log10(plotdata(:,pl.typeidx==3,:))),[],'all')];
        t.lims = [0 max(abs(log10(plotdata(:,i_fig,:))),[],'all')];

        % define colormaps
        t.pcriterion = abs(log10(0.05));
        if t.lims(2)<t.pcriterion
            t.colormap = repmat(linspace(1,0.3,1000)',1,3);
        else
            t.border = ceil((t.pcriterion / t.lims(2))*1000);
            % temp.colormap = [repmat([0.5 0.5 0.5],temp.border,1); [linspace(0,1,100-temp.border)' repmat(0,100-temp.border,1) linspace(1,0,100-temp.border)']];
            t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
        end
    end
    for i_spl = 1:size(plotdata,3)
        h.sp(i_spl)=subplot(t.row,t.col,i_spl);
        if pl.typeidx(i_fig) == 1 || pl.typeidx(i_fig) == 2
            topoplot(plotdata(:,i_fig,i_spl), EP.electrodes(1:64), ...
                'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',t.lims, ...
                'electrodes','off','colormap',flipud(cbrewer2('RdBu')),'whitebk','on');
        elseif pl.typeidx(i_fig) == 3
            topoplot( abs(log10(plotdata(:,i_fig,i_spl))), EP.electrodes(1:64), ...
                'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',t.lims, ...
                'electrodes','off','colormap',t.colormap,'whitebk','on');
        end
       
        title(sprintf('[%1.0f %1.0f]',t.time(i_spl,1),t.time(i_spl,2)),'FontSize',6)
        t.pos = get(h.sp(i_spl),'Position');
        set(h.sp(i_spl),'Position',[t.pos(1:2)-(t.pos(3:4).*((pl.posScale-1)/2)) t.pos(3:4).*pl.posScale])
    end
    h.sp(i_spl+1)=subplot(t.row,t.col,i_spl+2);
    topoplot( [], EP.electrodes(1:64),  ...
        'style','blank','whitebk','on');
    title(sprintf('%s',pl.con_label{i_fig}),'FontSize',8)
    t.pos = get(h.sp(i_spl+1),'Position');
    set(h.sp(i_spl+1),'Position',[t.pos(1:2)-(t.pos(3:4).*((pl.posScale-1)/2)) t.pos(3:4).*pl.posScale])
    
    t.pos2 = get(h.sp(i_spl+1),'Position');
    t.pos3 = get(h.sp(i_spl+1),'OuterPosition');
    h.a1 = axes('position',[t.pos3(1) t.pos2(2) t.pos3(3) t.pos2(4)],'Visible','off');
    if pl.typeidx(i_fig) == 3
        colormap(gca, t.colormap)
        clim(t.lims);
        h.c3 = colorbar();
        t.pos4 = get(h.c3,'Position');
        set(h.c3,'Position',[t.pos4(1)+0.065 t.pos4(2)+(t.pos4(4)*(1/6)) t.pos4(3) t.pos4(4)*(2/3)])
        t.yticks = get(h.c3,'YTick');
        set(h.c3,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.lims(end),1,'last')), ...
            'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.lims(end),1,'last')))
    else
        colormap(gca,flipud(cbrewer2('RdBu')))
        clim(t.lims);
        h.c3 = colorbar();
        t.pos4 = get(h.c3,'Position');
        set(h.c3,'Position',[t.pos4(1)+0.065 t.pos4(2)+(t.pos4(4)*(1/6)) t.pos4(3) t.pos4(4)*(2/3)])
    end
    axcopy(h.fig(i_fig))
end












