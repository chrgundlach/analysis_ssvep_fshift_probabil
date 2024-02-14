%% plot TFA images
clearvars
F.PathInEEG             = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2023_FShift_Probabil\eeg\erp'; 


F.Subs                  = arrayfun(@(x) sprintf('%02.0f',x),1:40,'UniformOutput',false)';
F.Subs2use              = [1 3 4 5 6 7 9 10 11 13 15 18 20 21 22 23 24 25]; 
                        % 2 and 8 are excluded as the didn't do the task properly, sub 11 has potentially low number of trials
                        % for subject 12, 14: eeg and behavior data don't match



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


%% plot ERPs exploratively for specified electrodes

% loop through conditions defined in contrasts
pl.con_contrast = {... % contrasts by 1st dim; averaged across second dim
    'trial_timing_type', {{'regular'}}; ...
    'cue_validity_label', {{'valid'};{'neutral'};{'invalid'}}};
pl.sub2plot = 1:numel(F.Subs2use);
pl.elec2plot = {'Pz';'POz';'PO3';'PO4'};
pl.elec2plot_i = ...
    any(cell2mat(cellfun(@(x) strcmp({EP.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1);

pl.time2plot = [-100 500]; % time in ms
pl.time2plot = [-100 700]; % time in ms
pl.time2plot_i = dsearchn(EP.time', pl.time2plot');

pl.concols = num2cell([255 133 4; 25 138 131; 255 194 64]'./255,1);


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

pl.con_label = [pl.con_contrast{2,2}{:}];

figure;
set(gcf,'Position',[100 100 500 200],'PaperPositionMode','auto')
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


h.ax2 = axes('Position', [0.785, .69, .19, .19],'Visible','off');
topoplot(find(pl.elec2plot_i{1}),ERP.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i{1}),'o','r',3,1});

h.ax3 = axes('Position', [0.65, .69, .19, .19],'Visible','off');
topoplot(find(pl.elec2plot_i),EP.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',3,1});
set(gcf, 'Color', [1 1 1]);
