%% script to analyze behavioral data
clearvars

p.path =                'D:\work\data\SSVEP_FShiftAlpha\logfiles\';
p.path =                '\\psyall2.misc.intern.uni-leipzig.de\experimental_data\2023_FShift_Probabil\behavior\raw\';
p.path =                'N:\AllgPsy\experimental_data\2023_FShift_Probabil\behavior\raw\';
p.pathout =             'C:\Users\psy05cvd\Dropbox\work\R-statistics\experiments\ssvep_fshiftprobabil\data\';
p.subs=                 arrayfun(@(x) sprintf('%02.0f',x),1:40,'UniformOutput',false)';
p.subs2use=             [1 3:6 7 8 9 10 11 12];% 
p.subs2use=             [1 3:6 7 9 10 11 12];%  % participant 2,8 measurement cancelled due to bad behavior
% p.subs2use =            [1 3 4 5];
% pl.sub2plot =           [1:14 16:35]; % sub 15 ith no SSVEP default = 1:numel(F.Subs2use)
p.responsewin =         [0.2 1.2]; % according to p.targ_respwin from run_FShiftAlpha





%% actual calculation
for i_sub = 1:numel(p.subs2use)
    % load data
    temp.files = dir(sprintf('%sVP%s_timing*.mat',p.path,p.subs{p.subs2use(i_sub)}));
    data_in.resp.experiment = repmat({[nan]},1,14);
    data_in.button_presses.experiment = repmat({[nan]},1,14);
    for i_fi = 1:numel(temp.files)
        temp.data_in{i_fi}=load(sprintf('%s%s',p.path,temp.files(i_fi).name));
        % extract relevant data
        try data_in.conmat.experiment = temp.data_in{i_fi}.conmat.experiment;
        end
        if any(strcmp(fieldnames(temp.data_in{i_fi}.resp),'experiment'))
            temp.index1 = find(~cellfun(@isempty,temp.data_in{i_fi}.resp.experiment));
            temp.index2 = cell2mat(cellfun(@(x) ~isempty(cell2mat({x(:).trialnumber})), temp.data_in{i_fi}.resp.experiment(temp.index1),'UniformOutput',false));
            data_in.resp.experiment(temp.index1(temp.index2))=temp.data_in{i_fi}.resp.experiment(temp.index1(temp.index2));
            data_in.button_presses.experiment(temp.index1(temp.index2))=temp.data_in{i_fi}.button_presses.experiment(temp.index1(temp.index2));
        end
    end
    
    
    %% loop for blocks
    for i_bl = 1:sum(cellfun(@isstruct,data_in.resp.experiment))
        t.fieldnames = fieldnames(data_in.resp.experiment{i_bl});
        t.fieldname_idx = [ones(1,24) 0 0 ones(1,2)]==1;
        data2append = data_in.resp.experiment{i_bl};
        %temp.cell = num2cell(repmat(str2num(p.subs{p.subs2use(i_sub)}),1,size(data_in.resp.experiment{i_bl},2)));
        temp.cell = repmat(p.subs(p.subs2use(i_sub)),1,size(data_in.resp.experiment{i_bl},2));
        [data2append.subject] = temp.cell{:};
        
        data2append = rmfield(data2append,t.fieldnames(~t.fieldname_idx));

        % add aditional information and deal with some errors
        % attended color label
        t.mat =  arrayfun(@(x) temp.data_in{1, 1}.RDK.RDK(x).col_label   , ...
            [data2append.cue], 'UniformOutput', false);
        [data2append.color_attended_label] = t.mat{:};
        t.mat = repmat(strcat(temp.data_in{1, 1}.RDK.RDK(1).col_label, {' + '}, temp.data_in{1, 1}.RDK.RDK(2).col_label), ...
            sum([data2append.cue]==3),1);
        [data2append([data2append.cue]==3).color_attended_label] = t.mat{:};
        % target color label
        t.mat =  arrayfun(@(x) temp.data_in{1, 1}.RDK.RDK(x).col_label   , ...
            [data2append.eventpos], 'UniformOutput', false);
        [data2append.eventRDK_color_label] = t.mat{:};
        % target type label
        t.mat1 = {'chroma_decrease';'chroma_increase'};
        t.mat2 = arrayfun(@(x) t.mat1{x}   , ...
            [data2append.event_type], 'UniformOutput', false);
        [data2append.event_type_label] = t.mat2{:};
        % post cue target time
        t.mat = num2cell([data2append.event_onset_times] - [data2append.pre_cue_times]);
        [data2append.event_onset_times_post_cue] = t.mat{:};


        
        % append data to response
        if i_sub == 1 & i_bl == 1
            response = data2append;
        else
            response = [response, data2append];
        end
    end
end

% responses counted as FAs will be recoded as misses
t.idx = strcmp({response.event_response_type},'FA');
[response(t.idx).event_response_type] = deal('miss');

%% intemediate statistics for subject
sub_2use = p.subs(p.subs2use);
% sub_2use = p.subs([5]);
clear avg_resp
for i_sub = 1:numel(sub_2use)
    sub_2use_idx = strcmp({response.subject},sub_2use(i_sub));
    t.idx = [sub_2use_idx & strcmp({response.cue_validity_label},'valid') & [response.eventpos]==1;...
        sub_2use_idx & strcmp({response.cue_validity_label},'neutral') & [response.eventpos]==1;...
        sub_2use_idx & strcmp({response.cue_validity_label},'invalid') & [response.eventpos]==1;...
        sub_2use_idx & strcmp({response.cue_validity_label},'valid') & [response.eventpos]==2;...
        sub_2use_idx & strcmp({response.cue_validity_label},'neutral') & [response.eventpos]==2;...
        sub_2use_idx & strcmp({response.cue_validity_label},'invalid') & [response.eventpos]==2];
    for i_con = 1:6
        avg_resp(i_sub).RT(i_con) = nanmean([response(t.idx(i_con,:)).event_response_RT]);
        avg_resp(i_sub).hitrate(i_con) = sum(strcmp({response(t.idx(i_con,:)).event_response_type},'hit'))/sum(t.idx(i_con,:));
        avg_resp(i_sub).missrate(i_con) = sum(strcmp({response(t.idx(i_con,:)).event_response_type},'miss'))/sum(t.idx(i_con,:));
        avg_resp(i_sub).errorrate(i_con) =  sum(strcmp({response(t.idx(i_con,:)).event_response_type},'error'))/sum(t.idx(i_con,:));
    end
    
    t.idx = [sub_2use_idx & strcmp({response.cue_validity_label},'valid') ;...
        sub_2use_idx & strcmp({response.cue_validity_label},'neutral');...
        sub_2use_idx & strcmp({response.cue_validity_label},'invalid') ];
    for i_con = 1:3
        avg_resp(i_sub).RT_bothcols(i_con) = nanmean([response(t.idx(i_con,:)).event_response_RT]);
        avg_resp(i_sub).hitrate_bothcols(i_con) = sum(strcmp({response(t.idx(i_con,:)).event_response_type},'hit'))/sum(t.idx(i_con,:));
        avg_resp(i_sub).missrate_bothcols(i_con) = sum(strcmp({response(t.idx(i_con,:)).event_response_type},'miss'))/sum(t.idx(i_con,:));
        avg_resp(i_sub).errorrate_bothcols(i_con) = sum(strcmp({response(t.idx(i_con,:)).event_response_type},'error'))/sum(t.idx(i_con,:));
    end
end

%% prepare data for export
t.fieldnames = fieldnames(response);
t.fieldname_idx = [1, 2, 3, 4, 6, 14, 17, 21, 23, 25, 26, 27, 28, 29, 30, 31];
t.fieldnames_idx2 = true(size(t.fieldnames));
t.fieldnames_idx2(t.fieldname_idx)=false;

% select only a few fields
response2export = rmfield(response,t.fieldnames(t.fieldnames_idx2));

% order fields
t.fieldnames = fieldnames(response2export);
t.idx = [12, 1, 2, 4, 3, 5, 7, 8, 13, 14, 6, 9, 16, 15, 10, 11];
response2export = orderfields(response2export, t.idx);

% export
response2export_t = struct2table(response2export);
t.datestr = datestr(now,'mm-dd-yyyy_HH-MM');
t.filename = 'behavior';
% writetable(response2export_t,fullfile(p.pathout,sprintf('%s_%s.csv',t.filename,t.datestr)),'Delimiter',';')
% writetable(response2export_t,fullfile(p.pathout,sprintf('%s.csv',t.filename)),'Delimiter',';')