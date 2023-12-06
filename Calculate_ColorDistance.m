% script to tranfer RGB colors and calculate distance
clearvars

p.path =              'C:\Users\psy05cvd\Dropbox\work\documents\Experimente\SSVEP_FShiftBase\Behavioral_data\';
p.subs=                 arrayfun(@(x) sprintf('%02.0f',x),1:30,'UniformOutput',false)';
p.subs2use=             [1:30];% 
% p.subs2use =            [1 3 4 5];
% p.subs2use =            [5 6 7 8 9 10 12 13 18 19 20 21 22 23 24 25];
p.responsewin =         [0.2 1]; % according to p.targ_respwin from run_posnalpha


p.cols =                [0    1.0000         0; 0    0.4000    1.0000; 1.0000    0.4000         0];
p.colnames =            {'green'; 'cyan'; 'orange'};




col_RGB = [0    1.0000         0; 0    0.4000    1.0000; 1.0000    0.4000         0];

% propixx color values
col.xy_r = [0.665 0.321]; col.xy_g = [0.172 0.726]; col.xy_b = [0.163 0.039]; col.xyz_w = [0.3127 0.3290];
% col_RGB = col_RGB.*[0.5 0.5 0.5; 0.3 0.3 0.3; 0.4 0.4 0.4];
% [xyz] = rgb2xyz_custom([1 0 0],col.xy_r, col.xy_g, col.xy_b, col.xyz_w);
% [RGB] = xyz2rgb_custom([0.6650    0.3210    0.0140],col.xy_r, col.xy_g, col.xy_b, col.xyz_w);


%% extract data for each subject and calculate distances
for i_sub = 1:numel(p.subs2use)
    % load data
    temp.files = dir(sprintf('%sVP%s_timing*.mat',p.path,p.subs{p.subs2use(i_sub)}));
    data_in.resp.experiment = repmat({[nan]},1,8);
    data_in.button_presses.experiment = repmat({[nan]},1,8);
    for i_fi = 1:numel(temp.files)
        temp.data_in{i_fi}=load(sprintf('%s%s',p.path,temp.files(i_fi).name));
        % extract relevant data
        try data_in.conmat.experiment = temp.data_in{i_fi}.randmat.experiment;
            data_in.RDK = temp.data_in{i_fi}.RDK.RDK;
        end
        if any(strcmp(fieldnames(temp.data_in{i_fi}.resp),'experiment'))
            temp.index=cell2mat(cellfun(@(x) ~isempty(cell2mat({x(:).trialnumber})), temp.data_in{i_fi}.resp.experiment,'UniformOutput',false));
            data_in.resp.experiment(temp.index)=temp.data_in{i_fi}.resp.experiment(temp.index);
            data_in.button_presses.experiment(temp.index)=temp.data_in{i_fi}.button_presses.experiment(temp.index);
        end
    end
    
    %% get color values for RDKs
    t.colmat = cell2mat(cellfun(@(x) x(1,1:3),{data_in.RDK.col},'UniformOutput',false)');
    t.colmat_base = round(bsxfun(@rdivide, t.colmat,max(t.colmat,[],2)).*10)/10;
    
    
    t.colmat_name = cellfun(@(x) p.colnames(ismember(p.cols, round((x(1:3)./max(x(1:3))).*10)./10,'rows')),...
        num2cell(t.colmat,2));
    t.RDKmat = {'RDK1';'RDK2';'RDK3'};
    t.RDKmat2 = reshape(repmat(t.RDKmat',2,1),numel(t.RDKmat)*2,[]);
    t.colmat_name2 = reshape(repmat(t.colmat_name',2,1),numel(t.colmat_name)*2,[]);
    t.colmat_XYZ = rgb2xyz_custom(t.colmat, col.xy_r, col.xy_g, col.xy_b, col.xyz_w );
    t.colmat_LAB = xyz2lab(t.colmat_XYZ);
    t.colmat_LAB_deg = rad2deg(atan2(t.colmat_LAB(:,3),t.colmat_LAB(:,2)));
    t.colmat_LAB_deg2 = reshape(repmat(t.colmat_LAB_deg',2,1),numel(t.colmat_name)*2,[]);
    t.colmat_base_XYZ = rgb2xyz_custom(t.colmat_base, col.xy_r, col.xy_g, col.xy_b, col.xyz_w );
    t.colmat_base_LAB = xyz2lab(t.colmat_base_XYZ);
    t.attmat_toRDK = repmat({'RDK1';'RDK2'},3,1); % which RDK is attended
    t.attmat_col = strcat('attend',{' '} ,repmat( t.colmat_name(1:2),3,1)); % which color is attended?
    t.attmat_RDK = {'attended';'unattended';'unattended';'attended';'irrelevant';'irrelevant'};
    colordist=[];
    colordist_angle = [];
    colordist_base_angle = [];
    colordist_base_angle2 = [];
    for i_src = 1:3
        for i_dest = 1:3
            % distance=squareroot((L1-L2)^2+(a1-a2)^2+(b1-b2)^2)
            colordist(i_src, i_dest) = sqrt(...
                (t.colmat_LAB(i_src,1)-t.colmat_LAB(i_dest,1))^2 +...
                (t.colmat_LAB(i_src,2)-t.colmat_LAB(i_dest,2))^2 +...
                (t.colmat_LAB(i_src,3)-t.colmat_LAB(i_dest,3))^2 ...
                );
%             % distance in hue angles h° = atan2(b,a)
%             colordist_angle(i_src, i_dest) = abs(rad2deg(...
%                 atan2(t.colmat_LAB(i_src,3),t.colmat_LAB(i_src,2)) - ...
%                 atan2(t.colmat_LAB(i_dest,3),t.colmat_LAB(i_dest,2))));
%             % distance in hue angles h° = atan2(b,a)
%             colordist_base_angle(i_src, i_dest) = abs(rad2deg(...
%                 atan2(t.colmat_base_LAB(i_src,3),t.colmat_base_LAB(i_src,2)) - ...
%                 atan2(t.colmat_base_LAB(i_dest,3),t.colmat_base_LAB(i_dest,2))));
            
            % distance in hue angles h° = atan2(b,a)
            t.dist = abs(atan2(t.colmat_LAB(i_src,3),t.colmat_LAB(i_src,2)) - atan2(t.colmat_LAB(i_dest,3),t.colmat_LAB(i_dest,2)));
            if t.dist >= pi
                colordist_angle(i_src, i_dest) = rad2deg((2*pi)-t.dist);
            else
                colordist_angle(i_src, i_dest) = rad2deg(t.dist);
            end
            
            t.dist = abs(atan2(t.colmat_base_LAB(i_src,3),t.colmat_base_LAB(i_src,2)) - ...
                atan2(t.colmat_base_LAB(i_dest,3),t.colmat_base_LAB(i_dest,2)));
            if t.dist >= pi
                 colordist_base_angle(i_src, i_dest) = rad2deg((2*pi)-t.dist);
            else
                colordist_base_angle(i_src, i_dest) = rad2deg(t.dist);
            end
            
            
 
            
            
        end
    end
    t.surr = [1 2 1];
    t.dist2attended = cellfun(@(x,y) colordist(strcmp(t.RDKmat,x),strcmp(t.RDKmat,y)), t.RDKmat2,t.attmat_toRDK);
    t.dist2unattended = cellfun(@(x,y) colordist(strcmp(t.RDKmat,x),t.surr(find(strcmp(t.RDKmat,y))+1)), t.RDKmat2,t.attmat_toRDK);
    t.dist2irrelevant = cellfun(@(x,y) colordist(strcmp(t.RDKmat,x),3), t.RDKmat2,t.attmat_toRDK);
    t.dist2attended_deg = cellfun(@(x,y) colordist_angle(strcmp(t.RDKmat,x),strcmp(t.RDKmat,y)), t.RDKmat2,t.attmat_toRDK);
    t.dist2unattended_deg = cellfun(@(x,y) colordist_angle(strcmp(t.RDKmat,x),t.surr(find(strcmp(t.RDKmat,y))+1)), t.RDKmat2,t.attmat_toRDK);
    t.dist2irrelevant_deg = cellfun(@(x,y) colordist_angle(strcmp(t.RDKmat,x),3), t.RDKmat2,t.attmat_toRDK);
    t.dist2attended_base_deg = cellfun(@(x,y) colordist_base_angle(strcmp(t.RDKmat,x),strcmp(t.RDKmat,y)), t.RDKmat2,t.attmat_toRDK);
    t.dist2unattended_base_deg = cellfun(@(x,y) colordist_base_angle(strcmp(t.RDKmat,x),t.surr(find(strcmp(t.RDKmat,y))+1)), t.RDKmat2,t.attmat_toRDK);
    t.dist2irrelevant_base_deg = cellfun(@(x,y) colordist_base_angle(strcmp(t.RDKmat,x),3), t.RDKmat2,t.attmat_toRDK);
    
    if i_sub == 1
        colmat = cell2struct([repmat({p.subs2use(i_sub)},numel(t.RDKmat2),1) t.RDKmat2 t.colmat_name2 num2cell(t.colmat_LAB_deg2)....
            t.attmat_toRDK t.attmat_col t.attmat_RDK ...
            num2cell([t.dist2attended t.dist2unattended t.dist2irrelevant t.dist2attended_deg t.dist2unattended_deg t.dist2irrelevant_deg ...
            t.dist2attended_base_deg t.dist2unattended_base_deg t.dist2irrelevant_base_deg])]',...
            {'participant', 'RDK_id', 'RDK_color', 'RDK_color_LAB_deg', 'RDK_attended', 'RDK_col_attended' , 'RDK_att_status'...
            'dist2attended','dist2unattended','dist2irrelevant','dist2attended_deg','dist2unattended_deg','dist2irrelevant_deg',...
            'dist2attended_base_deg','dist2unattended_base_deg','dist2irrelevant_base_deg'});
        colmat_RGB = t.colmat;
        colmat_XYZ = t.colmat_XYZ;
        colmat_lab = t.colmat_LAB;
        colmat_base_rgb = t.colmat_base;
        colmat_base_LAB =t.colmat_base_LAB;
        colmat_name = t.colmat_name;
    else
        colmat = [colmat; ...
            cell2struct([repmat({p.subs2use(i_sub)},numel(t.RDKmat2),1) t.RDKmat2 t.colmat_name2 num2cell(t.colmat_LAB_deg2)....
            t.attmat_toRDK t.attmat_col t.attmat_RDK ...
            num2cell([t.dist2attended t.dist2unattended t.dist2irrelevant t.dist2attended_deg t.dist2unattended_deg t.dist2irrelevant_deg ...
            t.dist2attended_base_deg t.dist2unattended_base_deg t.dist2irrelevant_base_deg])]',...
            {'participant', 'RDK_id', 'RDK_color', 'RDK_color_LAB_deg', 'RDK_attended', 'RDK_col_attended' , 'RDK_att_status'...
            'dist2attended','dist2unattended','dist2irrelevant','dist2attended_deg','dist2unattended_deg','dist2irrelevant_deg',...
            'dist2attended_base_deg','dist2unattended_base_deg','dist2irrelevant_base_deg'})];
        colmat_RGB = [t.colmat; colmat_RGB];
        colmat_XYZ = [t.colmat_XYZ; colmat_XYZ];
        colmat_lab = [colmat_lab; t.colmat_LAB];
        colmat_base_rgb = [colmat_base_rgb; t.colmat_base];
        colmat_base_LAB = [colmat_base_LAB; t.colmat_base_LAB];
        colmat_name = [colmat_name; t.colmat_name];
    end
    
end

%% graphical checking
% figure; histogram([colmat.dist2attended],20)

% plot lab space
% plot_Lab(1,colmat_lab',1,colmat_base_rgb,15,0)
figure; polarhistogram(atan2(colmat_lab(:,3),colmat_lab(:,2)),90)
% figure; histogram([colmat.dist2attended_deg],150)
% figure; histogram([colmat.dist2unattended_deg],150)

figure;
pl.data = [];
for i_col = 1:3
    pl.data(:,i_col) = colmat_lab(strcmp(colmat_name,p.colnames{ i_col}),1);
    histogram(pl.data(:,i_col),20:0.5:30,'FaceColor',p.cols(i_col,:)); hold on
end
title('L values of CIE L*a*b space for different colors');
figure;
subplot(1,3,1); scatter(pl.data(:,1),pl.data(:,2)); xlabel(p.colnames{1});ylabel(p.colnames{2})
[corres.r corres.p] = corr(pl.data(:,1),pl.data(:,2),'Type','Spearman');
title(sprintf('L values; Rho(%1.0f) = %1.3f; p = %1.3f', size(pl.data,1), corres.r, corres.p))
subplot(1,3,2); scatter(pl.data(:,1),pl.data(:,3)); xlabel(p.colnames{1});ylabel(p.colnames{3})
[corres.r, corres.p] = corr(pl.data(:,1),pl.data(:,3),'Type','Spearman');
title(sprintf('L values; Rho(%1.0f) = %1.3f; p = %1.3f', size(pl.data,1), corres.r, corres.p))
subplot(1,3,3); scatter(pl.data(:,2),pl.data(:,3)); xlabel(p.colnames{2});ylabel(p.colnames{3})
[corres.r, corres.p] = corr(pl.data(:,2),pl.data(:,3),'Type','Spearman');
title(sprintf('L values; Rho(%1.0f) = %1.3f; p = %1.3f', size(pl.data,2), corres.r, corres.p))


figure;
pl.data = [];
for i_col = 1:3
    pl.idx = strcmp(colmat_name,p.colnames{ i_col});
    pl.data(:,i_col) = rad2deg(atan2(colmat_lab(pl.idx,3),colmat_lab(pl.idx,2))-mean(atan2(colmat_lab(pl.idx,3),colmat_lab(pl.idx,2))));
    histogram(pl.data(:,i_col),linspace(-0.5,0.5,30),'FaceColor',p.cols(i_col,:)); hold on
end
title('hue angle differences of CIE L*a*b space for different colors');


figure;
pl.data = [];
for i_col = 1:3
    pl.idx = strcmp(colmat_name,p.colnames{ i_col});
    pl.data(:,i_col) = atan2(colmat_lab(pl.idx,3),colmat_lab(pl.idx,2));
    polarhistogram(pl.data(:,i_col),'BinWidth',2*pi/90,'FaceColor',p.cols(i_col,:)); hold on
end
title('hue angles of CIE L*a*b space for isoluminant colors');

figure;
pl.data = [];
for i_col = 1:3
    pl.idx = strcmp(colmat_name,p.colnames{ i_col});
    pl.data(:,i_col) = atan2(colmat_base_LAB(pl.idx,3),colmat_base_LAB(pl.idx,2));
    polarhistogram(pl.data(:,i_col),'BinWidth',2*pi/90,'FaceColor',p.cols(i_col,:)); hold on
end
title('hue angles of CIE L*a*b space for base colors');


pl.data_m = []; pl.data_std = [];
for i_col = 1:3
    pl.data_m(i_col,:) = mean(colmat_lab(strcmp(colmat_name,p.colnames{ i_col}),:));
    pl.data_std(i_col,:) = std(colmat_lab(strcmp(colmat_name,p.colnames{ i_col}),:));
    fprintf('\n%s Mean(Std) of CIE L*a*b values: %1.3f (%1.3f), %1.3f (%1.3f), %1.3f (%1.3f)',...
        p.colnames{i_col}, reshape([pl.data_m(i_col,:);pl.data_std(i_col,:)],1,[]))
end
fprintf('\n')

% calculate color differences
pl.data_m = [];
for i_col = 1:3
    pl.data_m(i_col,:) = mean(colmat_lab(strcmp(colmat_name,p.colnames{ i_col}),:));
end
pl.data_m_hue = atan2(pl.data_m(:,3),pl.data_m(:,2));
pl.data_m_hue(pl.data_m_hue<0) = 2*pi -abs(pl.data_m_hue(pl.data_m_hue<0));
pl.hue_diffs = abs([pl.data_m_hue(1)-pl.data_m_hue(2) pl.data_m_hue(1)-pl.data_m_hue(3) pl.data_m_hue(2)-pl.data_m_hue(3)]);
pl.hue_diffs(pl.hue_diffs>pi)=(2*pi)-pl.hue_diffs(pl.hue_diffs>pi);
pl.hue_diffs = rad2deg(pl.hue_diffs);
fprintf('range of differences between color hues: %1.3f° %1.3f° %1.3f°\n',pl.hue_diffs)


plot_colorwheel([0 1 0; 0 0.4 1; 1 0.4 0],'ColorSpace','propixxrgb',...
    'SavePath','C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_FShiftBase\figures',...
    'LAB_L',50,'NumSegments',60,'AlphaColWheel',1,'LumBackground',100)

plot_colorwheel([0 1 0; 0 0.4 1; 1 0.4 0],'ColorSpace','propixxrgb',...
    'LAB_L',50,'NumSegments',60,'AlphaColWheel',1,'LumBackground',100,'disp_LAB_vals',true)
plot_colorwheel([0    0.0580    0.1451;     0.0941    0.0376         0;          0    0.0784         0],'ColorSpace','propixxrgb',...
    'LAB_L',50,'NumSegments',60,'AlphaColWheel',1,'LumBackground',100,'disp_LAB_vals',true)

% plot_colorwheel([0 1 0; 0 0.4 1; 1 0.4 0],'ColorSpace','propixxrgb',...
%     'SavePath','C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_FShiftBase\figures',...
%     'LAB_L',75,'NumSegments',60,'AlphaColWheel',1,'LumBackground',100)


%% saving
colmat_table = struct2table(colmat);
% general summary mean RT and hit rate
% export data
p.path = 'C:\Users\psy05cvd\Dropbox\work\R-statistics\experiments\ssvep_FShiftBase\data';
p.filename1 = 'color_distances.csv';
writetable(colmat_table, fullfile(p.path,p.filename1))







