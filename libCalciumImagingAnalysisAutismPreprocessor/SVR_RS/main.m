% File processing
clc;clear;close all
PATH = '/Users/shong/dev/ASD/ASD_v10Feb2020/src/DB/DataASD_LabNCA/MableBurying/SVR/RS/';

FILE_LIST = dir(fullfile(PATH,'CNTNAP2KO*.xlsx'));
LIST_L_DATA = length(FILE_LIST(:, 1));
CNTNAP2KO = [];
for file = 1:LIST_L_DATA
    FILE_FULLNAME = FILE_LIST(file).name;
    opts = detectImportOptions([PATH FILE_FULLNAME]);
    CNTNAP2KO_new = readtable([PATH FILE_FULLNAME],opts);
    CNTNAP2KO = [CNTNAP2KO;CNTNAP2KO_new];
end
CNTNAP2KO_RS = CNTNAP2KO(:,1);
CNTNAP2KO_label = CNTNAP2KO(:,3);

FILE_LIST = dir(fullfile(PATH,'SHANK3KO*.xlsx'));
LIST_L_DATA = length(FILE_LIST(:, 1));
SHANK3KO = [];
for file = 1:LIST_L_DATA
    FILE_FULLNAME = FILE_LIST(file).name;
    opts = detectImportOptions([PATH FILE_FULLNAME]);
    SHANK3KO_new = readtable([PATH FILE_FULLNAME],opts);
    SHANK3KO = [SHANK3KO;SHANK3KO_new];
end
SHANK3KO_RS = SHANK3KO(:,1);
SHANK3KO_label = SHANK3KO(:,3);

FILE_LIST = dir(fullfile(PATH,'WT*.xlsx'));
LIST_L_DATA = length(FILE_LIST(:, 1));
WT = [];
for file = 1:LIST_L_DATA
    FILE_FULLNAME = FILE_LIST(file).name;
    opts = detectImportOptions([PATH FILE_FULLNAME]);
    WT_new = readtable([PATH FILE_FULLNAME],opts);
    WT = [WT;WT_new];
end
WT_RS = WT(:,1);
WT_label = WT(:,3);

key_behavior = {'Grooming', 'Digging', 'Rearing'};

for idx_behavior = 1:numel(key_behavior)
    [row_behavior,column_behavior] = find(strcmpi(CNTNAP2KO_label{:,:},key_behavior{idx_behavior}));
    CNTNAP2KO_label{row_behavior,column_behavior} = {idx_behavior};
    [row_behavior,column_behavior] = find(strcmpi(SHANK3KO_label{:,:},key_behavior{idx_behavior}));
    SHANK3KO_label{row_behavior,column_behavior} = {idx_behavior};
    [row_behavior,column_behavior] = find(strcmpi(WT_label{:,:},key_behavior{idx_behavior}));
    WT_label{row_behavior,column_behavior} = {idx_behavior};
end

shuffle_ratio = 0.70;

num_row = size(CNTNAP2KO_RS,1);
idx_row = randperm(num_row);

feat_CNTNAP2KO_tr = CNTNAP2KO_RS(idx_row(1:round(shuffle_ratio*num_row)),:);feat_CNTNAP2KO_tr = table2array(feat_CNTNAP2KO_tr);csvwrite('feat_CNTNAP2KO_tr.csv',feat_CNTNAP2KO_tr);
feat_CNTNAP2KO_te = CNTNAP2KO_RS(idx_row(round(shuffle_ratio*num_row)+1:end),:);feat_CNTNAP2KO_te = table2array(feat_CNTNAP2KO_te);csvwrite('feat_CNTNAP2KO_te.csv',feat_CNTNAP2KO_te);
target_CNTNAP2KO_tr = CNTNAP2KO_label(idx_row(1:round(shuffle_ratio*num_row)),:);target_CNTNAP2KO_tr = cell2mat(table2array(target_CNTNAP2KO_tr));csvwrite('target_CNTNAP2KO_tr.csv',target_CNTNAP2KO_tr);
target_CNTNAP2KO_te = CNTNAP2KO_label(idx_row(round(shuffle_ratio*num_row)+1:end),:);target_CNTNAP2KO_te = cell2mat(table2array(target_CNTNAP2KO_te));csvwrite('target_CNTNAP2KO_te.csv',target_CNTNAP2KO_te);


num_row = size(SHANK3KO_RS,1);
idx_row = randperm(num_row);

feat_SHANK3KO_tr = SHANK3KO_RS(idx_row(1:round(shuffle_ratio*num_row)),:);feat_SHANK3KO_tr = table2array(feat_SHANK3KO_tr);csvwrite('feat_SHANK3KO_tr.csv',feat_SHANK3KO_tr);
feat_SHANK3KO_te = SHANK3KO_RS(idx_row(round(shuffle_ratio*num_row)+1:end),:);feat_SHANK3KO_te = table2array(feat_SHANK3KO_te);csvwrite('feat_SHANK3KO_te.csv',feat_SHANK3KO_te);
target_SHANK3KO_tr = SHANK3KO_label(idx_row(1:round(shuffle_ratio*num_row)),:);target_SHANK3KO_tr = cell2mat(table2array(target_SHANK3KO_tr));csvwrite('target_SHANK3KO_tr.csv',target_SHANK3KO_tr);
target_SHANK3KO_te = SHANK3KO_label(idx_row(round(shuffle_ratio*num_row)+1:end),:);target_SHANK3KO_te = cell2mat(table2array(target_SHANK3KO_te));csvwrite('target_SHANK3KO_te.csv',target_SHANK3KO_te);


num_row = size(WT_RS,1);
idx_row = randperm(num_row);

feat_WT_tr = WT_RS(idx_row(1:round(shuffle_ratio*num_row)),:);feat_WT_tr = table2array(feat_WT_tr);csvwrite('feat_WT_tr.csv',feat_WT_tr);
feat_WT_te = WT_RS(idx_row(round(shuffle_ratio*num_row)+1:end),:);feat_WT_te = table2array(feat_WT_te);csvwrite('feat_WT_te.csv',feat_WT_te);
target_WT_tr = WT_label(idx_row(1:round(shuffle_ratio*num_row)),:);target_WT_tr = cell2mat(table2array(target_WT_tr));csvwrite('target_WT_tr.csv',target_WT_tr);
target_WT_te = WT_label(idx_row(round(shuffle_ratio*num_row)+1:end),:);target_WT_te = cell2mat(table2array(target_WT_te));csvwrite('target_WT_te.csv',target_WT_te);

% Support vector regression (SVR) simulation
[SVR_WT_target_te,SVR_WT_y_fit] = simulateSVR(feat_WT_tr,feat_WT_te,target_WT_tr,target_WT_te);%figure();plotconfusion(SVR_WT_target_te,SVR_WT_y_fit)
[SVR_CNTNAP2KO_target_te,SVR_CNTNAP2KO_y_fit] = simulateSVR(feat_CNTNAP2KO_tr,feat_CNTNAP2KO_te,target_CNTNAP2KO_tr,target_CNTNAP2KO_te);%figure();plotconfusion(SVR_CNTNAP2KO_target_te,SVR_CNTNAP2KO_y_fit)
[SVR_SHANK3KO_target_te,SVR_SHANK3KO_y_fit] = simulateSVR(feat_SHANK3KO_tr,feat_SHANK3KO_te,target_SHANK3KO_tr,target_SHANK3KO_te);%figure();plotconfusion(SVR_SHANK3KO_target_te,SVR_SHANK3KO_y_fit)