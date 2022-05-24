%% File processing
clc;clear;close all

addpath('') % SupplFunc necessitated
xlsmat = dir('*.xlsx');
csvmat = dir('*.csv'); 
numfiles_xlsx = length(xlsmat);
numfiles_csv = length(csvmat); 

if numfiles_xlsx == numfiles_csv
    numfiles = numfiles_xlsx;
else
    msgbox('Error : Difference between number of files');
    pause;
end

xlsdata = cell(1, numfiles);
csvdata = cell(1, numfiles);
raw = cell(1,numfiles);
raw_CA = cell(1,numfiles);
for k = 1:numfiles
    xlsdata{k} = readmatrix(xlsmat(k).name);
    csvdata{k} = readmatrix(csvmat(k).name);
    raw{k} = readcell(xlsmat(k).name);
    raw_CA{k} = readcell(csvmat(k).name);
end

RS_allmice_total = {};

for lazy_filename_int = 1:length(xlsmat)
    tagged_data = xlsdata{lazy_filename_int};
    neuron_data = csvdata{lazy_filename_int};
    
    time_shift = 0;

    tagged_data = tagged_data-time_shift;
    
    Positive_neuron_data = neuron_data;
    Positive_neuron_data(Positive_neuron_data <= 0) = 0;
    Negative_neuron_data = neuron_data;
    Negative_neuron_data(Negative_neuron_data >= 0) = 0;
    
    tagged_data_size = size(tagged_data);
    neuron_data_size = size(neuron_data);
    
    RS_total = [];
    Cell_id_increased = [];
    Cell_id_decreased = [];
    Cell_id_neutral = [];
    
    Behavior_Classification = zeros(3,3);
    fprintf('\nBEHAVIOR OPTIONS: \n');
    behavior_size = tagged_data_size(2)/2;
    for behavior_it = 1:behavior_size
        fprintf(int2str(behavior_it) + ". " + convertCharsToStrings(raw{lazy_filename_int}{1,2*behavior_it-1}) + '\n');
    end
    
    for lazy_behavior_int = 1:tagged_data_size(2)/2
        behavior_size = tagged_data_size(2)/2;
        Missing = cellfun(@(x) ismissing(x),raw{lazy_filename_int}(:,2*lazy_behavior_int-1),'UniformOutput',false);
        Missing_select = cell2mat(Missing(3));
        nb = 0;
        if Missing_select
           continue;
           nb = nb + 1;
        else
           behavior = lazy_behavior_int;
        end

        behavior_fullname = raw{lazy_filename_int}{1,2*behavior-1};
        behavior_name=string(behavior_fullname);
        fprintf("\nSelected Behavior: " + behavior_fullname + '\n');
        
        r = 2*behavior-1;

        behavior_data  = zeros(neuron_data_size(1), neuron_data_size(2));

        Time = neuron_data(:, 1);
        Time = Time - Time(1);

        pass_behavior = size(tagged_data);
        
        if pass_behavior(2) > r
            behavior_data_size = find(tagged_data(:, r) > 0, 1, 'last');
        else 
            continue;
        end        
        
        if isempty(behavior_data_size) 
            behavior_data_size = 0; end
        
        Behavior_start_stop = zeros(tagged_data_size(1), 2);

        epochs = zeros(neuron_data_size(1), 1);
        Behaviori = 1;
        
        for iterator1 = 1:tagged_data_size(1)
            Behavior_start = find(Time >= tagged_data(iterator1, r), 1, 'first');
            Behavior_stop = find(Time >= tagged_data(iterator1, r+1), 1, 'first');
            if (iterator1 <= behavior_data_size && ~isempty(Behavior_start) && ~isempty(Behavior_stop) && Behavior_start~=1 && Behavior_stop ~=1)
                Behavior_start_stop(Behaviori, 1) = Behavior_start;
                Behavior_start_stop(Behaviori, 2) = Behavior_stop;
                Behaviori = Behaviori+1;
            end
            behavior_data(Behavior_start:Behavior_stop, 1:neuron_data_size(2)) = 1; 
            epochs(Behavior_start:Behavior_stop, 1) = iterator1;
        end
        
        behavior_data_size = find(Behavior_start_stop(:,1) > 0, 1, 'last');
        
        if isempty(behavior_data_size)
            behavior_data_size = 0; 
        end
        
        Behavior_start_stop(behavior_data_size+1:end,:) = [];
        for i = 1:size(Behavior_start_stop,1)
            frame_start = Behavior_start_stop(i,1);
            frame_stop = Behavior_start_stop(i,2);
            SVM_input_data(lazy_filename_int,lazy_behavior_int,i)= {neuron_data(frame_start:frame_stop,:)};
        end        
    end
end
%% Data control
% Default setting in inscopix 20 frame/s 
% 1sec = 20 in Extraction time
% SVMdataExtraction function need Data,Device_Framerate and Time to see.
% SVM_data_input = SVMdataExtraction(Data,Device_Framerate,Time) 
Inputdata_trimmed = SVMdataExtraction(SVM_input_data,20,2);
% Inputdata_trimmed = SVMdataExtraction(SVM_input_data,20,0.5);
% Inputdata_trimmed = SVMdataExtraction(SVM_input_data,20,1);
longest_Dur = longestbehaviorsignal(size(SVM_input_data),SVM_input_data);

Inputdata_stretched = stretchInputSVM(SVM_input_data,20,1,0.7);
Inputdata_trimmed_signal_raw = Inputdata_trimmed(:,1:end-1);
Inputdata_trimmed_label_behavior = Inputdata_trimmed(:,end);
Inputdata_stretched_signal_raw = Inputdata_stretched(:,1:end-1);
Inputdata_stretched_label_behavior = Inputdata_stretched(:,end);

shuffle_ratio = 0.70;
num_row = size(Inputdata_trimmed_signal_raw,1);
idx_row = randperm(num_row);

feat_WT_tr_trimmed = Inputdata_trimmed_signal_raw(idx_row(1:round(shuffle_ratio*num_row)),:);
feat_WT_te_trimmed = Inputdata_trimmed_signal_raw(idx_row(round(shuffle_ratio*num_row)+1:end),:);
target_WT_tr_trimmed = Inputdata_trimmed_label_behavior(idx_row(1:round(shuffle_ratio*num_row)),:);
target_WT_te_trimmed = Inputdata_trimmed_label_behavior(idx_row(round(shuffle_ratio*num_row)+1:end),:);

shuffle_ratio = 0.70;
num_row = size(Inputdata_stretched_signal_raw,1);
idx_row = randperm(num_row);

feat_WT_tr_stretched = Inputdata_stretched_signal_raw(idx_row(1:round(shuffle_ratio*num_row)),:);
feat_WT_te_stretched = Inputdata_stretched_signal_raw(idx_row(round(shuffle_ratio*num_row)+1:end),:);
target_WT_tr_stretched = Inputdata_stretched_label_behavior(idx_row(1:round(shuffle_ratio*num_row)),:);
target_WT_te_stretched = Inputdata_stretched_label_behavior(idx_row(round(shuffle_ratio*num_row)+1:end),:);
%% Support vector regression (SVR) simulation
% [SVR_WT_target_te_trimmed,SVR_WT_y_fit_trimmed] = simulateSVR(feat_WT_tr_trimmed,feat_WT_te_trimmed,target_WT_tr_trimmed,target_WT_te_trimmed);
% [SVR_WT_target_te_stretched,SVR_WT_y_fit_stretched] = simulateSVR(feat_WT_tr_stretched,feat_WT_te_stretched,target_WT_tr_stretched,target_WT_te_stretched);
Mdl = fitcecoc(feat_WT_tr_stretched,categorical(target_WT_tr_stretched));
Value = predict(Mdl,feat_WT_te_stretched);
plotconfusion(Value,categorical(target_WT_te_stretched))