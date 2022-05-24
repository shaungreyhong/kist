% Extract behavior-labeled raw data
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

neuron_data_all = [];
for lazy_filename_int = 1:length(xlsmat)
    tagged_data = xlsdata{lazy_filename_int};
    neuron_data = csvdata{lazy_filename_int};
    
    time_shift = 0;
    tagged_data = tagged_data-time_shift;
    tagged_data_size = size(tagged_data);
    neuron_data_size = size(neuron_data);
    neuron_data(:,neuron_data_size(2)+1) = 4; % 0 or 4 for "Nothing"
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
%                 
%                 size_gold_time = size(find(Behavior_start_stop<12000),1);
%                 label_behavior = [];
%                 key_behavior = {'Grooming', 'Digging', 'Rearing', 'Nothing'};
%                 for idx_gold_time = 1:12000 % Gold time 10 miniutes: 600s * 20frames/s * 12000samples
%                    for idx_gold_time_2 = 1:size_gold_time
%                        if idx_gold_time < Behavior_start_stop
%                            for idx_behavior = 1:numel(key_behavior)
%                                [row_behavior,column_behavior] = find(strcmpi(neuron_data(:,:),key_behavior{idx_behavior}));
%                                neuron_data(row_behavior,end+1) = idx_behavior;
%                            end
%                        else
%                            idx_behavior = 4;
%                            [row_behavior,column_behavior] = find(strcmpi(neuron_data(:,:),key_behavior{idx_behavior}));
%                            neuron_data(row_behavior,end+1) = idx_behavior;
%                        end
%                    end
%                 end               
                Behaviori = Behaviori+1;
            end
            behavior_data(Behavior_start:Behavior_stop, 1:neuron_data_size(2)) = 1; 
            epochs(Behavior_start:Behavior_stop, 1) = iterator1;
        end
        
        % Labeling in whole CA data for bhea
        Behavior_matching = round(tagged_data*20);       

       B_start_stop = sum(~isnan(Behavior_matching(:,2*lazy_behavior_int-1)));
       Behavior_setting_number = size(neuron_data,1);
       for j = 1:Behavior_setting_number
           for k = 1:B_start_stop
               if j >= Behavior_matching(k,2*lazy_behavior_int-1) && j <= Behavior_matching(k,2*lazy_behavior_int)
                   neuron_data(j,end) = lazy_behavior_int;
               end
           end
       end

        behavior_data_size = find(Behavior_start_stop(:,1) > 0, 1, 'last');
        
        if isempty(behavior_data_size)
            behavior_data_size = 0; 
        end
        
        Behavior_start_stop(behavior_data_size+1:end,:) = [];
        for i = 1:size(Behavior_start_stop,1)
            frame_start = Behavior_start_stop(i,1);
            frame_stop = Behavior_start_stop(i,2);
            SVM_input_data(lazy_filename_int,lazy_behavior_int,i) = {neuron_data_all(:,:)};
            % SVM_input_data(lazy_filename_int,lazy_behavior_int,i) = {neuron_data(frame_start:frame_stop,:)};
        end
        save('WT_B629_MB_CA.mat','neuron_data')
        % neuron_data_all = [neuron_data_all;neuron_data(1:12000,2:end)'];
    end
end