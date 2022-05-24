%% Calcium normalization

clc;clear;close all

files = dir('inputSVM_raw_WT*.csv');
num_files = length(files);
data = cell(1,num_files);
for idx_files = 1:num_files
    inputSVM_raw{idx_files} = readcell(files(idx_files).name);
    inputSVM_normalized = inputSVM_raw{1,idx_files};
    num_neuron = size((inputSVM_normalized(:,2:end-1)),2);
    for idx_neuron = 1:num_neuron
        calcium_neuron = [inputSVM_normalized{:,idx_neuron+1}];
        ceiling_neuron = max(calcium_neuron);
        normalization_neuron = calcium_neuron/ceiling_neuron;
        normalization_neuron = normalization_neuron';
        inputSVM_normalized(:,idx_neuron+1) = num2cell(normalization_neuron);
    end
    inputSVM_normalized = cell2mat(inputSVM_normalized);
    csvwrite(strcat(files(idx_files).name(1:end-4),'_normalized.csv'),inputSVM_normalized)
end