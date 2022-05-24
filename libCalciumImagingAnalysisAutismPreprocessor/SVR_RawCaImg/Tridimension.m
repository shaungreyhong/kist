%% Tridimensional vector
% x: Temporal frame domain
% y: Calcium response domain
% z: Behavior label (1-4) domain

clc;clear;close all

files = dir('inputSVM_raw_WT*.csv');
num_files = length(files);
for idx_files = 1:num_files
    inputSVM3D{idx_files} = readcell(files(idx_files).name);
    inputSVM3D_mouse = inputSVM3D{1,idx_files};
    [num_frame,num_neuron] = size((inputSVM3D_mouse(:,2:end-1)));
    for idx_neuron = 1:num_neuron
        for idx_frame = 1:num_frame
            inputSVM3D_mouse_voxel = cell(1,3);
            inputSVM3D_mouse_voxel{1,1} = inputSVM3D_mouse{idx_frame,1};
            inputSVM3D_mouse_voxel{1,2} = inputSVM3D_mouse{idx_frame,idx_neuron+1};
            inputSVM3D_mouse_voxel{1,3} = inputSVM3D_mouse{idx_frame,end};
            inputSVM3D_mouse{idx_frame,idx_neuron+1} = inputSVM3D_mouse_voxel;
        end
    end
    inputSVM3D_mouse(:,1) = [];
    inputSVM3D_mouse(:,end) = [];
    save(strcat(files(idx_files).name(1:end-4),'_3D.mat'),'inputSVM3D_mouse')
end
%% Merge normalized WT data

clc;clear;close all

files = dir('*CA_normalized_3D.mat');
inputSVM3D_normalization = [];
num_files = length(files);
for idx_files = 1:num_files
    inputSVM3D{idx_files} = load(files(idx_files).name);
    inputSVM3D_normalization = [inputSVM3D_normalization inputSVM3D{idx_files}.inputSVM3D_mouse(1:10800,:)]; % 9*60*20 = 10800
end
%% Merge raw WT data

clc;clear;close all

files = dir('*CA_3D.mat');
inputSVM3D_raw = [];
num_files = length(files);
for idx_files = 1:num_files
    inputSVM3D{idx_files} = load(files(idx_files).name);
    inputSVM3D_raw = [inputSVM3D_raw inputSVM3D{idx_files}.inputSVM3D_mouse(1:10800,:)]; % 9*60*20 = 10800
end