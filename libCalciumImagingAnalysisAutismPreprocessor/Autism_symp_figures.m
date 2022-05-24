load('data_day1.mat');
load('filtered_data_day1.mat');
load('mb_start_stop_day1.mat');
load('r_start_stop_day1.mat');
load('Time.mat');
neuronNum = length(filtered_data(1,:))/2;

%% plot raster using only above or under threshold

R_h_labels = plotData3Threshold('Rearing',Rearing_start_stop, neuron_data, Time, neuronNum, 1);
R_m_labels = plotData3Threshold('Rearing',Rearing_start_stop, neuron_data, Time, neuronNum, 2);
R_l_labels = plotData3Threshold('Rearing',Rearing_start_stop, neuron_data, Time, neuronNum, 3);

%%

MB_h_labels = plotData3Threshold('Marble burying',MB_start_stop, neuron_data, Time, neuronNum, 1);
MB_m_labels = plotData3Threshold('Marble burying',MB_start_stop, neuron_data, Time, neuronNum, 2);
MB_l_labels = plotData3Threshold('Marble burying',MB_start_stop, neuron_data, Time, neuronNum, 3);

%% plot rater plot of signals from -1 to 1 seconds of epochs(behaviors)
% 20 ~70%
MB_labels =  plotFilteredData('Marble burying',MB_start_stop, neuron_data, Time, neuronNum, 0);
MB_flabels = plotFilteredData('Marble burying',MB_start_stop, filtered_data, Time, neuronNum, 1);


%%
% 40~70 %
R_labels = plotFilteredData('Rearing',Rearing_start_stop, neuron_data, Time, neuronNum, 0);
R_flabels = plotFilteredData('Rearing',Rearing_start_stop, filtered_data, Time, neuronNum, 1);


%%
% MB Response Strenths
Positive_neuron_data = MB_filtered_data_correlation_factors;
Positive_neuron_data(Positive_neuron_data <= 0) = 0;
Negative_neuron_data = MB_filtered_data_correlation_factors;
Negative_neuron_data(Negative_neuron_data >= 0) = 0;


figure
hold on
bar(Positive_neuron_data, 'r')
bar(Negative_neuron_data, 'b')

title('MB Response Strenths')
xlabel('Neuron Channels')
ylabel('Response Strength')
ylim([-1 1])

% Rearing Response Strenths
Positive_neuron_data = Rearing_filtered_data_correlation_factors;
Positive_neuron_data(Positive_neuron_data <= 0) = 0;
Negative_neuron_data = Rearing_filtered_data_correlation_factors;
Negative_neuron_data(Negative_neuron_data >= 0) = 0;

figure
hold on
bar(Positive_neuron_data, 'r')
bar(Negative_neuron_data, 'b')

title('Rearing Response Strenths')
xlabel('Neuron Channels')
ylabel('Response Strength')
ylim([-1 1])