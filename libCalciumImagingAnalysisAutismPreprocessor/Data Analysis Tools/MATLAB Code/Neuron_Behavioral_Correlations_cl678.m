
T_Filename = '1st Day Tagged Data';
N_Filename = '1st Day Calcium Cell Traces.xlsx';

% %Read in the data from the tagged data
% T_Filename = '2nd Day Tagged Data';
% N_Filename = '2nd Day Calcium Cell Traces.xlsx';
              
% Create raw data matrix
tagged_data = xlsread(T_Filename);
neuron_data = xlsread(N_Filename);

% Create positive and negative neuron data for positive and negative
% correlation
Positive_neuron_data = neuron_data;
Positive_neuron_data(Positive_neuron_data <= 0) = 0;
Negative_neuron_data = neuron_data;
Negative_neuron_data(Negative_neuron_data >= 0) = 0;

% Store the sizes of the data
tagged_data_size = size(tagged_data);
neuron_data_size = size(neuron_data);

%%

% Initialize Rearing_data and MB_data to the same size as the neuron data
Rearing_data = zeros(neuron_data_size(1), neuron_data_size(2));
MB_data = zeros(neuron_data_size(1), neuron_data_size(2));

% Create an array complete time data increments from neuron data
Time = neuron_data(:, 1);

% Initialize arrays that store the data point numbers for Rearing or MB
% start and stop
Rearing_data_size = find(tagged_data(:, 1) >= 0, 1, 'last');
MB_data_size = find(tagged_data(:, 3) >= 0, 1, 'last');
Rearing_start_stop = zeros(tagged_data_size(1), 2);
MB_start_stop = zeros(tagged_data_size(1), 2);

% Initialize an array, the same size as the neuron data, that stores
% increasing integer values for each behavior window
epochs = zeros(neuron_data_size(1), 1);

% Iterate through tagged data and fill in ones during rearing or MB in
% their respective data sets (Rearing_data & MB_data)
for iterator1 = 1:tagged_data_size(1)
    Rearing_start = find(Time >= tagged_data(iterator1, 1), 1, 'first');
    Rearing_stop = find(Time >= tagged_data(iterator1, 2), 1, 'first');
       
    MB_start = find(Time >= tagged_data(iterator1, 3), 1, 'first');
    MB_stop = find(Time >= tagged_data(iterator1, 4), 1, 'first');

    % Records the stop and start data points
    if (iterator1 <= Rearing_data_size)
        Rearing_start_stop(iterator1, 1) = Rearing_start;
        Rearing_start_stop(iterator1, 2) = Rearing_stop;
    end
    if (iterator1 <= MB_data_size)
        MB_start_stop(iterator1, 1) = MB_start;
        MB_start_stop(iterator1, 2) = MB_stop;
    end
    
    % Fills the cells with 1s in the given time interval
    Rearing_data(Rearing_start:Rearing_stop, 1:neuron_data_size(2)) = 1; 
    MB_data(MB_start:MB_stop, 1:neuron_data_size(2)) = 1; 
    
    % Fills the 'epochs' array with the intended values
    epochs(Rearing_start:Rearing_stop, 1) = iterator1*2 - 1;
    epochs(MB_start:MB_stop, 1) = iterator1*2;
end

% Normalize each neuron channel data
% for iterator = 2:neuron_data_size(2)
%    neuron_data(:, iterator) = neuron_data(:, iterator)/max(neuron_data(:, iterator));
% end

%%
% Plot all the neuron data channels on the same graph, normalized and
% separated
figure (1)
Normalizer = max(max(neuron_data));
plot((Time - Time(1)), neuron_data(:, 2)/Normalizer + 1)
hold on
for channels = 3:neuron_data_size(2)
    plot((Time - Time(1)), neuron_data(:, channels)/Normalizer + channels - 1)
end

% Set axis size based on size of data
axis([0 (Time(end) - Time(1)) 0 neuron_data_size(2)])

% Fill in the areas on the graph (Green --> MB, Red --> Rearing)
area((Time - Time(1)), MB_data(:, 1)*neuron_data_size(2),'FaceColor','g', 'FaceAlpha',.2,'EdgeAlpha', 0.001)
area((Time - Time(1)), Rearing_data(:, 1)*neuron_data_size(2),'FaceColor','r', 'FaceAlpha',.2,'EdgeAlpha', 0.001)
 
xlabel('Time (sec)')
ylabel('Neuron Channels')


%%
% Filter the neuron data using Fast Online Deconvolution of Calcium Imaging Data

% User need to install or setup to run this code
% Friedrich-deconvolveCa: https://github.com/zhoupc/OASIS_matlab
% go to sub directory called OASIS and run ">> setup" if it shows "Error in
% deconvolveCa (line 70) options.sn = GetSn(y);"

filtered_data = zeros(neuron_data_size(1), 2*(neuron_data_size(2) - 1));
for iterator2 = 1:(neuron_data_size(2) - 1)
   [filtered_data(:, iterator2*2 - 1), filtered_data(:, iterator2*2)] = ...
       deconvolveCa(neuron_data(:, iterator2 + 1), 'ar1', 'foopsi', 'optimize_pars');
end

%Plot the filtered data just like the original data
figure (2)
Normalizer = max(max(filtered_data(:, 1:2:end)));
plot((Time - Time(1)), filtered_data(:, 1)/Normalizer + 1)
hold on
for iterator3 = 2:(neuron_data_size(2) - 1)
   plot((Time - Time(1)),filtered_data(:, iterator3*2 - 1)/Normalizer + iterator3) 
end

area((Time - Time(1)), MB_data(:, 1)*neuron_data_size(2),'FaceColor','g', 'FaceAlpha',.2,'EdgeAlpha', 0.001)
area((Time - Time(1)), Rearing_data(:, 1)*neuron_data_size(2),'FaceColor','r', 'FaceAlpha',.2,'EdgeAlpha', 0.001)

% Set axis size based on size of data
axis([0 (Time(end) - Time(1)) 0 neuron_data_size(2)])

xlabel('Time (sec)')
ylabel('Neuron Channels')

%%
% Compute and plot the ROC curves for the orignal and filtered data

% Create arrays to store the auROC data
auROC_neuron_rearing_data = zeros(1, neuron_data_size(2) - 1);
auROC_neuron_MB_data = zeros(1, neuron_data_size(2) - 1);

auROC_neuron_rearing_filtered_data = zeros(1, neuron_data_size(2) - 1);
auROC_neuron_MB_filtered_data = zeros(1, neuron_data_size(2) - 1);

% neuron_data_size(2) looping for cells
for iterator4 = 2:neuron_data_size(2)
    [~, ~, ~, auROC_neuron_rearing_data(1, iterator4 - 1)]...
        = perfcurve(Rearing_data(:, 1), neuron_data(:, iterator4), 1);
    
    [~, ~, ~, auROC_neuron_MB_data(1, iterator4 - 1)]...
        = perfcurve(MB_data(:, 1), neuron_data(:, iterator4), 1);
    

    [~, ~, ~, auROC_neuron_rearing_filtered_data(1, iterator4 - 1)]...
        = perfcurve(Rearing_data(:, 1), filtered_data(:, (iterator4 - 1)*2 - 1), 1);
    
    [~, ~, ~, auROC_neuron_MB_filtered_data(1, iterator4 - 1)]...
        = perfcurve(MB_data(:, 1), filtered_data(:, (iterator4 - 1)*2 - 1), 1);
end

% ROC plot
    N_cell_number = 2; % iput cell number here to plot ROC curve (cell number ends at neuron_data_size(2) "34"
    figure(999)
    [x, y, ~, auROC_neuron_MB_filtered_data(1, N_cell_number - 1)]...
        = perfcurve(MB_data(:, 1), filtered_data(:, (N_cell_number - 1)*2 - 1), 1);
    plot(x,y)

% Calculate correlation factors
MB_data_correlation_factors = (auROC_neuron_MB_data - 0.5)*2;
Rearing_data_correlation_factors = (auROC_neuron_rearing_data - 0.5)*2;

Rearing_filtered_data_correlation_factors = (auROC_neuron_rearing_filtered_data - 0.5)*2;
MB_filtered_data_correlation_factors = (auROC_neuron_MB_filtered_data - 0.5)*2;

figure(3)
stem(MB_data_correlation_factors)
hold on
stem(MB_filtered_data_correlation_factors)

title('MB Response Strenths')
xlabel('Neuron Channels')
ylabel('Response Strength')
legend('Raw Data', 'Filtered Data')
axis([0 (neuron_data_size(2) - 1) -1 1])

figure(4)
stem(Rearing_data_correlation_factors)
hold on
stem(Rearing_filtered_data_correlation_factors)

title('Rearing Response Strengths')
xlabel('Neuron Channels')
ylabel('Response Strength')
legend('Raw Data', 'Filtered Data')
axis([0 (neuron_data_size(2) - 1) -1 1])


%%
% Apply Linear Discriminant Analysis to the neuron data. Then scatter plot
% the neuron channels based on the eigenvectors determined from the LDA

% Creating an array of simply the denoised filtered data
filtered = 0;
if(filtered)
    lda_data = filtered_data;
    lda_data(:, 2:2:end) = []; 
else
    lda_data = neuron_data(:, 2:end);
%     lda_data(:, 2:2:end) = [];
end
% Creating the "Labels" for the LDA (0 --> No behavior, 1 --> Rearing, 
% 2 --> MB)
Labels = (Rearing_data(:, 1) + 2*MB_data(:, 1));
Labels(Labels == 3) = 0;

% Apply the built in MATLAB LDA 
% Creates Figures 6 and 7 which are visualizations of the optimization
% parameters
Mdl = fitcdiscr(lda_data, Labels, 'OptimizeHyperparameters','auto', 'DiscrimType', 'pseudolinear');

% Determine the eigenvector and corresponding eigenvalues from the LDA
% model. Then sort the eigenvalues to find the two highest values and use
% the corresponding eigenvectors for the LDA plot
[V,Lambda] = eig(Mdl.BetweenSigma, Mdl.Sigma, 'qz'); %Computes eigenvectors in matrix V and eigenvalues in matrix Lambda
[Lambda, sorted] = sort(diag(Lambda), 'descend');


%Scatter plot the eigenvectors with the hightest corresponding eigenvalues
figure (7)
scatter(V(:, sorted(1)), V(:, sorted(2)))

for iterator5 = 1:neuron_data_size(2) - 1
    text(V(iterator5, sorted(1)), V(iterator5, sorted(2)), num2str(iterator5)) 
end

title('Linear Discriminant Analysis Plot (Filtered Data)')
xlabel('Linear Coordiate 1')
ylabel('Linear Coordiate 2')
axis([-1 1 -1 1])

%%
% Plot the LDA by groups relative to Rearing and MB correlation factors
correlation_groups = zeros(neuron_data_size(2) - 1, 5);
for iterator6 = 1:neuron_data_size(2) - 1
   if (Rearing_filtered_data_correlation_factors(iterator6) > 0)
       correlation_groups(iterator6, 1) = correlation_groups(iterator6) + 20;
       
       correlation_groups(iterator6, 2:3) = correlation_groups(iterator6, 2:3)...
           + 100*Rearing_filtered_data_correlation_factors(iterator6);
   end
   if (Rearing_filtered_data_correlation_factors(iterator6) < 0)
       correlation_groups(iterator6, 1) = correlation_groups(iterator6) + 1;
       
       correlation_groups(iterator6, 4:5) = correlation_groups(iterator6, 4:5)...
           - 100*Rearing_filtered_data_correlation_factors(iterator6);
   end
   if (MB_filtered_data_correlation_factors(iterator6) > 0)
       correlation_groups(iterator6, 1) = correlation_groups(iterator6) + 5;
       
      correlation_groups(iterator6, 3) = correlation_groups(iterator6, 3)...
           + 100*MB_filtered_data_correlation_factors(iterator6);
      correlation_groups(iterator6, 5) = correlation_groups(iterator6, 5)...
           + 100*MB_filtered_data_correlation_factors(iterator6);
   end
   if (MB_filtered_data_correlation_factors(iterator6) < 0)
       correlation_groups(iterator6, 1) = correlation_groups(iterator6) + 10;
       
       correlation_groups(iterator6, 2) = correlation_groups(iterator6, 2)...
           - 100*MB_filtered_data_correlation_factors(iterator6);
       correlation_groups(iterator6, 4) = correlation_groups(iterator6, 4)...
           - 100*MB_filtered_data_correlation_factors(iterator6);
   end
end

N_R__P_MB = find(correlation_groups(:, 1) == 6);
N_R__N_MB = find(correlation_groups(:, 1) == 11);
P_R__P_MB = find(correlation_groups(:, 1) == 25);
P_R__N_MB = find(correlation_groups(:, 1) == 30);

figure (8)
% gscatter(V(:, sorted(1)), V(:, sorted(2)), correlation_groups)

scatter(V(N_R__P_MB, sorted(1)), V(N_R__P_MB, sorted(2)), correlation_groups(N_R__P_MB, 5), 'filled')
hold on
scatter(V(N_R__N_MB, sorted(1)), V(N_R__N_MB, sorted(2)), correlation_groups(N_R__N_MB, 4), 'filled')
scatter(V(P_R__P_MB, sorted(1)), V(P_R__P_MB, sorted(2)), correlation_groups(P_R__P_MB, 3), 'filled')
scatter(V(P_R__N_MB, sorted(1)), V(P_R__N_MB, sorted(2)), correlation_groups(P_R__N_MB, 2), 'filled')


title('LDA with Correlation Factors (Filtered Data)')
xlabel('Linear Coordiate 1')
ylabel('Linear Coordiate 2')
legend('Rearing < 0, MB < 0', 'Rearing < 0, MB > 0', 'Rearing > 0, MB > 0', 'Rearing > 0, MB < 0')
axis([-1 1 -1 1])
% for iterator5 = 1:neuron_data_size(2) - 1
%     text(V(iterator5, sorted(1)), V(iterator5, sorted(2)), num2str(iterator5)) 
% end


%%


% For analyzing one specific time period on one specific channel
figure(9)
plot((Time - Time(1)), neuron_data(:, 2)/max(neuron_data(:, 2)))
hold on
plot((Time - Time(1)), filtered_data(:, 1)/max(filtered_data(:, 1)), 'Linewidth', 2)

axis([0 330 -0.5 1])

% Fill in the areas on the graph (Green --> MB, Red --> Rearing)
area((Time - Time(1)), MB_data(:, 1)*neuron_data_size(2),'FaceColor','g', 'FaceAlpha',.2,'EdgeAlpha', 0.001)
area((Time - Time(1)), Rearing_data(:, 1)*neuron_data_size(2),'FaceColor','r', 'FaceAlpha',.2,'EdgeAlpha', 0.001)

xlabel('Time (sec)')
ylabel('Calcium Imaging Amplitude')
legend('Raw Data', 'Filtered Data')

%%






% t = 0:1/20:10;
% X = sin(2.*pi.*t);
% 
% figure (8)
% plot(t, X)

% Y1 = fft(neuron_data);
% Y2 = fft(filtered_data_denoised);
% 
% fs = 20;
% n = 12034;          % number of samples
% 
% f = (0:n-1)*(fs/n);     % frequency range
% 
% power1 = abs(Y1).^2/n;    % power of the DFT
% power2 = abs(Y2).^2/n;    % power of the DFT
% 
% figure(8)
% semilogx(f(1:end/2),pow2db(power1(1:end/2, 2:end)))
% hold on
% semilogx(f(1:end/2),pow2db(power2(1:end/2, :)))
% xlabel('Frequency')
% ylabel('Power')
% legend('Raw Data', 'Filtered Data')
%axis([-1 10 -0.4 4e8])





% %%
% X = [V(:, sorted(1)), V(:, sorted(2))];
% 
% k = 3;
% Sigma = {'diagonal','full'};
% nSigma = numel(Sigma);
% SharedCovariance = {true,false};
% SCtext = {'true','false'};
% nSC = numel(SharedCovariance);
% d = 500;
% x1 = linspace(min(X(:,1)) - 2,max(X(:,1)) + 2,d);
% x2 = linspace(min(X(:,2)) - 2,max(X(:,2)) + 2,d);
% [x1grid,x2grid] = meshgrid(x1,x2);
% X0 = [x1grid(:) x2grid(:)];
% threshold = sqrt(chi2inv(0.99,2));
% options = statset('MaxIter',1000); % Increase number of EM iterations
% 
% figure;
% c = 1;
% for i = 1:nSigma
%     for j = 1:nSC
%         gmfit = fitgmdist(X,k,'CovarianceType',Sigma{i},...
%             'SharedCovariance',SharedCovariance{j},'Options',options);
%         clusterX = cluster(gmfit,X);
%         mahalDist = mahal(gmfit,X0);
%         subplot(2,2,c);
%         h1 = gscatter(X(:,1),X(:,2),clusterX);
%         hold on;
%             for m = 1:k
%                 idx = mahalDist(:,m)<=threshold;
%                 Color = h1(m).Color*0.75 + -0.5*(h1(m).Color - 1);
%                 h2 = plot(X0(idx,1),X0(idx,2),'.','Color',Color,'MarkerSize',1);
%                 uistack(h2,'bottom');
%             end
%         plot(gmfit.mu(:,1),gmfit.mu(:,2),'kx','LineWidth',2,'MarkerSize',10)
%         title(sprintf('Sigma is %s, SharedCovariance = %s',...
%             Sigma{i},SCtext{j}),'FontSize',8)
%         legend(h1,{'1','2','3'});
%         hold off
%         c = c + 1;
%     end
% end