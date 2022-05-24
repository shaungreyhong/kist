function [descend, sortedRaster] = plotFilteredData(behavior, start_stop, Alldata, Time, neuronNum, filtered)
%plot neural colormap raster plot at behavior epoch time
% 2018/11/27 Wooyeon Shin
a =2; b = 10;

twosecond = round(1/(Time(2) - Time(1)))*2+1;

data = zeros(twosecond, length(Alldata(1,:)));

% cut time -1 to -1 second from the epochs
for i= 1:length(start_stop)
    if start_stop(i,1) ~= 0
        start = find(Time>= Time(start_stop(i,1),1)-1.01 , 1, 'first');
        stop = find(Time<= Time(start_stop(i,1),1)+1.01 , 1, 'last');

        if (stop-start+1) == twosecond
            data = data+Alldata(start:stop, :);
        end
    end
end

time = Time(start:stop);

if filtered
    data = data(:, 1:2:end);
end

% min max for normalization
Min = min(min(data));
Max = max(max(data));
epoch = round(length(time)/2)-1;

% figure unshorted firings
raster = zeros(length(time), neuronNum);
sortedRaster = zeros(length(time), neuronNum);
afterEpochSig = zeros(neuronNum, 2);

for i = 1:neuronNum
    raster(:,i) = data(:, i) - mean(data(1:epoch, i));    
%     raster(:,i) = data(:, i);
    %data for sorting(comparing signal before and after epoch)
    afterEpochSig(i,1) = i;
    afterEpochSig(i,2) = sum(data(epoch+1:end,i)) - sum(data(1:epoch, i));
end

% sort data
Min = min(min(raster));
Max = max(max(raster));
descend = sortrows(afterEpochSig,2,'descend');


for i = 1:neuronNum
    j = descend(i,1);
    sortedRaster(:,i) = raster(:,j);
end

figure('Position',[1,1,1400,500])
subplot(1,2,1)
imagesc(raster')
colormap('jet')
% colormap('default')
% caxis([Min Max])
caxis([Min+(Max-Min)*a/10 Min+(Max-Min)*b/10])
% caxis([Min Max/2])

colorbar('southoutside', 'Ticks', [])
hold on
line([epoch epoch], [0 neuronNum-1], 'Color', 'w', 'LineStyle', '--', 'LineWidth', 3)
xticks({})
title({behavior,'start'})



subplot(1,2,2)
imagesc(sortedRaster')
colormap('jet')
% colormap('default')

% caxis([Min Max])
caxis([Min+(Max-Min)*a/10 Min+(Max-Min)*b/10])
% caxis([Min Max/2])

colorbar('southoutside', 'Ticks', [])
hold on
line([epoch epoch], [0 neuronNum], 'Color', 'w', 'LineStyle', '--',  'LineWidth', 3)
xticks({})
title({behavior,'start'})

end
