function [descend] = plotData3Threshold(behavior, start_stop, Alldata, Time, neuronNum,threshold)
%plot neural colormap raster plot at behavior epoch time
% 2018/11/27 Wooyeon Shin
twosecond = find(Time<= Time(start_stop(1,1),1)+1 , 1, 'last') - find(Time>= Time(start_stop(1,1),1)-1 , 1, 'first')+1;

data = zeros(twosecond, length(Alldata(1,:)));

% cut time -1 to -1 second from the epochs
for i= 1:length(start_stop)
    if start_stop(i,1) ~= 0
        start = find(Time>= Time(start_stop(i,1),1)-1 , 1, 'first');
        stop = find(Time<= Time(start_stop(i,1),1)+1 , 1, 'last');

        if (stop-start+1) == twosecond
            data = data+Alldata(start:stop, :);
        end
    end
end

time = Time(start:stop);


% min max for normalization
Min = min(min(data));
Max = max(max(data));
epoch = round(length(time)/2)-1;

% set data by threshold
if threshold == 1
    data = (data>5000).*data;
    Min = 0;
    Min = Min+(Max-Min)*0/10;
    Max = Min+(Max-Min)*7/10;
elseif threshold == 2
    data = (data>100).*data;
    data = (data<=5000).*data;
    Min = 0;
    Min = Min+(Max-Min)*0/10;
    Max = Min+(Max-Min)*7/10;
elseif threshold == 3
    data = (data<0).*data;
    Max = 0;  
    Min = Min+(Max-Min)*6/10;
    Max = Min+(Max-Min)*10/10;
end

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
descend = sortrows(afterEpochSig,2,'descend');


for i = 1:neuronNum
    j = descend(i,1);
    sortedRaster(:,i) = raster(:,j);
end

figure('Position',[1,1,1400,450])
subplot(1,2,1)
imagesc(raster')
colormap('jet')
% colormap('default')
caxis([Min Max])
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
caxis([Min Max])

colorbar('southoutside', 'Ticks', [])
hold on
line([epoch epoch], [0 neuronNum], 'Color', 'w', 'LineStyle', '--',  'LineWidth', 3)
xticks({})
title({behavior,'start'})

end
