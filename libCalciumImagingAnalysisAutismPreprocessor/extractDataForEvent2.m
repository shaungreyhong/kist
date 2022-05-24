% 6/15/18
% Joey Costello

% extracts all calcium data between start/end of event
% return a vector values
% ('extractDataForEvent' only returns mean values, not all)

% - times is a 1xm vector of times
% - data is a 1xm vector of the calcium value for each time
% - eventTimes is a nx2 matrix, 1st col start times, 2nd col end times


% - eventVals returns a 1xn vector of the values



function [eventVals,index] = extractDataForEvent2(data,times,eventTimes)

eventTimes(any(isnan(eventTimes), 2), :) = []; % remove nan values (row-wise)
eventVals = [];

% create an index of which time points are events (use all events)
index = zeros(1,length(data));

% loop over each event epoch
for i = 1:size(eventTimes,1)
    tmin = eventTimes(i,1);
    tmax = eventTimes(i,2);
    
    % make sure this is within the times
    if (tmax > max(times))
        continue;
    end
    
    % find closest time point to min
    [~,idxMin] = min(abs(times-tmin));
    
    % find closest time point to max
    [~,idxMax] = min(abs(times-tmax));
    
    % get points between the times
    thisValues = data(idxMin:idxMax)';
    eventVals = [eventVals, thisValues];
    
    % mark these points in the overall index
    index(idxMin:idxMax) = 1;
end

eventVals = eventVals(~isnan(eventVals)); % remove nan values












%% Plotting (optional)
% % plot all events
% time = -(numPts*.05):.05:(numPts*.05);
% figure; hold on;
% for i = 1:size(e1,1)
%     plot(time,e1(i,:),'b');
% end
% % plot average
% figure, 
% plot(time,mean(e1),'r','LineWidth',2);



end