

%% load data
load('calciumData.mat'); % 'calciumData'
load('eventTimes_3events.mat'); % 'eventTimes' and 'eventTimesAll'

times = calciumData(:,1);
% numEvents = size(eventTimes,2)/2;
numEvents = 2;


%% Setup input parameters

% tmax = 10800;
% cells = 1:37;
tmin = 1;
tmax = length(calciumData);

goodCells = find(~sum(isnan(calciumData(tmin:tmax,(1:50)+1))));
cells = goodCells;

inputData = calciumData(tmin:tmax,cells+1); % ***need to +1 (time first column)
inputTimes = times(tmin:tmax);



%% filter signal - optional
% deconvolve and smooth
% for i = 1:length(cells)
%     [~,inputData(:,i),~] = deconvolveCa(inputData(:,i));
%     inputData(:,i) = smooth(inputData(:,i),60);
% end
% 
% for i = 1:length(cells)
%     [c,s,~] = deconvolveCa(inputData(:,i));
%     %inputData(:,i) = movsum(s>1,10);
%     inputData(:,i) = s;
% end

% %% downsample
% n = 10;
% s1 = size(inputTimes, 1);
% M  = s1 - mod(s1, n);
% inputData2 = zeros(M/n,length(cells));
% for i = 1:length(cells)
%     x = inputData(:,i);
%     y  = reshape(x(1:M), n, []);
%     result = transpose(sum(y, 1) / n);
%     inputData2(:,i) = result;
% end
% inputTimes2 = inputTimes(1:n:end);
% 
% inputTimes = inputTimes2;
% inputData = inputData2;


%% RUN LDA - find index of each event, use as input

dummy = zeros(size(inputData,1),1);
fullIdx = zeros(1,size(inputData,1));
eventMats = {};
for eventNum = 1:numEvents
    % find the indices of each event, store in fullIndex
    eventMatrixIdx = (eventNum-1)*2 + 1; % 1, 3, 5...
    eventMatrixIdx = [eventMatrixIdx eventMatrixIdx+1]; % [1 2], [3 4]...
    [~,eventTimeIdx] = extractDataForEvent2(dummy,inputTimes,eventTimes(:,eventMatrixIdx));
    
    % mark this event in the main index
    fullIdx(logical(eventTimeIdx)) = eventNum;
end

% format lda inputs
X = inputData;
Y = fullIdx;

% make baseline non-zero
Y(Y==0) = numEvents+1;

% run and plot LDA
[Xnew,W] = jcLDA(X,Y,true,3);
Xnew = real(Xnew);


% Variables we have now:
% - fullIdx    (mx1 index of all events 1,2,3..., baseline is 0)
% - Xnew       (mx3 matrix of the lda coordinates




%% MAKE COLOR PLOT OF EACH EVENT



% Create an index of time within each event

% set all values to 1
timeIdx = double(logical(fullIdx));
% loop over, counting up for each event
for i = 2:length(timeIdx)
    if timeIdx(i)>0
        timeIdx(i) = timeIdx(i-1)+1;
    end
end

% now normalize so each event ranges 0-100
timeIdx_norm = timeIdx;
for i = 2:length(timeIdx_norm)
    if timeIdx_norm(i) > 0
        
        % check if on end
        if (i+1) <= length(timeIdx_norm)
            % not end, check if further number
            if timeIdx_norm(i+1) > timeIdx_norm(i)
                continue;
            end
        end
        
        % on the end of an event -> normalize
        thisNum = timeIdx_norm(i);
        timeIdx_norm((i-thisNum+1):i) = round(linspace(1,100,thisNum));
        
    end
end


% Variables we have now:
% - fullIdx      (use to select an event)
% - Xnew         (coordinate values)
% - *timeIdx      (counting index of each event occurance)
% - *timeIdx_norm (normalized time index of each event occurance)




%%%%%%%%%%%%%%%%%%%%
% % create a matrix of colors based on the time within each event

% matrix to store all colors
colors = zeros(length(fullIdx),3);

ncolors = 100;
cmap = jet(ncolors);
temp = timeIdx_norm(timeIdx_norm~=0);
colors(fullIdx~=0,:) = cmap(temp,:);


% Variables we have now:
% - fullIdx      (use to select an event)
% - Xnew         (coordinate values)
% - timeIdx_norm (normalized time index of each event occurance)
% - *colors      (colors based on the normalized time, blue->red)




% %%%%%%%%%%%%%%%%%%%%
% % % Make a simpler color matrix - start,middle,end
% r = [1 0 0];
% g = [0 1 0];
% b = [0 0 1];
% or = [1 0.6 0];
% 
% r = [247 71 49]/255;
% g = [20 150 77]/255;
% b = [21 147 255]/255;
% or = [245 169 49]/255;
% 
% colorsSimple = zeros(length(fullIdx),3);
% eventPartIdx = zeros(1,length(fullIdx));
% 
% % set how long the start/end are (seconds)
% startTime = 0.3;
% endDuration = round(startTime/0.05);
% 
% % iterate over each timept, stop when reach end of an event
% for i = 2:length(timeIdx)
%     if timeIdx(i) > 0
%         
%         % check if on end
%         if (i+1) <= length(timeIdx)
%             % not end, check if further number
%             if timeIdx(i+1) > timeIdx(i)
%                 continue;
%             end
%         end
%         
%         % on the end of an event -> setup event colors
%         thisNum = timeIdx(i);
%         idx = (i-thisNum+1):i;
%         firstHalf = idx(1) : idx(1)+round((idx(end)-idx(1)-0.5)/2);
%         secondHalf = (firstHalf(end)+1) : idx(end);
%         
%         % fill in colors/mark in index
%         colorsSimple(firstHalf,:) = repmat(b,length(firstHalf),1); % blue
%         colorsSimple(secondHalf,:) = repmat(r,length(secondHalf),1); % red
%         eventPartIdx(firstHalf) = 1;
%         eventPartIdx(secondHalf) = 4;
%         if (length(firstHalf)>endDuration)
%             % green
%             colorsSimple(firstHalf((endDuration+1):end),:) = repmat(g,length(firstHalf((endDuration+1):end)),1);
%             eventPartIdx(firstHalf((endDuration+1):end)) = 2;
%         end
%         if (length(secondHalf)>endDuration)
%             % orange
%             colorsSimple(secondHalf(1:(end-endDuration)),:) = repmat(or,length(secondHalf(1:(end-endDuration))),1);
%             firstHalf(secondHalf(1:(end-endDuration))) = 3;
%         end
%         
%     end
% end
% 
% 
% 



%%%%%%%%%%%%%%%%%%
% % Make plots of paths of each event over time
% - can use 'colors' or 'colorsSimple'

%figure, hold on; axis vis3d; rotate3d on;

for event = 1:numEvents
    % create a plot with all points
    jcLDA(X,Y,true,4);
    
    % plot this event points
%     figure;
    pts = Xnew(fullIdx==event,:);
%     scatter3(pts(:,1),pts(:,2),pts(:,3), 20,colors(fullIdx==event,:),'filled','MarkerFaceAlpha',1);
    scatter(pts(:,1),pts(:,2), 10,colors(fullIdx==event,:),'filled','MarkerFaceAlpha',1);
%     b = [143/255 170/255 220/255];
%     db = [100/255 140/255 200/255];
%     r = [220/255 96/255 42/255];
%     dr = [190/255 80/255 32/255];
%     scatter(pts(:,1),pts(:,2), 20, 'MarkerEdgeColor', dr,'MarkerFaceColor', r,'MarkerFaceAlpha',1);
%     title(['Event ' num2str(event)]);
%     % axis lim fit to LDA result - WY
%     ylim([-2000 2000])
%     xlim([-3000 2000])
    
end





% -------------------------------------------------------------------------
%% Calculate mahalanobis distance to baseline over time


for event = 1:numEvents
    
    % get the distances and times of this event
    distances = mahal(Xnew(fullIdx==event,:), Xnew(fullIdx==0,:))';
    thisTimes = timeIdx_norm(fullIdx==event);
    
    mahalDataBins = zeros(1,100); % to sum up data bins
    mahalDataNums = zeros(1,100); % num of vals summed in each bin
    
    % sort into bins and count
    for n = 1:100
        mahalDataBins(n) = mahalDataBins(n) + sum(distances(thisTimes==n));
        mahalDataNums(n) = mahalDataNums(n) + sum(thisTimes==n);
    end
    
    % plot
    figure, hold on;
    smoothing = 5;
    y = smooth(mahalDataBins ./ mahalDataNums,smoothing);
    plot(1:100,y);
    title(['Mahalanobis Distance to Baseline, Event ', num2str(event)]);
    xlabel('Time, event length normalized to 100');
    ylabel('average mahalanobis distance to Baseline');
    
end









%% calculate Discriminability Index (mahalanobis)

disp(' ')
disp(' ')
disp('Distances for normal input data: ');
for event = 1:numEvents
    % find clusters
    baseMat = inputData(fullIdx==0,:);
    eventMat = inputData(fullIdx==event,:);
    
    % calculate mahal distances
    M1 = mahal(eventMat,baseMat);
    M2 = mahal(baseMat,eventMat);
    avgDist = mean([M1; M2]);
    
    disp(['Avg dist event ' num2str(event) ' to baseline is ' num2str(avgDist,3)]);

end

disp(' ')
disp(' ')
disp('Distances for LDA transformed data: ');
for event = 1:numEvents
    % find clusters
    baseMat = Xnew(fullIdx==0,:);
    eventMat = Xnew(fullIdx==event,:);
    
    % calculate mahal distances
    M1 = mahal(eventMat,baseMat);
    M2 = mahal(baseMat,eventMat);
    avgDist = mean([M1; M2]);
    
    disp(['Avg dist event ' num2str(event) ' to baseline is ' num2str(avgDist,3)]);

end















%%




%% LDA (OLD) - events separated into different matrices, then recombined
%
% dummy = zeros(size(inputData,1),1);
% eventMats = {};
% for eventNum = 1:numEvents
%     % extract event/baseline data
%     eventMatrixIdx = (eventNum-1)*2 + 1; % 1, 3, 5...
%     eventMatrixIdx = [eventMatrixIdx eventMatrixIdx+1]; % [1 2], [3 4]...
%     [~,eventTimeIdx] = extractDataForEvent2(dummy,inputTimes,eventTimes(:,eventMatrixIdx));
%
%     % form matrix: each row is a timepoint containing values from all cells
%     eventMats{eventNum} = inputData(logical(eventTimeIdx'),:);
% end
%
% % get baseline
% [~,baseIdx] = getBaseline(dummy,inputTimes,eventTimesAll);
% eventMats{eventNum+1} = inputData(logical(baseIdx'),:);
%
% % format LDA inputs
% X = []; Y = [];
% for i = 1:length(eventMats)
%     X = [X; eventMats{i}];
%     Y = [Y; i*ones(size(eventMats{i},1),1)];
% end
%
% % run and plot LDA
% [Xnew,W] = jcLDA(X,Y,true);
%
%
%
%
% %% highlight hindlimb
%
% rearing = eventMats{1};
% hind = eventMats{3};
%
%  [C,idxRearing,~] = intersect(rearing,hind,'rows');
%
%
% %% set color order
%
% colorOrder = get(groot,'defaultAxesColorOrder');
%
% % swap yellow for green
% temp = colorOrder(3,:);
% colorOrder(3,:) = colorOrder(5,:);
% colorOrder(5,:) = temp;
%
% set(groot,'defaultAxesColorOrder',colorOrder);
%
% co = colorOrder(randperm(size(colorOrder,1)),:);
% set(groot,'defaultAxesColorOrder',co);
% 

