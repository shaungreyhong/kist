% 6/21/18
% Joey Costello

% plots a points color based on how far into the event it is

% shows how the ensemble might be dynamic, if certains colors tend to
% cluster together

% timepoints early on are dark blue, after about 3 seconds are red, then
% brown



idx = double(event1Idx);
idx2 = [];
for i = 2:length(idx)
    if idx(i)>0
        idx(i) = idx(i-1)+1;
        idx2 = [idx2 idx(i-1)+1];
    end
end
time1Idx = idx;
time1Idx_short = idx2;


idx = double(event2Idx);
idx2 = [];
for i = 2:length(idx)
    if idx(i)>0
        idx(i) = idx(i-1)+1;
        idx2 = [idx2 idx(i-1)+1];
    end
end
time2Idx = idx;
time2Idx_short = idx2;


idx = double(baseIdx);
idx2 = [1];
for i = 2:length(idx)
    if idx(i)>0
        idx(i) = idx(i-1)+1;
        idx2 = [idx2 idx(i-1)+1];
    end
end
timebaseIdx = idx;
timebaseIdx_short = idx2;



%% create colors based on the time
ncolors = 100;
cmap = jet(ncolors);
time1Idx_short = time1Idx_short+0;
time1Idx_short(time1Idx_short>ncolors) = ncolors;
c1 = cmap(time1Idx_short,:);


ncolors = 80;
cmap = jet(ncolors);
time2Idx_short = time2Idx_short+0;
time2Idx_short(time2Idx_short>ncolors) = ncolors;
c2 = cmap(time2Idx_short,:);


ncolors = 150;
cmap = jet(ncolors);
timebaseIdx_short(timebaseIdx_short>ncolors) = ncolors;
c3 = cmap(timebaseIdx_short,:);

%%
figure, hold on, axis vis3d, rotate3d on;
% *** keep in mind the indexes are not necessarily in order ***
i = 1;
scatter3(Xnew(Y==i,1),Xnew(Y==i,2),Xnew(Y==i,3),30,c1,'filled');

title('event2');





