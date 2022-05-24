% 6/20/18
% Joey Costello

% implementation of LDA based on http://sebastianraschka.com/Articles/2014_python_lda.html
% (could be vectorized for faster computation)

% - X is NxFeatures matrix
% - Y is a Nx1 vector of class assignments (1...numClasses)

% returns a Nx2 matrix of the points transformed to the LDA space



function [Xnew,W] = jcLDA(X,Y, plotOn,markerSize)

numFeatures = size(X,2);
numClasses = max(Y);


%% Step 1: compute mean vectors
m = [];
for i = 1:numClasses
    idx = (Y==i);
    m(:,i) = mean(X(idx,:));
end


%% Step 2: compute scatter matrices

%  within-class matrix
Sw = zeros(numFeatures,numFeatures);
for cl = 1:numClasses
    clScatterMat = zeros(numFeatures,numFeatures);
    
    % get values for this class
    idx = (Y==cl);
    x = X(idx,:); % Nx4 matrix
    
    % iterate over each example
    for r = 1:size(x,1)
        row = x(r,:)'; 
        clScatterMat = clScatterMat + (row-m(:,cl))*(row-m(:,cl))';
    end
    
    % add to the full Sw matrix
    Sw = Sw + clScatterMat;
end



%% between-class matrix

overallMean = mean(X,1)';
Sb = zeros(numFeatures,numFeatures);

for i = 1:size(m,2)
    meanVec = m(:,i);
    % get values for this class
    idx = (Y==i);
    x = X(idx,:);
    % size of this class
    n = size(x,1);
    
    Sb = Sb +  n * (meanVec-overallMean)*(meanVec-overallMean)';
end



%% compute eigenvectors

[V,D] = eig(inv(Sw)*Sb);

% sort so most important vectors are first
[D,I] = sort(diag(D),'descend');
V = V(:, I);

% select the first two eigenvectors
W = -V(:, 1:3);



%% Transform to new space, plot

Xnew = X*W;
C = [[0 0 255] [255 255 0] [255 0 0]];

%**************** 2D Plotting *******************
if (plotOn)
    figure, hold on;
    for i = 1:max(Y)
        markerSize = markerSize;
        scatter(Xnew(Y==i,1),Xnew(Y==i,2),markerSize,'filled','MarkerFaceAlpha', 0);
        
    end
    xlabel('LDA Feature 1');
    ylabel('LDA Feature 2');
end
legend({'Rearing','Marble','Baseline'});




%**************** 3D Plotting *******************
% if (plotOn)
%     figure, hold on;
%     for i = 1:max(Y)
%         scatter3(Xnew(Y==i,1),Xnew(Y==i,2),Xnew(Y==i,3),markerSize,'filled');%,...
%             %'MarkerFaceAlpha',.3);
%             %'MarkerEdgeColor','w',...
%             %'MarkerEdgeAlpha',.3); 
%     end
%     xlabel('LDA Feature 1');
%     ylabel('LDA Feature 2');
%     legend({'Rearing','Marble','Hindleg','Baseline'});
%     axis vis3d;
%     rotate3d on;
% end



% % boundary
% c = {'red','blue','green','red'};
% P = Xnew(Y==i,:);
% k = boundary(P);
% hold on
% trisurf(k,P(:,1),P(:,2),P(:,3),'Facecolor',c{i},'FaceAlpha',0.3)





end









