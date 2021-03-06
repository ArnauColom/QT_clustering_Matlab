clear;

%% Variables used to customize the behaviour of QT Clustering
showData  = true; % Used to show the initial data
lowRes    = true; % Boolean used to determine whether to use low resolution data or not
step_data = 5;    % value used to undersample the data to be clustered (reduce resolution)
nCenters  = Inf;  % Number of potential centers considered at each iteration over the data to be distributed
threshold = 0.1;  % Threshold distance within the cluster
threshold_dot = cos(pi/10);   % Threshold distance within the cluster

%% Data (for now, only the 3D positions will be used)
if(lowRes)
    Cornell_Arnau_positions_low_res;
else
    Cornell_Arnau_positions;
end

% Normalize the 3D position data
min_x = min(min(data(:,:,1)));
min_y = min(min(data(:,:,2)));
min_z = min(min(data(:,:,3)));
data(:,:,1) = data(:,:,1) - min_x;
data(:,:,2) = data(:,:,2) - min_y;
data(:,:,3) = data(:,:,3) - min_z; 

max_x = max(max(data(:,:,1)));
max_y = max(max(data(:,:,2)));
max_z = max(max(data(:,:,3)));
data(:,:,1) = data(:,:,1) / max_x;
data(:,:,2) = data(:,:,2) / max_y;
data(:,:,3) = data(:,:,3) / max_z; 

% Undersample according to step and reshape
dataUR = reshape( data(1:step_data:end, 1:step_data:end, :), [], 3 ); 

if showData
    figure;
    scatter3(dataUR(:,1), dataUR(:,2), dataUR(:,3), 'fill');
    axis equal;
end

%% QT-Cluster the data
tic

% Initialize the cluster index at zero
clusterIndex = zeros(size(dataUR, 1), 1);

% Initialize the candidate cluster
candidateCluster = [];
c_index = 1; % Index of the current cluster

% Initialize the structure holding the points yet to be assigned to a
% cluster
dataToCluster      = dataUR;
dataToClusterIndex = 1:1:size(dataToCluster, 1); % Keep the index

while ~isempty(dataToCluster)
    fprintf('\nAssigning points to cluster number %d \n', c_index);
    fprintf('There are %d points to assign\n', size(dataToCluster, 1));
    
    % For each datapoint to be clustered
    step = max(1, round(size(dataToCluster, 1) / nCenters));
    for i = 1:step:size(dataToCluster, 1)
        p_i = dataToCluster(i,:);
        % Build a candidate cluster for the datapoint
        potentialCandidateCluster = [];
        for j = 1:size(dataToCluster,1)
            c_j = dataToCluster(j,:);
            if dist(p_i, c_j) < threshold
                potentialCandidateCluster(end+1) = j;
            end
        end
        
        % If the current candidate cluster has more elements than those in
        % candidate cluster, then keep the current candidate cluster
        if length(potentialCandidateCluster) > length(candidateCluster)
            candidateCluster = potentialCandidateCluster;
        end
    end
    
    fprintf('%d points have been selected to form the cluster %d \n', length(candidateCluster), c_index);
    
    % Save the selected candidate cluster
    clusterIndex(dataToClusterIndex(candidateCluster)) = c_index;
    c_index = c_index + 1;
    
    % Take the point already clustered out of the data
    dataToCluster(candidateCluster, :) = [];
    dataToClusterIndex(candidateCluster) = [];
    
    candidateCluster = [];
end
toc

if showData
    figure;
    for i=1:c_index-1
        indexes = find(clusterIndex == i);
        markerSize = ones(size(indexes)) * 20;
        rnd_color = ones(size(indexes))*rand(1,3);
        scatter3(dataUR(indexes,1), dataUR(indexes,2), dataUR(indexes,3), markerSize, rnd_color, 'fill');
        hold on
        axis equal;
    end
end

save('test.mat');

function d = dist(p1, p2)
    d = sqrt(sum((p1-p2).^2));
end

% function d = dist(p1, p2, n1, n2)
%     if (dot(n1, n2) > threshold_dot)
%         d = sqrt(sum((p1-p2).^2));
%     else
%         d = Inf;
%     end
% end
