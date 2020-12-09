clear;

%% Variables used to customize the behaviour of QT Clustering
showData  = true; % Used to show the initial data
lowRes    = true; % Boolean used to determine whether to use low resolution data or not
step_data = 5;   % value used to undersample the data to be clustered (reduce resolution)
% nCenters  = 20;    % Number of potential centers considered at each iteration over the data to be distributed
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
dataU  = data(1:step_data:end, 1:step_data:end, :);
dataUR = reshape( dataU, [], 3 ); 

if showData
    figure;
    scatter3(dataUR(:,1), dataUR(:,2), dataUR(:,3), 'fill');
    axis equal;
end

%% QT-Cluster the data
tic

% Initialize the structure containing the cluster to which each data point 
%   belongs, with value zero. In this version this is a 2D matrix with the
%   same resolution as the undersampled data (dataU).
clusterIndex = zeros( size(dataU, 1), size(dataU, 2) );
nRows        = size(clusterIndex, 1);
nCols        = size(clusterIndex, 2);
sz = [nRows, nCols];

% Initialize the members of the candidate cluster, initially set to empty.
candidateCluster = [];
c_index          = 1; % Index of the current cluster

% Initialize the structure holding the points yet to be assigned to a
% cluster. In this case, the format is a vector with data points
% (dataToCluster) and another vector with the index (1D) of this point.
% These vectors must keep the same size throughout the algorithm (that is,
% when we eliminate an element from dataToCluster, this same element should
% also be eliminated from dataToClusterIndex). These two vectors can be
% thought of as two members of the same data structure.
%dataToCluster      = dataUR;
%dataToClusterIndex = 1:1:size(dataToCluster, 1); % Keep the index
nPointsToCluster = nRows * nCols;

% Test to verify indexes (this does the same as the previous line).
%rows = repmat(1:nRows, 1, nCols);
%cols = reshape(repmat(1:nCols, nRows, 1), 1, []);
%test = sub2ind([nRows, nCols], rows, cols);

% Search for potential/candidate custers while there are still datapoints
% to assign to clusters.
while (nPointsToCluster > 0)
    fprintf('\nAssigning points to cluster number %d \n', c_index);
    fprintf('There are %d points to assign\n', nPointsToCluster);
    
    % Compute the step used to search for potential cluster centers. This
    % allows to skip points in the dataToCluster vector, meaning that part
    % of those data points will not be used to find the candidate cluster
    % with more datapoints. Setting step < size(dataToCluster, 1) reduces 
    % the quality of the final result and reduces execution time.
%     step = max(1, round(nPointsToCluster / nCenters));
    step = 1;
    % Initialize a candidate cluster for the datapoint
    
    candidateCluster = []; % The largest one found so far
    for i = 1:step:nRows
        for j = 1:step:nCols
            potentialCandidateCluster = []; % The current one (centered in p_i)
            % If the current point (lin, col) hasn't been clustered yet
            if clusterIndex(i, j) <= 0
                % Select the current data point as a potential cluster center
                p_ij = reshape(dataU(i,j,:), 1, 3);
                
                
                % Search (over the neighbours in the image plane) the potential
                % data points for this hypothetical cluster. This
                % is done so that we can make the search for potential cluster
                % elements (centered on p_i) over a neighbourhood in the original
                % image undersampled (dataU). The idea is to avoid having to search
                % over all elements in dataToCluster which is very costly
                % (especially at the beginning when there are many points to be
                % clustered).
                % Search (over the neighbours in the image plane) the potential
                % data points for this hypothetical cluster
                step_local = 7;
                for k = i-step_local:i+step_local
                    if (k>0 && k<=nRows)
                        for l = j-step_local:j+step_local
                            if(l>0 && l<=nCols)
                                % Ensure that this datapoint has NOT already been
                                % clustered
                                if (clusterIndex(k, l) <= 0)
                                    c_kl = reshape(dataU(k,l, :), 1, 3);
                                    if dist(p_ij, c_kl) < threshold
                                        %potentialCandidateCluster(end+1,:) = [l, k];
                                        potentialCandidateCluster(end+1,:) = [k, l];
                                    end
                                end
                            end
                        end
                    end
                end
            end
            % If the current candidate cluster has more elements than those in
            % candidate cluster, then keep the current candidate cluster
            if length(potentialCandidateCluster) > length(candidateCluster)
                candidateCluster = potentialCandidateCluster;
            end
        end
    end
        
    fprintf('%d points have been selected to form the cluster %d \n', length(candidateCluster), c_index);
    
    % Save the selected candidate cluster
    %[row, col] = ind2sub(sz, dataToClusterIndex(candidateCluster)); 
    clusterIndex(sub2ind([nRows, nCols], candidateCluster(:,1), candidateCluster(:,2))) = c_index;
%     clusterIndex(dataToClusterIndex(candidateCluster)) = c_index;
    c_index = c_index + 1;
    
    % Take the point already clustered out of the data
    %dataToCluster(candidateCluster, :) = [];
    %dataToClusterIndex(candidateCluster) = [];
    
    nPointsToCluster = nPointsToCluster - length(candidateCluster);
    candidateCluster = [];
    
end
toc

if showData
    figure;
    for i=1:c_index-1
        %indexes = find(clusterIndex == i);
        [row, col] = find(clusterIndex == i);
        markerSize = ones(size(row)) * 20;
        rnd_color = ones(size(row))*rand(1,3);
        %scatter3(dataUR(indexes,1), dataUR(indexes,2), dataUR(indexes,3), markerSize, rnd_color, 'fill');
        auxX = dataU(:,:,1);
        auxY = dataU(:,:,2);
        auxZ = dataU(:,:,3);
        ind = sub2ind([nRows,nCols],row,col);
        scatter3(auxX(ind), auxY(ind), auxZ(ind), markerSize, rnd_color, 'fill');
        hold on
        axis equal;
    end
end

save('test.mat');

% % Convert from 1D (array) coordinates to 2D (matrix) coordinates
% function [lin, col] = lineToMatrix(x1D, nCols)
%     lin = fix( (x1D-1) / nCols) + 1; % Get the integer part
%     col = rem(x1D-1, nCols) + 1;  % Get the remainder part
% end

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
