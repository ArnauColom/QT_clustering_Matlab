%% QT CLUSTERING PARAMETERS
thr = 2;
min_cluster_size = 1;

%% DATA GENERATION

Cornell_Arnau_normals_low_res;
normals = data;

num_points = numel(normals(:,:,1));%Number of points
% Show mitsuba file with normals
%figure
%imshow(normals);

normals = 2 * normals - 1;
normals_array = reshape(normals, [num_points 3]);

low_cent_norm = normals(1:10:end,1:10:end,:);
low_cent_norm_array = reshape(low_cent_norm, [num_points/100 3]);
% Read mitsuba file with positions (normalized between zero and one in each
% dimension)

Cornell_Arnau_positions_low_res;
points = data;
points_array = reshape(points, [num_points 3]);

low_cent_point =  points(1:10:end,1:10:end,:);
low_cent_points_array = reshape(low_cent_point, [num_points/100 3]);
% Reuslts

position_arr = 1:num_points;
position_arr = position_arr';

%% NECESSARY VARIABLES


label = [];
points_out = zeros([num_points 1]);

num_points_remaining = num_points;
points_selected = 0;
cluster_num = 1;


tic
while(num_points_remaining>0)    
    
    max_points = 0;
    pos_cluster = 0;
    
    pos_clust = [];
    
    for i=1:(num_points/100)
        
        pos_aux = [];
        total_points = 0; 
        normal_i = low_cent_norm_array(i,:);
        point_i = low_cent_points_array(i,:);
        
        for j=1:num_points_remaining 
            if(dot(normal_i,normals_array(j,:))>0)
                distance = norm(point_i-points_array(j,:));
                if(distance<thr)
                  total_points = total_points + 1; 
                  pos_aux(total_points) = j;
                end    
            end                    
        end 
        
        
        if(total_points > max_points)
            max_points = total_points;
            pos_cluster = i;
            pos_clust = pos_aux;
        end    
        
    end
    
    
    count = 0;
    for i=1:max_points
        
       points_out(i+points_selected,:) = position_arr(pos_clust(i),:);
       label(i+points_selected) = cluster_num;
       
       %Delete Points
       points_array(pos_clust(i)-count,:) = [];       
       normals_array(pos_clust(i)-count,:) = [];
       position_arr(pos_clust(i)-count,:) = [];
       count = count +1;
    end
    
    cluster_num = cluster_num +1;
    num_points_remaining = num_points_remaining - max_points;
    points_selected = points_selected + max_points;
end
toc

PlotClusters(points_out',label+1);