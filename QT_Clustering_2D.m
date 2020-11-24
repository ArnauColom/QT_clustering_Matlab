clear all
close all

%% QT CLUSTERING SIMPLE


thr = 2;
min_cluster_size = 1;

%Data generation
num_points = 1000;

x=rand(1,num_points)*10;
y=rand(1,num_points)*10;
%z=rand(1,100)*10;
scatter(x,y)

points = [x;y];

selected = zeros(1,num_points);
label = [];

num_points_remaining = num_points;
points_e = points;

points_out = [];

cluster_num = 1;

points_selected = 0;
tic
while(num_points_remaining>0)    
    
    max_points = 0;
    pos_cluster = 0;
    
    pos_clust = [];
    for i=1:num_points_remaining    
        pos_aux = [];
        total_points = 0;
        for j=1:num_points_remaining            
            distance = norm(points(:,i)-points(:,j));
            if(distance<thr)
               total_points = total_points + 1; 
               pos_aux(total_points) = j;
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
       points_out(:,i+points_selected) = points(:,pos_clust(i)-count);
       label(i+points_selected) = cluster_num;
       points(:,pos_clust(i)-count) = [];
       count = count +1;
    end
    
    cluster_num = cluster_num +1;
    num_points_remaining = num_points_remaining - max_points;
    points_selected = points_selected + max_points;

end
toc

PlotClusters(points_out',label+1);



