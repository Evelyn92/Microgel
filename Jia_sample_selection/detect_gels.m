% cloud = read_point_cloud('LHA_FOCTS_NIPMAM_ALL.3d');
% [centers_x, centers_y, x_range, y_range] = detect_gels(cloud);
% gel1 = get_bounding_box(cloud, centers_x, centers_y, x_range, y_range, 1); 
%
function [centers_x, centers_y, x_range, y_range] = detect_gels(cloud, Flag)

grid_size = 500;

[bandwidth,density,X,Y]= kde2d(cloud.Location(:,1:2), grid_size);

density             = min_max_normalise(density);
density_filtered    = medfilt2(density);
density_thresholded = density_filtered>graythresh(density_filtered);

clusters = bwlabel(density_thresholded);
numbers  = unique(clusters);

centers_x = []; centers_y = [];
x_range   = []; y_range   = [];

for(i=2:length(numbers))
    [row,col] = find(clusters==numbers(i));
    
    %for PEG/ODS
    temp=row; row=col; col=temp;
    %
   
    centers_x = cat(1, centers_x, X(round(median(row)), round(median(row))) ); 
    centers_y = cat(1, centers_y, Y(round(median(col)), round(median(col))) );
    x_range   = cat(2, x_range, [X(min(row), min(row)); X(max(row), max(row))] );
    y_range   = cat(2, y_range, [Y(min(col), min(col)); Y(max(col), max(col))] );
end
if Flag
figure, imagesc(flipud(density)); 
figure, imagesc(flipud(clusters));
end
%for PEG/ODS
temp=centers_x; centers_x = centers_y; centers_y=temp;
temp=x_range; x_range=y_range; y_range=temp;
%

