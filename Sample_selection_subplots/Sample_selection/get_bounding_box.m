% gel1 = get_bounding_box(cloud, centers_x, centers_y, x_range, y_range, 1);
%
function bounding_box = get_bounding_box(cloud, centers_x, centers_y, x_range, y_range, gel_nr, flag_show, flag_photon)

factor = 1.2; 
dist_threshold = 700;

x_lower = find(cloud.Location(:,2) >= x_range(1, gel_nr) * 1/factor ); 
x_upper = find(cloud.Location(:,2) <= x_range(2, gel_nr) * factor );
y_lower = find(cloud.Location(:,1) >= y_range(1, gel_nr) * 1/factor ); 
y_upper = find(cloud.Location(:,1) <= y_range(2, gel_nr) * factor );

y = intersect(y_upper, y_lower);
x = intersect(x_upper, x_lower);
points = intersect(y,x);

%remove points that are closer to another center/cloud (can occur if two clouds are close neighbours):
remove_list=[];
for(i=1:length(points))
    d = zeros(1,length(centers_x));
    for(j=1:length(centers_x))
        d(j) = sqrt(sum((cloud.Location(points(i),[2,1]) - [centers_x(j) centers_y(j)]).^2));
    end
    [minimum,min_pos] = min(d);
    if(min_pos~=gel_nr || d(gel_nr)>dist_threshold) 
        remove_list=cat(1,remove_list,i);
    end
end
points(remove_list)=[];


bounding_box = [cloud.Location(points,1), cloud.Location(points,2), cloud.Location(points,3)];

bounding_box = pointCloud(bounding_box);
if flag_photon
    photon  = cloud.Intensity(points,1);
    bounding_box.Intensity = photon;
end
if flag_show
f=figure();
scatter3(cloud.Location(:,1), cloud.Location(:,2), cloud.Location(:,3), ones(1,size(cloud.Location,1)), 'black', 'filled'); 
hold on;
scatter3(cloud.Location(points,1),cloud.Location(points,2),cloud.Location(points,3), ones(1,length(points))*3,'red', 'filled');
daspect([1 1 1]);
end


