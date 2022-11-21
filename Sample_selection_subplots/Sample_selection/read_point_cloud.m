% cloud = read_point_cloud('LHA_FOCTS_NIPMAM_SINGLE.csv');
% pcshow(cloud)
%
function cloud = read_point_cloud(filename, photon)

data  = dlmread(filename);

if(size(data,2)>3) 
    data = data(:,1:3);
end

data(:,3) = 1-data(:,3); %flip z dimension

cloud = pointCloud(data);
pcshow(cloud); 















