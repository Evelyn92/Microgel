clear all;
clc;
cd D:\forClone\Hiwi\Microgel\Jia_sample_selection
%1) read the entire slide: (the three .3d files contain the data from the Nano Letters paper)
cloud_add{1} = 'LHA_FOCTS_NIPMAM_ALL.3d';
cloud_add{2} = 'LHA_ODS_NIPMAM_ALL.3d';
cloud_add{3} = 'LHA_PEG_NIPMAM_ALL.3d';

save_add{1} = './dataset/focts/'; 
save_add{2} = './dataset/ods/';
save_add{3} = './dataset/peg/';


for id_mg = 1:3
cloud = read_point_cloud(cloud_add{id_mg});    
%2) detect individual microgels:
[centers_x, centers_y, x_range, y_range] = detect_gels(cloud, 0);

%3) extract the first of these microgels as a point cloud:

for idx=1 : min(length(centers_x), 20) % we select 5 sample for training and 2 samples for test
    gel1 = get_bounding_box(cloud, centers_x, centers_y, x_range, y_range, idx, 0, 0); 
    x = gel1.XLimits(2)-gel1.XLimits(1);
    y = gel1.YLimits(2)-gel1.YLimits(1);
    z = gel1.ZLimits(2)-gel1.ZLimits(1);
%     disp(gel1.Count)
    add = strcat(save_add{id_mg}, 'tr/', int2str(idx),'.ply');
%     pcwrite(gel1,add,'PLYFormat','binary');
    pcshow(gel1)
end

end






