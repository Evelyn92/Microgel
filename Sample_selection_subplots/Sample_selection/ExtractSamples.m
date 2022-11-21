function ExtractSamples(class_name, cloud_add, save_add, max_samples)

    class_num = length(class_name);
    for id_mg = 1:class_num
    mkdir (save_add{id_mg}+"/tr/")
    addpath(save_add{id_mg}+"/tr/")
    cloud = read_point_cloud(cloud_add{id_mg});    
    %2) detect individual microgels:
    [centers_x, centers_y, x_range, y_range] = detect_gels(cloud, 0);

    %3) extract the first of these microgels as a point cloud:

        for idx=1 : min(length(centers_x), max_samples) % we select 5 sample for training and 2 samples for test
            gel1 = get_bounding_box(cloud, centers_x, centers_y, x_range, y_range, idx, 0, 0); 
            x = gel1.XLimits(2)-gel1.XLimits(1);
            y = gel1.YLimits(2)-gel1.YLimits(1);
            z = gel1.ZLimits(2)-gel1.ZLimits(1);
        %     disp(gel1.Count)
            add = strcat(save_add{id_mg}, 'tr/', int2str(idx),'.ply');
            pcwrite(gel1,add,'PLYFormat','binary');
            %pcshow(gel1)
        end

    end

end