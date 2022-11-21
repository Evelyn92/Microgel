function SampleSelection(save_add, processed_path, class_name,min_inx, max_inx, selected_sample_num)
class_num = length(class_name);


for idx = 1:class_num
    density_vector = [];
    className = class_name{idx};
    save_add_path = save_add{idx}+"tr/"; 
  
    mkdir (processed_path+"/"+className+'/')
    addpath (processed_path+"/"+className+'/')
    
    selected_path = 'selected_samples/';
    mkdir (selected_path+"/"+className+'/')
    addpath (selected_path+"/"+className+'/')
    sample_files = dir(save_add_path+'*.ply');
    num_sample = length(sample_files);
    max_index = min(max_inx,num_sample);
    if strcmp(className, 'FOCTS')
        
        h_bias = 200;
        d_max = 600;
        h_max = 400;
    elseif strcmp(className, 'ODS')
%         max_index =  min(max_inx,78);
        h_bias = 400;
        d_max = 600;
        h_max = 800;
    else 
%         max_index = min(max_inx,8);
        h_bias = 400;
        d_max = 420;
        h_max = 800;
    end
    plot_config = {h_bias, d_max, h_max, className};

    for index = min_inx:max_index %inside each dataset index_sample
        
       density_vector = VisualizeDensityPlot(density_vector, save_add_path,processed_path, index, plot_config);

    end


    % Start to calculate the correlation scores and select the samples
    close all;
%     index_max=0;
    corr_mat = zeros(max_index,max_index);
    for cor_i=1:max_index
%         sum_cur=0;
        for cor_j=1:max_index
            r = corrcoef(density_vector(:,cor_i), density_vector(:,cor_j));
            corr_mat(cor_i,cor_j)=r(1,2);
        end
    end
    [sorted_corr_sum, sample_index] = sort(sum(corr_mat),'descend');
    % Write the result of sample selection into txt
    Sample_sel_info=strcat(selected_path,className,'-sample_selection_info.txt');
%    % open file identifier
%     fid=fopen(Sample_sel_info,'w');
%     fprintf(fid, ['The sorted samples w.r.t. the sum of correlation'  '\n']);
%     fclose(fid);
    T = table(sample_index', sorted_corr_sum', 'VariableNames', { 'Sample Index', 'Sum of Correlation'} );
    % Write data to text file
    writetable(T, Sample_sel_info)
    disp('The ranked correlation scores are saved in the path:')
    disp(Sample_sel_info)
    
    %Start generate the plots for selected samples
    for inx = 1:selected_sample_num
        % add subplot here for selected ones
        % add subplot here for unselected ones
       index = sample_index(inx);      
       density_vector = VisualizeDensityPlot(density_vector, save_add_path,selected_path, index, plot_config);
    end
    disp('The density plots and point cloud plots are saved in the path:')
    disp(selected_path)
    
    
end




end