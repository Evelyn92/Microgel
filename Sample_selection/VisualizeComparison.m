clc;
clear all;
cd ('D:\forClone\Hiwi\Sample_selection')
class_name = {'FOCTS','ODS','PEG'};
%VisualizeComparison(class_name)

%function VisualizeComparison(class_name)
selected_path = 'selected_samples/';

ranking_info = dir('./selected_samples/*.txt');
class_num = length(ranking_info);

for idx = 1:class_num
    save_add = strcat('./samples/', class_name{idx}, '/');
    save_add_path = save_add+"tr/"; 
    className = class_name{idx};
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
    
    txt_fname = ranking_info(idx).name;
    txt_dir = strcat(selected_path, txt_fname);
    f_id = fopen(txt_dir);
    C_text = textscan(f_id,"%s",2,'Delimiter',',');
    info = textscan(f_id,"%d,%f");
    fclose(f_id);
    
    rank_inx_arr = info{1};
    sam_num = length(rank_inx_arr);
    if sam_num<8
        cols = 3;
    elseif sam_num>8
        cols = 4;  
    end
    
    for k = 1:2*cols
        ax(k) = subplot(2,cols,k);
    end
    
    for j = 1:cols
       % subplot(ax(j))
        index = rank_inx_arr(j);
        density_vector=[];
        VisualizeDensitySubPlot(ax(j), density_vector, save_add_path, selected_path, index, plot_config)
        
        ax2 = subplot(8,2,j);
        density_vector=[];
        index = rank_inx_arr(sam_num-cols+j);
        VisualizeDensitySubPlot(ax2, density_vector, save_add_path, selected_path, index, plot_config)
        
      
    end
    
    fname_2d = strcat(selected_path,className,'/',className,'Selected_VS_Discarded','-density.png');
    saveas(gcf,fname_2d)
   
%end

end