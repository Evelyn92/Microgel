function SelectedVsDiscarded(path)
cd (path)
clear all;

class_name = {'FOCTS','ODS','PEG'};
selected_path = 'selected_samples/';

ranking_info = dir('./selected_samples/*.txt');
class_num = length(ranking_info);
Visible = false;
for idx = 1:class_num
    clear img;
    save_add = strcat('./samples/', class_name{idx}, '/');
    save_add_path = save_add+"tr/"; 
    className = class_name{idx};
     
    
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
   
    %https://de.mathworks.com/matlabcentral/answers/163904-how-to-display-multiple-image-in-one-figure-window
    %https://de.mathworks.com/help/images/display-multiple-images.html
    for j = 1:cols
         
        index = rank_inx_arr(j);
        img_path = strcat('./processed_samples/',className,'/',className,'-',string(index),'-density.png');
        img(:,:,:,j) = imread(img_path);
      
        
        index = rank_inx_arr(sam_num-cols+j);
        img_path = strcat('processed_samples/',className,'/',className,'-',string(index),'-density.png');
        img(:,:,:,2*cols-j+1) = imread(img_path);  
        
    end
    
    if Visible
      figure("Visible", 'on');
    else
      figure("Visible", 'off');
    end
    montage(img,'size',[2 NaN]);
    
%     if className ~= "FOCTS"
%     gca.Position = [10 10 1200 1200];
%     end
    title("The upper row: SELECTED | The lower row: DROPPED")
    fname_2d = strcat(selected_path,className,'_Selected_VS_Discarded','-density.png');
    saveas(gcf,fname_2d)
    disp(strcat('The Selected VS Discarded figure is saved in: ',fname_2d));
   
end

end
