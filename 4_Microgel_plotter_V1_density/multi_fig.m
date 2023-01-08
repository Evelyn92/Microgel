clc;
clear all;

path = "D:\forClone\Hiwi\results_mtemp_hilo";
data_type = '_go';%'_go' or 'gg'
Visible = true;
for tmp_inx = 1:6
    cols=4;
    
if tmp_inx ==1
    tmp_info = '21C';
elseif tmp_inx ==2
    tmp_info = '33C';
elseif tmp_inx ==3
    tmp_info = '35C';
elseif tmp_inx ==4
    tmp_info = '38C';
elseif tmp_inx ==5
    tmp_info = '43C';
else
    tmp_info = '53C';
end
    %https://de.mathworks.com/matlabcentral/answers/163904-how-to-display-multiple-image-in-one-figure-window
    %https://de.mathworks.com/help/images/display-multiple-images.html
    for sample_inx = 1:cols
        img_path = strcat(path,'/2d_',data_type,tmp_info,'_sample',string(sample_inx),'.png');
        %'./processed_samples/'
        img(:,:,:,sample_inx) = imread(img_path); 
        
    end
    
    if Visible
      figure("Visible", 'on');
    else
      figure("Visible", 'off');
    end
    montage(img,'size',[1 NaN]);
    
%     if className ~= "FOCTS"
%     gca.Position = [10 10 1200 1200];
%     end
    fname_2d = strcat(path,'/2d_',data_type,tmp_info,'.png');
    saveas(gcf,fname_2d)
   
end


