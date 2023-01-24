%save_add  = './dataset/timelife/'; 
clc;
clear all;
userpath('clear')

cloud_add = '16-FM-raw-data.csv';
save_add = {'results_lifetime/'};
Visible = 'on';
% M = readtable(cloud_add);
% data = table2array(M(:,[3,4,5,10]));%10 col: tau_fast
% cloud = pointCloud(data(:,1:3));
% pcshow(cloud); 

class_num = 1;
className = '16-FM';
min_inx = 1;
max_inx = 6;
sample_type = {'go', 'gg'};
for idx = 1:class_num
    density_vector = [];
    %className = class_name{idx};
    save_add_path = string(save_add{1}); 
    sample_files = dir(save_add_path+'*.ply');
    num_sample = length(sample_files);
    
    h_bias = 2;
    
    
    d_max = 3.5;
    h_max = 5;

    clear count V d_range h_range;
    d_dist=0.1;
    h_dist=0.1;

    labelx_dist=1;
    labely_dist=1;

    d_range=0:d_dist:d_max;
    h_range=0:h_dist:h_max;
    
    for index = 1:1
        add = strcat(save_add_path,'sample_',sample_type{idx},string(index),'.ply');
        all_localiz  = dlmread(add);
        add_info = strcat(save_add_path,'\sample_','info',string(index),".ply");
        info  = dlmread(add_info);
        all_localiz  = all_localiz  .* (info(2,:)-info(3,:)) + info(3,:)+ info(1,:);

        zvec=all_localiz(:,3);
        zvec = zvec-median(zvec);
        xcoord=all_localiz(:,1);
        deltax = xcoord-median(xcoord);
        ycoord=all_localiz(:,2);
        deltay = ycoord-median(ycoord);
        %data = all_localiz;
        data = [deltax, deltay, zvec];
%         scatter3(deltax,deltay,zvec, 'filled')
%         Visualization the pd
        cloud_sample = pointCloud(data);
        pcshow(cloud_sample, 'MarkerSize', 500);
        fname_2d = save_add_path+"PC-S"+string(index)+'.png';
%         saveas(gcf,fname_2d)
        
        for i=1:length(d_range)-1
            V(i)=((d_range(i)+d_dist)^2-(d_range(i))^2)*pi*h_dist;
        end
        
        d_axis = sqrt(deltax.^2+deltay.^2);
        h = (zvec-median(zvec))+h_bias;
        data=[d_axis h];
        for i=1:length(d_range)-1
            data((data(:,1)>d_range(i))&(data(:,1)<=d_range(i+1)),3)=i;
        end    
        
         for i=1:length(h_range)-1
            data((data(:,2)>h_range(i))&(data(:,2)<=h_range(i+1)),4)=i;  
         end

        count=zeros(length(d_range)-1,length(h_range)-1);
        data=data(data(:,3)>0,:); % if a data point is out of the x range, throw it away
        data=data(data(:,4)>0,:);% if a data point is out of the y range, throw it away

        for i=1:size(data,1)
            try
                count(data(i,3),data(i,4))=count(data(i,3),data(i,4))+1; 
                %count{(data(i,3),data(i,4))}=data(i,3),data(i,4),data(i,5); 
               %count(data(i,5))=count(data(i,5))+1
            catch
                i;
            end
        end
        
         % visualization
        figure_height_in_pixel=500;
        if Visible
            h=figure("Visible", 'on');
        else
            h=figure("Visible", 'off');
        end
        set(h,'Position',[20 20 figure_height_in_pixel max(h_range)/max(d_range)*figure_height_in_pixel]);
        count=count'./repmat(V,size(count,2),1);
    %     countbin = count;
        density_vector = [density_vector reshape(count,[],1)];
        imagesc(count); % To plot wrt density values
        %imagesc(medint); % To plot wrt I_ratios

        axis equal % suggestion by Eric to prevent elonagation in z
        set(gca,'Ydir','Normal');
        set(gca,'FontSize',22);
        xlim([0.5 max(d_range/d_dist)+1]);
        ylim([0.5 max(h_range/h_dist)+1]);


        set(gca,'Ydir','Normal');
        set(gca,'FontSize',22);
        xlim([0.5 max(d_range/d_dist)+1]);
        ylim([0.5 max(h_range/h_dist)+1]);
        title('\fontsize{30}\fontname{Calibri} 16-FM-'+string(index))
        xlabel('Dist. from symmetry axis','FontSize',22);
        ylabel('Relative z','FontSize',22);
        set(gca,'XTick',d_range/d_dist*labelx_dist/d_dist+0.5);
        set(gca,'XTickLabel',0:labelx_dist:max(d_range));
        set(gca,'YTick',h_range/h_dist*labely_dist/h_dist+0.5);
        set(gca,'YTickLabel',-max(h_range)/2:labely_dist:max(h_range)/2);
        hcb = colorbar;
        constant_colorbar = 1;
%         maximum_limcolor = 6e-03;
        maximum_limcolor = 2e+03;
        if constant_colorbar == 1 % keeping the colorbar constant
            cmax = maximum_limcolor;
        else 
            cmax = max(max(count));
            %cmax = max(max(medint));
        end
        cmin=0;
      
       caxis([cmin cmax]);  
       [map] = ColormapAll(colormap);
       colormap(map); 
       density_pth = strcat(save_add_path,"den-S",string(index),sample_type{idx},'.png');
       saveas(gcf,density_pth)
 

        
        
    end
%     % Start to calculate the correlation scores and select the samples
% 
% %     index_max=0;
%     corr_mat = zeros(max_inx,max_inx);
%     for cor_i=1:max_inx
% %         sum_cur=0;
%         for cor_j=1:max_inx
%             r = corrcoef(density_vector(:,cor_i), density_vector(:,cor_j));
%             corr_mat(cor_i,cor_j)=r(1,2);
%         end
%     end
%     [sorted_corr_sum, sample_index] = sort(sum(corr_mat),'descend');
%     % Write the result of sample selection into txt
%     Sample_sel_info=strcat(save_add_path,className,'-sample_selection_info.txt');
%    % open file identifier
%     fid=fopen(Sample_sel_info,'w');
%     fprintf(fid, ['The sorted samples w.r.t. the sum of correlation'  '\n']);
%     fclose(fid);
%     T = table(sample_index', sorted_corr_sum', 'VariableNames', { 'Sample Index', 'Sum of Correlation'} );
%     % Write data to text file
%     writetable(T, Sample_sel_info)
%     disp('The ranked correlation scores are saved in the path:')
%     disp(Sample_sel_info)
     
    
end

