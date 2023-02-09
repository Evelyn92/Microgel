% Sample selection for the original diff_temperature
cd ('D:\forClone\Hiwi\Microgel\4_Microgel_plotter_V1_density')
clc
clear all
close all
%%
for tmp_inx = 1:6
tmp_list = [21 33 35 38 43 53];

if tmp_inx == 1
        tmp_info = '21C';
        samp_list  = [16,11,6,5,10,2];
elseif tmp_inx == 2
    tmp_info = '33C';
    samp_list  = [18,15,7,14,1,11];
elseif tmp_inx == 3
    tmp_info = '35C';
    samp_list  = [20,17,16,1,2,18];
elseif tmp_inx == 4
    tmp_info = '38C';
    samp_list  = [16,3,6,8,7,11];
elseif tmp_inx == 5
    tmp_info = '43C';
    samp_list  = [3,14,12,7,10,13];
else
    tmp_info = '53C';
    samp_list  = [10,11,7,18,6,19];
end
% FOCTS ODS PEG GLAS
% className = 'PEG';
density_vector = [];
max_index = 20;
% if strcmp(className, 'GLAS')
%     max_index = 10;
% elseif strcmp(className, 'FOCTS')
%     max_index = 8;
% elseif strcmp(className, 'PEG')
%     max_index = 8;
% else
%     max_index = 16;
% end

for index = 1:6
% micro_index = '2';
micro_index = string(samp_list(index));
% add = strcat("../24/",className,'/NiPMAM/24_',className,'_NiPMAM_',string(index),'.3d');
add = strcat("./PAINT_DiffTemp/Core-shell/",tmp_info,"/Microgel_plotter_v2_solvatochromism/",micro_index,".3d"); 
% all_localiz  = dlmread(add);
%add = "./results/results_mtemp/gg5.ply";

all_localiz  = dlmread(add);
% add_info = "./results/results_mtemp/info5.ply";
% info  = dlmread(add_info);


% all_localiz  = all_localiz  .* (info(2,:)-info(3,:)) + info(3,:)+ info(1,:);
% intens=all_localiz(:,4);

zvec=all_localiz(:,3);
% intens=all_localiz(:,4);
xcoord=all_localiz(:,1);
deltax = xcoord-median(xcoord);
ycoord=all_localiz(:,2);
deltay = ycoord-median(ycoord);
% visualize one microgel in 3D space
data = all_localiz;
if(size(data,2)>3) 
    data = data(:,1:3);
end

data(:,3) = 1-data(:,3); %flip z dimension

cloud = pointCloud(data);
% pcshow(cloud);

% 
% fname = strcat("./results/visualization/",temperature,"/",temperature,'C_',micro_index,'.png');
% saveas(gcf,fname)

    % calculations and preparation of histogram ranges
d_axis = sqrt(deltax.^2+deltay.^2);
h = (zvec-median(zvec))+400;
clear count V d_range h_range;
    
    % create test data (not needed for real data)
%     n=100000;
% 
%     x=rand(n,1)*200-100;
%     y=rand(n,1)*200-100;
%     z=rand(n,1)*200;
% 
%     centerx=0;
%     centery=0;
%     centerz=100;

    % calculations and preparation of histogram ranges
%     d_axis=sqrt((x-centerx).^2+(y-centery).^2);
%     h=z;

    d_dist=10;
    h_dist=10;
    
    labelx_dist=100;
    labely_dist=100;

    d_range=0:d_dist:400;
    h_range=0:h_dist:800;
    

   % hist_h=hist(h,h_range);
   % hist_d=hist(d_axis,d_range);

    % calculation of the 2D histogram
    for i=1:length(d_range)-1
        V(i)=((d_range(i)+d_dist)^2-(d_range(i))^2)*pi*h_dist;
    end

    data=[d_axis h];
    for i=1:length(d_range)-1
        data((data(:,1)>d_range(i))&(data(:,1)<=d_range(i+1)),3)=i;
    end
    
    data2 = [d_axis h];
    for i=1:length(d_range)-1
        for j=1:length(data)
            if data2(j,1)>d_range(i)&&(data2(j,1)<=d_range(i+1))
             data2(j,3)=i;
            end
        end
    end    
    
    
    
    
    for i=1:length(h_range)-1
        data((data(:,2)>h_range(i))&(data(:,2)<=h_range(i+1)),4)=i;  
    end
    
    
         
%     data(:,5)= intens; % for including the solvatochromism
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
    
    



% Plotting the localization density wrt median of I_ratio
%    medint=zeros(length(d_range)-1,length(h_range)-1); %defining the size of "medint" to determine the medians of intensity ratios in specific pixels 
%   
%    for d=1:size(d_range',1)-1
%    for h=1:size(h_range',1)-1
%        ProcessData=data; % just copying the data to ProcessData variable
%        
%        % deleting all the values which are out of specific range
%        for k=1:size(data,1)
%             if (ProcessData(k,1)<=d_range(d)||ProcessData(k,1)>d_range(d+1)) || (ProcessData(k,2)<=h_range(h)||ProcessData(k,2)>h_range(h+1))
%                 ProcessData(k,5)=0;
%             end
%        end
%        % copying all I_ratios to a new vector
%        for i=size(ProcessData,5)
%            OnlyInt = ProcessData(:,5);
%        end
%        OnlyInt(OnlyInt==0)=[];%% deleting the I_ratios = 0
%        medint(d,h)=median(OnlyInt);% Replace the zeros with either NAN or median of I_ratio
%    end
%    end
% %medint(isnan(medint))=0;% Replacing NAN values with zeros
% medint=medint';
    % visualization
    figure_height_in_pixel=600;
    
    figure;
    set(gcf,'Position',[20 20 figure_height_in_pixel max(h_range)/max(d_range)*figure_height_in_pixel]+90);
    count=count'./repmat(V,size(count,2),1);
%     countbin = count;
    density_vector = [density_vector reshape(count,[],1)];
    imagesc(count); % To plot wrt density values
    %imagesc(medint); % To plot wrt I_ratios
    
    axis equal % suggestion by Eric to prevent elonagation in z
    
   
   
    set(gca,'Ydir','Normal');
    set(gca,'FontSize',30);
    xlim([0.5 max(d_range/d_dist)+1]);
    ylim([0.5 max(h_range/h_dist)+1]);
%     if strcmp(temperature,'38')
%          title('\fontsize{70}\fontname{Calibri} g','position',[-7.424 82.9050 0])
%     elseif strcmp(temperature,'21')
%         title('\fontsize{70}\fontname{Calibri} d','position',[-7.424 82.9050 0])
%     elseif strcmp(temperature,'33')
%         title('\fontsize{70}\fontname{Calibri} e','position',[-7.424 82.9050 0])
%     elseif strcmp(temperature,'35')
%         title('\fontsize{70}\fontname{Calibri} f','position',[-7.424 82.9050 0])
%     elseif strcmp(temperature,'43')
%         title('\fontsize{70}\fontname{Calibri} h','position',[-7.424 82.9050 0])
%     else
%         title('\fontsize{70}\fontname{Calibri} i','position',[-7.424 82.9050 0])
%     end
   
    %title('\fontsize{70}\fontname{Calibri} 21 C', 'position',[20.424 82.9050 0]); 
    xlabel('Dist. from symmetry axis / nm','FontSize',30);
    ylabel('Relative z / nm','FontSize',30);
    set(gca,'XTick',d_range/d_dist*labelx_dist/d_dist+0.5);
    set(gca,'XTickLabel',0:labelx_dist:max(d_range));
    set(gca,'YTick',h_range/h_dist*labely_dist/h_dist+0.5);
    set(gca,'YTickLabel',-400:labely_dist:400);
%     set(gcf,'color','w');%Ashvini test
%         set(gca,'YTick',-300:10:300);
%         set(gca,'YTickLabel',0:labely_dist:max(h_range));
hcb = colorbar;

constant_colorbar = 1;
maximum_limcolor = 6e-03;
if constant_colorbar == 1 % keeping the colorbar constant
    cmax = maximum_limcolor;
else 
    cmax = max(max(count));
    %cmax = max(max(medint));
end
cmin=0;
% cmin = -(cmax/(interval+200))+0.4; % correction for colorbar to start the blue color at almost 0 
caxis([cmin cmax]);
%caxis([0.4 0.8]);
LowerIntensity = 0;
UpperIntensity = 1;
if LowerIntensity == 0 && UpperIntensity == 2000
    % colormap gray;
%       [map] = ColormapAll(colormap); %For jet
%       colormap(map);
%       
      [map] = Colormapblue2red(colormap);% Use this for Radial solvato values
   colormap(map); 

elseif LowerIntensity == 0 && UpperIntensity == 200
%     [map] = ColormapAll(colormap);
%      colormap(map);
        colormap gray;
        
        
else
   % colormap hot;
    %[map] = Colormapblue2red(colormap);% Use this for Radial solvato values
   [map] = ColormapAll(colormap);
   colormap(map); 
%    [map] = ColormapAll(colormap);
%    colormap(map);
end

 %Namemat = baseNamecsv;
  %save(Namemat, 'data','count', 'V');
   
save_fd = 'D:\forClone\Hiwi\Microgel\4_Microgel_plotter_V1_density\PAINT_DiffTemp\';
fname_2d = strcat(save_fd,tmp_info,'_',string(index),'samp',micro_index,'_2d-density-distr.png');
saveas(gcf,fname_2d)
    %set(gca,'XTickLabel', d_range-d_dist/2);
    %set(gca,'YTickLabel', h_range+h_dist/2);
end
    % end Ashvini code ends
    close all;
end
%     %% Calculate the correlation
close all;
index_max=0;
corr_mat = zeros(20,20);
for cor_i=1:20
    sum_cur=0;
    for cor_j=1:20
        r = corrcoef(density_vector(:,cor_i), density_vector(:,cor_j));
        corr_mat(cor_i,cor_j)=r(1,2);
%         sum_cur=sum_cur+r[0, 1]
%         if(sum_cur>sum_hist):
%             sum_hist=sum_cur
%             index_max=i
%         end
    end
end
[sorted_corr_sum, sample_index] = sort(sum(corr_mat),'descend');
% Write the result of sample selection into txt
Sample_sel_info=strcat(string(temperature),'C_sample_selection_info.txt');
% %open file identifier
% fid=fopen(Sample_sel_info,'w');
% fprintf(fid, ['The sorted samples w.r.t. the sum of correlation'  '\n']);
% fclose(fid);

T = table(sample_index', sorted_corr_sum', 'VariableNames', { 'Sample Index', 'Sum of Correlation'} );
% Write data to text file
cd D:\forClone\Hiwi\Microgel\4_Microgel_plotter_V1_density
writetable(T, Sample_sel_info)