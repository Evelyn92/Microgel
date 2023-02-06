%
clear all;
clc;
close all;
userpath('clear')

path = "D:\forClone\Hiwi\Microgel\4_Microgel_plotter_V1_density\results_mtemp\results_mtemp_hilo";
for tmp_inx = 1:6

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
data_type = 'gg';%'gg' or '_go'

for sample_inx = 1:4

add = strcat(path,'\sample_',string(sample_inx),data_type,string(tmp_inx),".3d");
all_localiz  = dlmread(add);
add_info = strcat(path,'\sample_',string(sample_inx),'info',string(tmp_inx),".3d");
info  = dlmread(add_info);
all_localiz  = all_localiz  .* (info(2,:)-info(3,:)) + info(3,:)+ info(1,:);

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
d_axis = sqrt(deltax.^2+deltay.^2);
h = (zvec-median(zvec))+400;
clear count V d_range h_range;

d_dist=10;
h_dist=10;

labelx_dist=100;
labely_dist=100;

d_range=0:d_dist:400;
h_range=0:h_dist:800;

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
    
    % visualization
    figure_height_in_pixel=600;
    
    figure;
    set(gcf,'Position',[20 20 figure_height_in_pixel max(h_range)/max(d_range)*figure_height_in_pixel]);
    count=count'./repmat(V,size(count,2),1);
%     countbin = count;

    imagesc(count); % To plot wrt density values
    %imagesc(medint); % To plot wrt I_ratios
    
    axis equal % suggestion by Eric to prevent elonagation in z
    
    set(gca,'Ydir','Normal');
    set(gca,'FontSize',22);
    xlim([0.5 max(d_range/d_dist)+1]);
    ylim([0.5 max(h_range/h_dist)+1]);
    title_info = strcat(tmp_info,'-',string(sample_inx));
    %title('\fontsize{30}\fontname{Calibri}'+title_info);

    %title('\fontsize{70}\fontname{Calibri} 21 C', 'position',[20.424 82.9050 0]); 
    
    set(gca,'XTick',d_range/d_dist*labelx_dist/d_dist+0.5);
    set(gca,'XTickLabel',0:labelx_dist:max(d_range),'fontsize',28);
    set(gca,'YTick',h_range/h_dist*labely_dist/h_dist+0.5);
    set(gca,'YTickLabel',-400:labely_dist:400, 'fontsize',28);
    xlabel('Dist. from symmetry axis / nm','FontSize',28);
    ylabel('Relative z / nm','FontSize',28);
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

        colormap gray;
        
else

   [map] = ColormapAll(colormap);
   colormap(map); 

end

fname_2d = strcat(path,'/single_2d','/2d_',data_type,tmp_info,'_sample',string(sample_inx),'.png');
% saveas(gcf,fname_2d)

end
end