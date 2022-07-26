cd D:\NutCloud\我的坚果云\RWTH-Study\Hiwi\20.07\4_Microgel_plotter_V1_density
clc
clear all;
% close all
% 21 33 35 38 43 53
type = 'gg';
index = 5;
if index==1
    temperature = '21';
elseif index==2
    temperature = '33';
elseif index==3
    temperature = '35';
elseif index==4
    temperature = '38';
elseif index==5
    temperature = '43';
else
    temperature = '53';
end


add = strcat("./PAINT_DiffTemp/toy_example/results/results_mtemp/",type,string(index),".ply");
all_localiz  = dlmread(add);
add_info = strcat("./PAINT_DiffTemp/toy_example/results/results_mtemp/info",string(index),".ply");

info  = dlmread(add_info);


all_localiz  = all_localiz  .* (info(2,:)-info(3,:)) + info(3,:)+ info(1,:);


% all_localiz  = all_localiz  .* (info(2,:)-info(3,:)) + info(3,:)+ info(1,:);
intens=all_localiz(:,4);
% Yiwei: Find what is the outlier
UpperSolva = 400;
%all_localiz = all_localiz(intens<UpperSolva,:);% if a data point has the Solva>400, throw it away
all_localiz(intens>UpperSolva,4) = UpperSolva ;
zvec=all_localiz(:,3);
intens=all_localiz(:,4);
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

% cloud = pointCloud(data);
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
    
    
         
    data(:,5)= intens; % for including the solvatochromism
    % Prepare the density
%     count=zeros(length(d_range)-1,length(h_range)-1);
%     data=data(data(:,3)>0,:); % if a data point is out of the x range, throw it away
%     data=data(data(:,4)>0,:);% if a data point is out of the y range, throw it away
%     
%     for i=1:size(data,1)
%         try
%             count(data(i,3),data(i,4))=count(data(i,3),data(i,4))+1; 
%             %count{(data(i,3),data(i,4))}=data(i,3),data(i,4),data(i,5); 
%            %count(data(i,5))=count(data(i,5))+1
%         catch
%             i;
%         end
%     end
%     
    



% Plotting the localization density wrt median of I_ratio
   medint=zeros(length(d_range)-1,length(h_range)-1); %defining the size of "medint" to determine the medians of intensity ratios in specific pixels 
  
   for d=1:size(d_range',1)-1
   for h=1:size(h_range',1)-1
       ProcessData=data; % just copying the data to ProcessData variable
       
       % deleting all the values which are out of specific range
       for k=1:size(data,1)
            if (ProcessData(k,1)<=d_range(d)||ProcessData(k,1)>d_range(d+1)) || (ProcessData(k,2)<=h_range(h)||ProcessData(k,2)>h_range(h+1))
                ProcessData(k,5)=0;
            end
       end
       % copying all I_ratios to a new vector
%        for i=size(ProcessData,5)
%            OnlyInt = ProcessData(:,5);
%        end
       OnlyInt = ProcessData(:,5);
       OnlyInt(OnlyInt==0)=[];%% deleting the I_ratios = 0
       medint(d,h)=median(OnlyInt);% Replace the zeros with either NAN or median of I_ratio
   end
   end
%medint(isnan(medint))=0;% Replacing NAN values with zeros
medint=medint';
    % visualization
    figure_height_in_pixel=600;
    
    figure;
    set(gcf,'Position',[20 20 figure_height_in_pixel max(h_range)/max(d_range)*figure_height_in_pixel]);
%     count=count'./repmat(V,size(count,2),1);
%     countbin = count;

    %imagesc(count); % To plot wrt density values
    %Normalize the medint
    medint = (medint-min(min(medint)))./(max(max(medint))-min(min(medint))).*0.4+0.4;
    imagesc(medint); % To plot wrt I_ratios
    
    axis equal % suggestion by Eric to prevent elonagation in z
    
   
   
    set(gca,'Ydir','Normal');
    set(gca,'FontSize',22);
    xlim([0.5 max(d_range/d_dist)+1]);
    ylim([0.5 max(h_range/h_dist)+1]);
    if strcmp(temperature,'38')
         title('\fontsize{70}\fontname{Calibri} g','position',[-7.424 82.9050 0])
    elseif strcmp(temperature,'21')
        title('\fontsize{70}\fontname{Calibri} d','position',[-7.424 82.9050 0])
    elseif strcmp(temperature,'33')
        title('\fontsize{70}\fontname{Calibri} e','position',[-7.424 82.9050 0])
    elseif strcmp(temperature,'35')
        title('\fontsize{70}\fontname{Calibri} f','position',[-7.424 82.9050 0])
    elseif strcmp(temperature,'43')
        title('\fontsize{70}\fontname{Calibri} h','position',[-7.424 82.9050 0])
    else
        title('\fontsize{70}\fontname{Calibri} i','position',[-7.424 82.9050 0])
    end
   
    %title('\fontsize{70}\fontname{Calibri} 21 C', 'position',[20.424 82.9050 0]); 
    xlabel('Dist. from symmetry axis / nm','FontSize',22);
    ylabel('Relative z / nm','FontSize',22);
    set(gca,'XTick',d_range/d_dist*labelx_dist/d_dist+0.5);
    set(gca,'XTickLabel',0:labelx_dist:max(d_range));
    set(gca,'YTick',h_range/h_dist*labely_dist/h_dist+0.5);
    set(gca,'YTickLabel',-400:labely_dist:400);
%     set(gcf,'color','w');%Ashvini test
%         set(gca,'YTick',-300:10:300);
%         set(gca,'YTickLabel',0:labely_dist:max(h_range));
hcb = colorbar;

constant_colorbar = 1;
maximum_limcolor = 0.8;
if constant_colorbar == 1 % keeping the colorbar constant
    cmax = maximum_limcolor;
else 
    %cmax = max(max(count));
    cmax = max(max(medint));
end
cmin=0.4;
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
        
        
% elseif LowerIntensity == 0 && UpperIntensity == 50
%    [map] = ColormapAll(mapblue);
%    colormap(map);
% elseif LowerIntensity == 50 && UpperIntensity == 100
%    [map] = ColormapAll(mapgreen);
%    colormap(map);
% elseif LowerIntensity == 100 && UpperIntensity == 150
%    [map] = ColormapAll(mapgreenyellow);
%    colormap(map);
% elseif LowerIntensity == 150 && UpperIntensity == 200
%    [map] = ColormapAll(mapred);
%    colormap(map);
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
   



   
%fname_2d = strcat("./results/visualization/",temperature,"/",temperature,'C','_2d-solv-distr.png');
fname_2d = strcat("./results/visualization/",temperature,"/",temperature,'C','_gg_2d-solv-distr.png');
saveas(gcf,fname_2d)
    %set(gca,'XTickLabel', d_range-d_dist/2);
    %set(gca,'YTickLabel', h_range+h_dist/2);

    %% end Ashvini code ends