% function [labelx_dist,labely_dist, d_dist, h_dist, d_range, h_range,enzcount,Ind_crossecn,Ind_crossecn_mean,ergebnis,hist_binmids,hist_binwidth, LowerIntensity,UpperIntensity,allcount]=localizationdensity2D_v2(baseNamecsv, folder,Rh)
close all;
clc;
clear all;

%3D Analysis of Point Cloud Data, taken from VisP. 18-06-15 included:
%%Offset of Datapoints
%Loading the data in 3d-Format

%% Input for simpler testing
  
%   add = "./M33/NR_1.3d";
%   folder = "./M33/";
%   baseNamecsv = 'NR_1.3d';
Rh = 300;
temperature = '53';
micro_index = '15';
add = strcat("./PAINT_DiffTemp/Core-shell/",temperature,"C/Microgel_plotter_v2_solvatochromism/",micro_index,".3d");
%% User Input 

HightOfPlot = 1000; %Size of relative z-values in px 

interval = 1126; % interval is the amount of colors inserted for the color map
%maximum_limcolor= 0.8; %The maximum limit you want to set with constant colorbar
maximum_limcolor= 1*10^-4; %for density
constant_colorbar=0; % if this is 1 then constant maximum value for the colorbar 

LowerIntensity = 0;
UpperIntensity = 1;

PlotActiveCenter = 0;   %Only for plotting active centers e.g. Enzymes - set to 1


%% User Input end! 
%% 

%fullFileNamecsv = fullfile(folder, baseNamecsv);

%txt = (['Reading ',baseNamecsv,' . This can take some time ...']);
txt = (['Reading ',add,' . This can take some time ...']);
disp(txt);
% fullbaseNamecsv = fullfile(folder, baseNamecsv);
% [~,baseNamecsv,~] = fileparts(fullbaseNamecsv);
%all_localiz = dlmread (fullbaseNamecsv,'	');% reads numeric data from ASCII
% all_localiz = dlmread (fullbaseNamecsv,',');% reads numeric data from ASCII

%add = "./../M33/NR_1.3d";

all_localiz  = dlmread(add);

%% Extracting x,y,z Data to single vectors
xvec = all_localiz(:,1);
yvec = all_localiz(:,2);
zvec = all_localiz(:,3);
intvec = all_localiz(:,4);
xyz = all_localiz(:,1:3);

data = all_localiz;





%% start To density distribution in different regions of each microgel using 2D plot--Dominik und Ashvini
% Part1 To determine distance from the axis
xcoord=all_localiz(:,1);
deltax = xcoord-median(xcoord);
ycoord=all_localiz(:,2);
deltay = ycoord-median(ycoord);
zcoord = all_localiz(:,3);
deltaz = zcoord-median(zcoord);

RelCoords = [deltax(:) deltay(:) deltaz(:) intvec(:)];

    %% calculations and preparation of histogram ranges
    d_axis = sqrt(deltax.^2+deltay.^2);
    h = (zvec-median(zvec))+HightOfPlot/2;

    clear count V d_range h_range;

    %% calculations and preparation of histogram ranges

    d_dist=20; %Size per pixel in x (nm/px)
    h_dist=20; %Size per pixel in y (nm/px)
    
    labelx_dist=100;
    labely_dist=100;

    d_range=0:d_dist:400;
    h_range=0:h_dist:HightOfPlot;



    %% calculation of the 2D histogram
    for i=1:length(d_range)-1
        V(i)=((d_range(i)+d_dist)^2-(d_range(i))^2)*pi*h_dist;
    end

    data=[d_axis h]; %Creating Localization data with coordinates with respect to center of microgel (median)
    for i=1:length(d_range)-1
        data((data(:,1)>d_range(i))&(data(:,1)<=d_range(i+1)),3)=i;
    end
    for i=1:length(h_range)-1
        data((data(:,2)>h_range(i))&(data(:,2)<=h_range(i+1)),4)=i;  
    end
%     
    data(:,5) = intvec;
%     
% %     save([folder, baseNamecsv,'_XY-Radial.mat'],'data');       %% Intermediatesave 
% %     writematrix(data,[folder, baseNamecsv,'_XY-Radial.csv']);
%          
%     %data(:,5)= intens; % for including the solvatochromism
    data=data(data(:,3)>0,:); % if a data point is out of the x range, throw it away
    data=data(data(:,4)>0,:);% if a data point is out of the y range, throw it away

    enzcount=zeros(length(d_range)-1,length(h_range)-1);
    for i=1:size(data,1)
        if data(i,5) == 2
            try
                enzcount(data(i,3),data(i,4))=enzcount(data(i,3),data(i,4))+1; 
            catch
                i;
            end
        else
            enzcount(data(i,3),data(i,4))=0;
        end
    end

    enzcount=enzcount'./repmat(V,size(enzcount,2),1);
    



%% Calculation of the localization density wrt median of I_ratio
   allcount=zeros(length(d_range)-1,length(h_range)-1); %defining the size of "medint" to determine the medians of intensity ratios in specific pixels 
  
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
       for i=size(ProcessData,5)
           OnlyInt = ProcessData(:,5);
       end
       OnlyInt(OnlyInt==0)=[];%% deleting the I_ratios = 0
       allcount(d,h)=median(OnlyInt);% Replace the zeros with either NAN or median of I_ratio
   end
   end
    allcount(isnan(allcount))=0;% Replacing NAN values with zeros
    allcount(allcount > 0) = -1;

allcount=allcount';

%% For a custom colorbar 
colormapgray = [0 0 0
    1 1 1];
[mapgray] = colormapgray;
[map] = ColormapAll(colormap);
close all;

    %% Setting figure data (pos, size,...)
    figure_height_in_pixel=500;
    
    figure(1);
%    set(gcf,'Position',[20 20 figure_height_in_pixel+100 max(h_range)/max(d_range+100)*figure_height_in_pixel]);
    set(gcf,'Position',[20 20 figure_height_in_pixel+100 max(h_range)/max(d_range+100)*figure_height_in_pixel]);
    


    %% Plotting graphs in figure
    ax1 = axes;
    title('\fontsize{40}\fontname{Arial} Densityplot')
    xlabel('Dist. from cylinder axis / nm', 'FontSize',22);
    ylabel('Relative z / nm', 'FontSize',22);
    
    Image1 = imagesc(ax1,allcount); % To plot all pixel with localization in gray
    Image1.AlphaData = 1;
    axis equal 
    set(gca,'Ydir','Normal');
    set(gca,'FontSize',22);
    xlim([0.5 max(d_range/d_dist)+1]);
    ylim([0.5 max(h_range/h_dist)+1]);
    hold all;
    

    
    
    ax2 = axes;
    Image2 = imagesc(ax2,enzcount); % To plot wrt density values

    linkaxes([ax1,ax2]);
    ax2.Visible = 'off'; 
    ax2.XTick = []; 
    ax2.YTick = []; 
    
    hold off;
    axis equal 
    set(gca,'Ydir','Normal');
    xlim([0.5 max(d_range/d_dist)+1]);
    ylim([0.5 max(h_range/h_dist)+1]);
    

    
    %% Setting Colorbar
   
    cmin = 0; 
     
    if constant_colorbar == 1 % keeping the colorbar constant
        cmax = maximum_limcolor;
    else 
        %cmax = max(max(count));
        cmax = max(max(enzcount));
    end
    
    colormap(ax1,mapgray);
    colormap(ax2,map);
 %   Image2.AlphaData = 0.8;
    hcb = colorbar(ax2);
    hcb1 = colorbar(ax1);
    set(hcb1, 'Visible', 'off');
    
    caxis([cmin cmax]);
    if PlotActiveCenter == 1
        set(hcb, 'Ticks', 0, 'TickLabels' , ' ');
        set(get(hcb,'label'),'string','No localization found         microgel                 enzyme           ');
    else 
        set(get(hcb,'label'),'string','Localization density / nm³');
        set(hcb,'FontSize',22,'Fontname','Arial');
    end

    

    linkaxes([ax1,ax2]);
    axis equal;
    set(ax1,'XTick',d_range/d_dist*labelx_dist/d_dist+0.5);
    set(ax1,'XTickLabel',0:labelx_dist:max(d_range));
    set(ax1,'YTick',h_range/h_dist*labely_dist/h_dist+0.5);
    set(ax1,'YTickLabel',-(HightOfPlot/2):labely_dist:(HightOfPlot/2)); 
    set(ax2,'XTick',d_range/d_dist*labelx_dist/d_dist+0.5);
    set(ax2,'XTickLabel',0:labelx_dist:max(d_range));
    set(ax2,'YTick',h_range/h_dist*labely_dist/h_dist+0.5);
    set(ax2,'YTickLabel',-(HightOfPlot/2):labely_dist:(HightOfPlot/2));
    axis equal;
    hold on
    
    if PlotActiveCenter == 1
        title('\fontsize{40}\fontname{Arial} Active Centers')
    else    
        title('\fontsize{40}\fontname{Arial} Densityplot')
    end

    
    
    %% Plotting final image with both plots

    
    
    if PlotActiveCenter == 1
        saveas(figure(1),[folder, baseNamecsv,'_ActiveCenters.tif']);
        saveas(figure(1),[folder, baseNamecsv,'_ActiveCenters.fig']);        
    else
        saveas(figure(1),[folder, baseNamecsv,'_2D-density-distr.tif']);
        saveas(figure(1),[folder, baseNamecsv,'_2D-density-distr.fig']);
    end
 close all;

