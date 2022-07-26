% function [labelx_dist,labely_dist, d_dist, h_dist, d_range, h_range,count,Ind_crossecn,Ind_crossecn_mean,ergebnis,hist_binmids,hist_binwidth, LowerIntensity,UpperIntensity,medint]=localizationdensity2D_v2(baseNamecsv, folder,Rh)
% close all
cd D:\NutCloud\我的坚果云\RWTH-Study\Hiwi\20.07\4_Microgel_plotter_V1_density
clear all
clc
%3D Analysis of Point Cloud Data, taken from VisP. 18-06-15 included:
%%Offset of Datapoints
%Loading the data in 3d-Format

%% Input for simpler testing
%add = "./../M33/NR_2.3d";
temperature = '21';
micro_index = '1';
add = strcat("./PAINT_DiffTemp/Core-shell/",temperature,"C/Microgel_plotter_v2_solvatochromism/",micro_index,".3d");
Rh = 300;
  
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
% 
% txt = (['Reading ',baseNamecsv,' . This can take some time ...']);
% disp(txt);
% fullbaseNamecsv = fullfile(folder, baseNamecsv);
% [~,baseNamecsv,~] = fileparts(fullbaseNamecsv);
%all_localiz = dlmread (fullbaseNamecsv,'	');% reads numeric data from ASCII %% uncomment
% all_localiz = dlmread (fullbaseNamecsv,',');% reads numeric data from ASCII %% DELETE!


all_localiz  = dlmread(add);

all_localiz(:,4) = [1];

%% Extracting x,y,z Data to single vectors
xvec = all_localiz(:,1);
yvec = all_localiz(:,2);
zvec = all_localiz(:,3);
intvec = all_localiz(:,4);
xyz = all_localiz(:,1:3);





%% start To density distribution in different regions of each microgel using 2D plot--Dominik und Ashvini
% Part1 To determine distance from the axis
xcoord=all_localiz(:,1);
deltax = xcoord-median(xcoord);
ycoord=all_localiz(:,2);
deltay = ycoord-median(ycoord);

%% Remove outlaiers

points = all_localiz(:,1:3);
meanvecs = mean(points);
points_dis = points-meanvecs;
points(:, 4) = sqrt(points_dis(:,1).^2+points_dis(:,2).^2);
points = sortrows(points, 4);
points = points(1:floor(0.99*length(points)), :);


NewRange = (2 - (-2));  
x = points(:,1) - mean(points(:,1));
y = points(:,2) - mean(points(:,2));
z = points(:,3) - mean(points(:,3));
x = (((x) * NewRange) / (max(x)-min(x)));
y = (((y) * NewRange) / (max(y)-min(y)));
z = (((z) * NewRange) / (max(z)-min(z)));
cloud = pointCloud([x y z]);
figure()
[X,Y,Z] = sphere(15);
surf (X, Y, Z,'FaceAlpha',0.0)
axis equal
% hold on
pcshow(cloud)
set(gca,'Color','w')
set(gcf,'color','w');

% points = points - mean(points);
d_axis = points(:, 4);
h = (points(:, 3)-median(points(:, 3)))+HightOfPlot/2;




cloud = pointCloud(points(:,1:3));
pcshow(cloud)
fname = strcat("./results/visualization/",temperature,"/",micro_index,'.png');
saveas(gcf,fname)

    %% calculations and preparation of histogram ranges
%     d_axis = sqrt(deltax.^2+deltay.^2);
%     h = (zvec-median(zvec))+HightOfPlot/2;

     clear count V d_range h_range;

    %% calculations and preparation of histogram ranges

    d_dist=20; %Size per pixel in x (nm/px)
    h_dist=20; %Size per pixel in y (nm/px)
    
    labelx_dist=100;
    labely_dist=100;

    d_range=0:d_dist:600;
    h_range=0:h_dist:HightOfPlot;
    


    %% calculation of the 2D histogram
    for i=1:length(d_range)-1
        V(i)=((d_range(i)+d_dist)^2-(d_range(i))^2)*pi*h_dist;
    end

    data=[d_axis h];
    for i=1:length(d_range)-1
        data((data(:,1)>d_range(i))&(data(:,1)<=d_range(i+1)),3)=i;
    end
    for i=1:length(h_range)-1
        data((data(:,2)>h_range(i))&(data(:,2)<=h_range(i+1)),4)=i;  
    end
    
    
         
    %data(:,5)= intens; % for including the solvatochromism
    count=zeros(length(d_range)-1,length(h_range)-1);
    if PlotActiveCenter == 1       %Only for plotting active centers e.g. Enzymes
        data(:,5) = intvec;
    end
    data=data(data(:,3)>0,:); % if a data point is out of the x range, throw it away
    data=data(data(:,4)>0,:);% if a data point is out of the y range, throw it away
    
    for i=1:size(data,1)
        try
            count(data(i,3),data(i,4))=count(data(i,3),data(i,4))+1; 
        catch
            i;
        end
    end
    
    



%% Calculation of the localization density wrt median of I_ratio
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
       for i=size(ProcessData,5)
           OnlyInt = ProcessData(:,5);
       end
       OnlyInt(OnlyInt==0)=[];%% deleting the I_ratios = 0
       medint(d,h)=median(OnlyInt);% Replace the zeros with either NAN or median of I_ratio
   end
   end
medint(isnan(medint))=0;% Replacing NAN values with zeros
if PlotActiveCenter == 1    %Only for plotting active centers e.g. Enzymes
    medint(medint > 1) = 2;
end
medint=medint';

%% For a custom colorbar 



if LowerIntensity == 0 && UpperIntensity == 2000
    [map] = ColormapAll(colormap); %For jet
    colormap(map);     
elseif LowerIntensity == 0 && UpperIntensity == 50
    [map] = ColormapAll(mapblue);
    colormap(map);
elseif LowerIntensity == 50 && UpperIntensity == 100
    [map] = ColormapAll(mapgreen);
    colormap(map);
elseif LowerIntensity == 100 && UpperIntensity == 150
    [map] = ColormapAll(mapgreenyellow);
    colormap(map);
elseif LowerIntensity == 150 && UpperIntensity == 200
    [map] = ColormapAll(mapred);
    colormap(map);
else
    [map] = ColormapAll(colormap);
    colormap(map);
end
% close all;
    
if PlotActiveCenter == 1
    colormapnew = [1 1 1  
                   1 0 0
                   0 0 1];
[map] = colormapnew;
end
    

    %% visualization
    figure_height_in_pixel=500;
    
    figure();
    set(gcf,'Position',[20 20 figure_height_in_pixel+100 max(h_range)/max(d_range+100)*figure_height_in_pixel]);
    count=count'./repmat(V,size(count,2),1);
    
    if PlotActiveCenter == 1
        imagesc(medint); % To plot wrt I_ratios
    else
        imagesc(count); % To plot wrt density values
    end

    cmin = 0; 
    
    
    
    
    if constant_colorbar == 1 % keeping the colorbar constant
        cmax = maximum_limcolor;
    else 
        %cmax = max(max(count));
        cmax = max(max(medint));
    end

        
        
    
    axis equal 
    set(gca,'Ydir','Normal');
    set(gca,'FontSize',22);
    xlim([0.5 max(d_range/d_dist)+1]);
    ylim([0.5 max(h_range/h_dist)+1]);
    if PlotActiveCenter == 1
        title('\fontsize{40}\fontname{Arial} Active Centers')
    else    
        title('\fontsize{40}\fontname{Arial} Densityplot')
    end
    xlabel('Dist. from cylinder axis / nm', 'FontSize',22);
    ylabel('Relative z / nm', 'FontSize',22);
    set(gca,'XTick',d_range/d_dist*labelx_dist/d_dist+0.5);
    set(gca,'XTickLabel',0:labelx_dist:max(d_range));
    set(gca,'YTick',h_range/h_dist*labely_dist/h_dist+0.5);
    set(gca,'YTickLabel',-(HightOfPlot/2):labely_dist:(HightOfPlot/2));

    caxis([cmin cmax]);
    colormap(map);
    hcb = colorbar;
    if PlotActiveCenter == 1
        set(hcb, 'Ticks', 0, 'TickLabels' , ' ');
        set(get(hcb,'label'),'string','No localization found         microgel                 enzyme           ');
    else 
        set(get(hcb,'label'),'string','Localization density / nm³');
    end
    

    
    if PlotActiveCenter == 1
        saveas(figure(1),[folder, baseNamecsv,'_ActiveCenters.tif']);
        saveas(figure(1),[folder, baseNamecsv,'_ActiveCenters.fig']);        
    else
        saveas(figure(1),[folder, baseNamecsv,'_2D-density-distr.tif']);
        saveas(figure(1),[folder, baseNamecsv,'_2D-density-distr.fig']);
    end
 close all;

