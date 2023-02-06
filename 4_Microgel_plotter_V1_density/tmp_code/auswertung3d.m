function [labelx_dist,labely_dist, d_dist, h_dist, d_range, h_range,count,Ind_crossecn,Ind_crossecn_mean,ergebnis,hist_binmids,hist_binwidth, LowerIntensity,UpperIntensity,medint]=auswertung3d(baseNamecsv, folder,Rh)
close all
%3D Analysis of Point Cloud Data, taken from VisP. 18-06-15 included:
%%Offset of Datapoints
%Loading the data in 3d-Format

%% Input for simpler testing
  folder = 'D:\Lukas\Sciebo\Promotion\Publikationen\LocalizingEnzymes\Data\Measurement-MicrogelsData\20200108_Julia\';
  baseNamecsv = 'All_M33_MGE1-8+M36.3d';
  Rh = 500;


%% Input end

%% 



%fullFileNamecsv = fullfile(folder, baseNamecsv);

txt = (['Reading ',baseNamecsv,' . This can take some time ...']);
disp(txt);
fullbaseNamecsv = fullfile(folder, baseNamecsv);
[~,baseNamecsv,~] = fileparts(fullbaseNamecsv);
%All localizations saved in all_localiz
all_localiz = dlmread (fullbaseNamecsv,'	');% reads numeric data from ASCII

%% User Input

% prompt = 'What do you expect to be Rh? (nm)'; % uncomment this line if you what only 1 microgel plot
% Rh = input(prompt);  % uncomment this line if you what only 1 microgel plot



LowerIntensity=0; % Lower border for Intensity cutoff (e.g. between 0-2000)
UpperIntensity=2001; % Upper border for Intensity cutoff (e.g. between 0-2000)



LowerLimit = 0;  % Lower Limit for climits - standard set to 0
UpperLimit = 1;  % Upper Limit for climits - standard set to 20

%interval = 257; % interval is the amount of colors inserted for the color map-redtoblue
%interval = 1126; % interval is the amount of colors inserted for the color map-jet
maximum_limcolor= 2; %The maximum limit you want to set with constant colorbar
%maximum_limcolor= 0.8; %The maximum limit you want to set with constant colorbar
constant_colorbar=0; % if this is 1 then constant maximum value for the colorbar else the maximum value changes for every data set
%%User Input end

% %% Sort out all unwanted data
% if LowerIntensity==0 && UpperIntensity==0
%     all_localiz = all_localiz;
% elseif LowerIntensity~=0 && UpperIntensity==0
%     all_localiz = all_localiz;
% else 
% 
% all_localiz(all_localiz(:,4)<LowerLimit*100,4)= 20; % To assign data < =0.4
% all_localiz(all_localiz(:,4)>UpperLimit*100,4)= 120; % To assign data >=0.8
% 
% end




%% Extracting x,y,z Data to single vectors
xvec = all_localiz(:,1);
%xvec=flip(xvec);
yvec = all_localiz(:,2);
%yvec=flip(yvec);
zvec = all_localiz(:,3);
%zvec=flip(zvec);
intens = all_localiz(:,4);
xyz = all_localiz(:,1:3);




N_microgel_local=length(all_localiz);
%% 
%Subtracting the offset in plot

if isempty(zvec(zvec<0))==0  
   disp('z contains negative elements.');
   zvec = zvec+abs(min(zvec))+1;
   disp('z-values shifted to offset = 0 .');
else
  disp('z: all positive.');
end
xvec=(xvec-min(xvec)+1)-50;
yvec=yvec-min(yvec)+1;

%% start To density distribution in different regions of each microgel using 2D plot--Dominik und Ashvini
% Part1 To determine distance from the axis
xcoord=all_localiz(:,1);
deltax = xcoord-median(xcoord);
ycoord=all_localiz(:,2);
deltay = ycoord-median(ycoord);
intens=all_localiz(:,4)/100;
%d_axis = sqrt(deltax.^2+deltay.^2);
%h = zvec;
% dh = [];
% dh = [daxis,h];
% scatter(daxis,h,'filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
% xlabel('d_{axis}');
% ylabel('h');

    %% calculations and preparation of histogram ranges
    d_axis = sqrt(deltax.^2+deltay.^2);
    h = (zvec-median(zvec))+400;

     clear count V d_range h_range;
    
    %% create test data (not needed for real data)
%     n=100000;
% 
%     x=rand(n,1)*200-100;
%     y=rand(n,1)*200-100;
%     z=rand(n,1)*200;
% 
%     centerx=0;
%     centery=0;
%     centerz=100;

    %% calculations and preparation of histogram ranges
%     d_axis=sqrt((x-centerx).^2+(y-centery).^2);
%     h=z;

    d_dist=10;
    h_dist=10;
    
    labelx_dist=100;
    labely_dist=100;

    d_range=0:d_dist:600;
    h_range=0:h_dist:1000;
    

   % hist_h=hist(h,h_range);
   % hist_d=hist(d_axis,d_range);

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
    
    
         
    data(:,5)= intens; % for including the solvatochromism
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
    
    



% %% Plotting the localization density wrt median of I_ratio
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
medint=medint';
    %% visualization
    figure_height_in_pixel=600;
    
    figure;
    set(gcf,'Position',[20 20 figure_height_in_pixel max(h_range)/max(d_range)*figure_height_in_pixel]);
    count=count'./repmat(V,size(count,2),1);
%     countbin = count;

    imagesc(count); % To plot wrt density values
    %imagesc(medint); % To plot wrt I_ratios
    
    axis equal % suggestion by Eric to prevent elonagation in z
    
   %paramters.count = count;
   
    set(gca,'Ydir','Normal');
    set(gca,'FontSize',22);
    xlim([0.5 max(d_range/d_dist)+1]);
    ylim([0.5 max(h_range/h_dist)+1]);
%    title('\fontsize{70}\fontname{Arial} c','position',[-7.424 82.9050 0])
    xlabel('Dist. from cylinder axis / nm', 'FontSize',22);
    ylabel('Relative z / nm', 'FontSize',22);
    set(gca,'XTick',d_range/d_dist*labelx_dist/d_dist+0.5);
    set(gca,'XTickLabel',0:labelx_dist:max(d_range));
    set(gca,'YTick',h_range/h_dist*labely_dist/h_dist+0.5);
    set(gca,'YTickLabel',-500:labely_dist:500);
%     set(gcf,'color','w');%Ashvini test
%         set(gca,'YTick',-300:10:300);
%         set(gca,'YTickLabel',0:labely_dist:max(h_range));
hcb = colorbar;


if constant_colorbar == 1 % keeping the colorbar constant
    cmax = maximum_limcolor;
else 
    cmax = max(max(count));
    %cmax = max(max(medint));
end
%cmin= -(cmax/(interval+1));
%cmin = -(cmax/(interval+1))+0.4; % correction for colorbar to start the blue color at almost 0 
%caxis([cmin cmax]);

if LowerIntensity == 0 && UpperIntensity == 2000
    % colormap gray;
      [map] = ColormapAll(colormap); %For jet
       colormap(map);
%       
     % [map] = Colormapblue2red(colormap);% Use this for Radial solvato values
   %colormap(map); 

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
   % [map] = Colormapblue2red(colormap);% Use this for Radial solvato values
  % colormap(map); 
   [map] = ColormapAll(colormap);
   colormap(map);
end

 %Namemat = baseNamecsv;
  %save(Namemat, 'data','count', 'V');
   



   
 saveas(figure(1),[folder, baseNamecsv,'_2D-density-distr.tif']);
 saveas(figure(1),[folder, baseNamecsv,'_2D-density-distr.fig']);
    %set(gca,'XTickLabel', d_range-d_dist/2);
    %set(gca,'YTickLabel', h_range+h_dist/2);

    %% end Ashvini code ends
%% 
%2-dimensional Matrices to calculate the distances between individual
%localizations
xmatrix1=ones(length(xvec),1)*xvec';
xmatrix2=xvec*ones(1,length(xvec));

ymatrix1=ones(length(yvec),1)*yvec';
ymatrix2=yvec*ones(1,length(yvec));

zmatrix1=ones(length(zvec),1)*zvec';
zmatrix2=zvec*ones(1,length(zvec));

%every distance between every localization is saved to distancematrix
%(Mirrored Matrix like a Table of city distances)
distancematrix=sqrt((xmatrix2-xmatrix1).^2+(ymatrix2-ymatrix1).^2+(zmatrix2-zmatrix1).^2);

%Getting a Histogram for the distribution of distances within a filtered
%range 
Distances=distancematrix(1:end);
Distances(Distances==0)=[];
Distances(Distances<10)=[];
Distances(Distances>=3000)=[];
%figure
%H = histogram(Distances,100); 

%In which radius around each localization should I search for neighbors?
%[nm]
filt_radius = 20; % [nm]

%This loop appends to the main matrix of all distances the indiviual Number
%of Next Neighbours within x nm. And additional the value of x for
%reproducibility and comparison.
for localiz = 1:length(distancematrix)
temp_auswahl = distancematrix(:,localiz);
temp_auswahl(temp_auswahl>=filt_radius)=[];
num_localiz = length(temp_auswahl);
all_localiz(localiz,6) = num_localiz;
all_localiz(localiz,7) = filt_radius;
end

figure;
colormap jet;
colorbar;
set(gcf,'Color','w')
axis equal;
axis square;
%This creates the Microgel Plot
%set(gcf, 'Color', [1,1,1]);
xmean = median(xvec);
ymean = median(yvec);
zmean = median(zvec);
%scatter3(xvec,yvec,zvec,10,all_localiz(:,6),'filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
%scatter3(xvec-(xmean),yvec-(ymean),-(zvec-(zmean)),10,all_localiz(:,4)/100,'filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6);%Ashvini plot wrt intensity ratio
scatter3(xvec,yvec,zvec,10,all_localiz(:,4)/100,'filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6);%Ashvini plot wrt intensity ratio
hold on
% zxy = zeros(length(xvec),1);

xlabel('x [nm]');
ylabel('y [nm]');
zlabel('z [nm]');
% xlim([-300 300])
% ylim([-300 300])
% zlim([-300 300])
% maketics(150);
az = 0;
el = 90;
view(az, el);
axis square
axis equal




%This will create the title of the colorbar, indicating also the radius for
%checking for neighbors
hcb = colorbar;
%txt = (['N_{loc}(',num2str(filt_radius),' nm)']);
txt = (['\fontname{calibri} \it I \rm\fontname{calibri}_{ratio}(L/S)']);% Ashvini to plot wrt I_ratio
title(hcb,txt);




% scatter3(0,0,0,50,'x','MarkerEdgeColor','k','LineWidth',10)
% set(gca,'linewidth',2)
% set(gca,'FontSize',18) 
% createsc_rh(0,0,Rh)
% load ('colormap_32.mat');
% colormap(c);
% %colormap(imcomplement(c));
% set(gca,'linewidth',2)
% set(gcf,'Color','k')
% hold off
% set(gcf,'Color','k')
% ax = gca % Get handle to current axes.
% ax.XColor = 'w'; % Red
% ax.YColor = 'w'; % Blue
% ax.ZColor = 'w'; % Blue
% set(gca,'Color','k')
% set(hcb,'Color','w')
% box on
% grid off

%If this part is uncommented, it will produce a short Video in *.mp4 File
%format. The numbers represent the track of the camera. (viewing positions
%and angles)
% OptionZ.FrameRate=120;OptionZ.Duration=5.5;OptionZ.Periodic=true;
% CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'Microgel',OptionZ)



% %Finding the "middle" of each microgel (careful, its just the median!)
xmean = median(xvec);
ymean = median(yvec);
zmean = median(zvec);
climits=([LowerLimit UpperLimit]);%Ashvini
caxis(climits); %Ashvini
scatter3(xmean,ymean,zmean,100,'x','MarkerEdgeColor','k','LineWidth',2);
hold off


%% EXPORT CROSSECTION PLOTS
hold on
xlim([xmean-Rh/2-50 xmean+Rh/2+50]);
ylim([ymean-Rh/2-50 ymean+Rh/2+50]); 
zlim([zmean-15 zmean+15]);
saveas(figure(2),[folder, baseNamecsv,'_crossection.tif']);
txt = ([folder, baseNamecsv,'_crossection.png']);
screen2png(txt);
%% To save the individual cross section -Ashvini 
Ind_crossecn = all_localiz;
Ind_crossecn(:,1) = xvec; 
Ind_crossecn(:,2) = yvec;
Ind_crossecn(:,3) = zvec;
Ind_crossecn(Ind_crossecn(:,3)>(zmean+15),:)=[];
Ind_crossecn(Ind_crossecn(:,3)<(zmean-15),:)=[];
%   figure
%   scatter3(Ind_crossecn(:,1)-xmean,Ind_crossecn(:,2)-ymean,Ind_crossecn(:,3)-zmean,10,Ind_crossecn(:,4)/100,'filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6);
Ind_crossecn_mean = Ind_crossecn;
Ind_crossecn_mean(:,1)= Ind_crossecn(:,1)-xmean;
Ind_crossecn_mean(:,2)= Ind_crossecn(:,2)-ymean;
Ind_crossecn_mean(:,3)= Ind_crossecn(:,3)-zmean;
%   saveas(figure(2),[folder, baseNamecsv,'_crossection.tif']);
%  txt = ([folder, baseNamecsv,'_crossection.png']);
%  screen2png(txt);

%% HISTOGRAM OF DISTANCES WITHIN THE MICROGEL
for localiz = 1:length(xvec)
all_localiz(localiz,8)=sqrt((xvec(localiz)-xmean).^2+(yvec(localiz)-ymean).^2+(zvec(localiz)-zmean).^2);
end
hist_binwidth = 10;
%hist_binwidth = 5;% Ashvini test
hist_binedges = [0:hist_binwidth:500];
figure;
H = histogram(all_localiz(:,8),hist_binedges);
hist_values = H.Values;
%Defines Y as 
Y = (4/3)*pi*(((H.BinEdges(2:end)).^3)-((H.BinEdges(1:end-1)).^3));
rdf = hist_values ./ Y;
hist_binmids = H.BinEdges(1:end-1) + hist_binwidth/2;
figure %===============================================HISTOGRAM: RADIAL DISTRIBUTION OF LOCALIZATIONS
bar(hist_binmids,rdf);
xlabel('r [nm]');
ylabel('Normalized Occurance g(r)');
ergebnis = [rdf];
ergebnis = ergebnis';
txt = ([folder, baseNamecsv,'_bins.csv']);
csvwrite (txt, hist_binmids);
 climits=([LowerLimit UpperLimit]);
 %climits=([0 2]);%Ashvini
plotter4(xvec,yvec,zvec,all_localiz,filt_radius,climits,folder,baseNamecsv,xmean,ymean,zmean,Rh);
end

