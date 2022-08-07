
clc
clear all;
% close all
% 21 33 35 38 43 53

%experiment='_hi-lo_';
experiment='_hi-lo_';
if strcmp(experiment,'_hi-lo_')
    cd D:\NutCloud\我的坚果云\RWTH-Study\Hiwi\07.08\rez_hi_lo
else
    cd D:\NutCloud\我的坚果云\RWTH-Study\Hiwi\07.08\result_drop35
end



% micro_index = '12';
type = 'gg';%or "_go"
%sample_index = 1;
tmp_index = 1;
if tmp_index ==1
    temperature = '21';
elseif tmp_index ==2
    temperature = '33';
elseif tmp_index ==3
    temperature = '35';
elseif tmp_index ==4
    temperature = '38';
elseif tmp_index ==5
    temperature = '43';
else
    temperature = '53';
end


add1 = strcat("./results/results_mtemp/",'sample_',string(1),type,string(tmp_index),".ply");
add2 = strcat("./results/results_mtemp/",'sample_',string(2),type,string(tmp_index),".ply");
add3 = strcat("./results/results_mtemp/",'sample_',string(3),type,string(tmp_index),".ply");
add4 = strcat("./results/results_mtemp/",'sample_',string(4),type,string(tmp_index),".ply");
% add = strcat("./results/results_mtemp/",'gg1',".ply");
% all_localiz_test  = dlmread(add);
% add_info_test = strcat("./PAINT_DiffTemp/toy_example/results/results_mtemp/",'info1',".ply");
% info_test = dlmread(add_info_test);
% all_localiz_test  = all_localiz_test  .* (info_test(2,:)-info_test(3,:)) + info_test(3,:)+ info_test(1,:);
% all_localiz  = dlmread(add);

all_localiz1  = dlmread(add1);
all_localiz2  = dlmread(add2);
all_localiz3  = dlmread(add3);
all_localiz4  = dlmread(add4);

add_info1 = strcat("./results/results_mtemp/",'sample_',string(1),'info',string(tmp_index),".ply");
add_info2 = strcat("./results/results_mtemp/",'sample_',string(2),'info',string(tmp_index),".ply");
add_info3 = strcat("./results/results_mtemp/",'sample_',string(3),'info',string(tmp_index),".ply");
add_info4 = strcat("./results/results_mtemp/",'sample_',string(4),'info',string(tmp_index),".ply");

%
info1  = dlmread(add_info1);
info2  = dlmread(add_info2);
info3  = dlmread(add_info3);
info4  = dlmread(add_info4);


all_localiz1  = all_localiz1  .* (info1(2,:)-info1(3,:)) + info1(3,:)+ info1(1,:);
all_localiz2  = all_localiz2  .* (info2(2,:)-info2(3,:)) + info2(3,:)+ info2(1,:);
all_localiz3  = all_localiz3  .* (info3(2,:)-info3(3,:)) + info3(3,:)+ info3(1,:);
all_localiz4  = all_localiz4  .* (info4(2,:)-info4(3,:)) + info4(3,:)+ info4(1,:);

all_localiz = [all_localiz1;all_localiz2;all_localiz3;all_localiz4];

% all_localiz  = all_localiz  .* (info(2,:)-info(3,:)) + info(3,:)+ info(1,:);
% intens=all_localiz(:,4);
% Yiwei: Find what is the outlier
UpperSolva = 400;
%all_localiz = all_localiz(intens<UpperSolva,:);% if a data point has the Solva>400, throw it away
% all_localiz(intens>UpperSolva,4) = UpperSolva ;
zvec1=all_localiz1(:,3);
% intens1=all_localiz1(:,4);
xcoord1=all_localiz1(:,1);
deltax1 = xcoord1-median(xcoord1);
ycoord1=all_localiz1(:,2);
deltay1 = ycoord1-median(ycoord1);

zvec2=all_localiz2(:,3);
% intens2=all_localiz2(:,4);
xcoord2=all_localiz2(:,1);
deltax2 = xcoord2-median(xcoord2);
ycoord2=all_localiz2(:,2);
deltay2 = ycoord2-median(ycoord2);

zvec3=all_localiz3(:,3);
% intens3=all_localiz3(:,4);
xcoord3=all_localiz3(:,1);
deltax3 = xcoord3-median(xcoord3);
ycoord3=all_localiz3(:,2);
deltay3 = ycoord3-median(ycoord3);

zvec4=all_localiz4(:,3);
% intens4=all_localiz4(:,4);
xcoord4=all_localiz4(:,1);
deltax4 = xcoord4-median(xcoord4);
ycoord4=all_localiz4(:,2);
deltay4 = ycoord4-median(ycoord4);

deltax = [deltax1; deltax2; deltax3; deltax4];
deltay = [deltay1; deltay2; deltay3; deltay4];
zvec = [zvec1; zvec2; zvec3; zvec4];
% visualize one microgel in 3D space
% %data = all_localiz;
% data = [deltax deltay zvec];
% if(size(data,2)>3) 
%     data = data(:,1:3);
% end
% 
% data(:,3) = 1-data(:,3); %flip z dimension
% 
% cloud = pointCloud(data);
% pcshow(cloud);

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
    if strcmp(temperature,'38')
         title('\fontsize{70}\fontname{Calibri} g')
    elseif strcmp(temperature,'21')
        title('\fontsize{70}\fontname{Calibri} d','position',[0 78.9050 0])
% %         title('\fontsize{70}\fontname{Calibri} d','position',[-7.424 82.9050 0])

    elseif strcmp(temperature,'33')
        title('\fontsize{70}\fontname{Calibri} e')
    elseif strcmp(temperature,'35')
        title('\fontsize{70}\fontname{Calibri} f')
    elseif strcmp(temperature,'43')
        title('\fontsize{70}\fontname{Calibri} h')
    else
        title('\fontsize{70}\fontname{Calibri} i')
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
% maximum_limcolor = 6e-03;
maximum_limcolor = 0.01;
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
if strcmp(experiment,'_hi-lo_')
    cd D:\NutCloud\我的坚果云\RWTH-Study\Hiwi\07.08\visualization_hi_lo
else
    cd D:\NutCloud\我的坚果云\RWTH-Study\Hiwi\07.08\visualization_drop35
end
fname_2d = strcat(temperature,"/average_2d-density-distr",experiment,temperature,'C','.png');
% fname_2d = strcat("./PAINT_DiffTemp/toy_example/results/yiwei_fig/new_try",type,string(tmp_index),".png");
saveas(gcf,fname_2d)
    %set(gca,'XTickLabel', d_range-d_dist/2);
    %set(gca,'YTickLabel', h_range+h_dist/2);

    %% end Ashvini code ends