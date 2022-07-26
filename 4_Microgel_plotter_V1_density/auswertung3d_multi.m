%Analysis combined for several Microgels, given as *.3d data
 interval = 1126; % interval is the amount of colors inserted for the color map
%interval = 257;
%maximum_limcolor= 0.8; %The maximum limit you want to set with constant colorbar
%maximum_limcolor= 1*10^-3; %for density
constant_colorbar=0; % if this is 1 then constant maximum value for the colorbar 

%% Analysis combined for several Microgels, given as *.3d data
prompt = 'What do you expect to be Rh? (nm)';
Rh = input(prompt);
prompt = 'Do you want to analyze more than 1 microgel? (1/0)';
morethanone = input(prompt);
if morethanone == 1  
   disp('Okay, I will analyze more than 1 microgel.');
   prompt = 'Largest increment for filename?';
   numberruns = input(prompt);
   disp(['Ok. Analyzing ' , num2str(numberruns),' microgels. Give the first file.']);
%    [baseNamecsv, folder] = uigetfile('*.3d','Select the positioning 3d-File','D:\sciebo\AK_Woell_ERIC\12_3D\matlab_plots\allmicrogels');   
    [baseNamecsv, folder] = uigetfile('*.3d','Select the positioning 3d-File','C:\Users\Lukas\Dropbox\Studium\Master\Fluorescence FA\Bilder\Final Talk\Intensity Ratios\Contrast Test');   
   fullFileNamecsv = fullfile(folder, baseNamecsv);
[pathstr,filename,ext] = fileparts(fullFileNamecsv);
    radlocdistr = [];
    All_crossectn = []; 
    All_crossectn=cell(numberruns,1); % declaring Cell array for storing all the individual crosssections of the microgels
    All_countbin=zeros(80,40,numberruns);
    All_medint=zeros(80,40,numberruns);
   for n=1:numberruns
      newNamecsv = [num2str(n),ext];
      [labelx_dist, labely_dist, d_dist, h_dist, d_range, h_range,count,Ind_crossecn,Ind_crossecn_mean,ergebnis,hist_binmids,hist_binwidth,LowerIntensity,UpperIntensity,medint] = auswertung3d(newNamecsv, folder,Rh); %Passing the output from auswertung3d to the main function
      radlocdistr(:,n) = ergebnis; 
      All_crossectn{n,1} = Ind_crossecn_mean; % saving the individual cross section values for each microgel
      All_countbin(1:end,1:end,n) = count;
      All_medint(1:end,1:end,n) = medint;
%       All_crossectn{n} = Ind_crossecn_mean;
   end
   count_avg = mean(All_countbin,3);
   %medint_med = nanmedian(All_medint,3);
   medint_avg = nanmean(All_medint,3);
   %% visualization average density profile
   figure_height_in_pixel=600;
    
    figure;
    set(gcf,'Position',[20 20 figure_height_in_pixel max(h_range)/max(d_range)*figure_height_in_pixel]);
%     count=count'./repmat(V,size(count,2),1);
%     countbin = count;
    imagesc(count_avg);
    %imagesc(medint_med);
   % imagesc(medint_avg);
    axis equal % suggestion by Eric to prevent elonagation in z
    set(gca,'Ydir','Normal');
    set(gca,'FontSize',22);
    xlim([0.5 max(d_range/d_dist)+1]);
    ylim([0.5 max(h_range/h_dist)+1]);
    title('\fontsize{70}\fontname{Arial} c','position',[-7.424 82.9050 0]) % change text for every temperature
    
    xlabel('Dist. from cylinder axis / nm','FontSize',22);
    ylabel('Relative z / nm','FontSize',22);
    set(gca,'XTick',d_range/d_dist*labelx_dist/d_dist+0.5);
    set(gca,'XTickLabel',0:labelx_dist:max(d_range));
    set(gca,'YTick',h_range/h_dist*labely_dist/h_dist+0.5);
    set(gca,'YTickLabel',-400:labely_dist:400);
%     set(gcf,'color','w');%Ashvini test
%         set(gca,'YTick',-300:10:300);
%         set(gca,'YTickLabel',0:labely_dist:max(h_range));
hcb = colorbar;

if constant_colorbar == 1 % keeping the colorbar constant
    cmax = maximum_limcolor;
else 
    cmax = max(max(count_avg));
  % cmax = max(max(medint_med)); % for median values of I_ratio
   %cmax = max(max(medint_avg)); % for mean values of I_ratio
end
%cmin=0.4;
cmin= -(cmax/(interval+1)); % correction for colorbar to start the blue color at almost 0 
caxis([cmin cmax]);

if LowerIntensity == 0 && UpperIntensity == 2000
    % colormap gray;
    
%        [map] = Colormapblue2red(colormap); %Change here
%    colormap(map); 
   caxis([cmin cmax]);
      [map] = ColormapAll(colormap);
       colormap(map);

elseif LowerIntensity == 0 && UpperIntensity == 200
     [map] = ColormapAll(colormap);
    colormap(map);
        colormap gray;
        
        

else
  % colormap hot;
  %  [map] = Colormapblue2red(colormap); %Change here
  % colormap(map); 
   [map] = ColormapAll(colormap);
    colormap(map);
end
saveas(figure(6),[folder, baseNamecsv,'_Average_2D-density-distr.tif']);
 saveas(figure(6),[folder, baseNamecsv,'_Average_2D-density-distr.fig']);
%% end
%%



    
   hold off
   figure
   for n=1:numberruns
   radius = [1:length(ergebnis)];
   hold on
   plot(radius*hist_binwidth,radlocdistr(:,n),'x');
   end
   radlocdistr(:,n+1)=mean(radlocdistr,2);
   plot(radius*hist_binwidth,radlocdistr(:,n+1),'--');
  hold on
   li = linspace(0,max(radlocdistr(1,:)),10);
   plot(repmat(Rh,[1 10]),li,'--','MarkerFaceColor','red')
   
  
    



    
else
  disp('Ok. Analyzing 1 microgel.');
  [baseNamecsv, folder] = uigetfile('*.3d','Select the positioning 3d-File','D:\sciebo\AK_Woell_ERIC\12_3D\matlab_plots\allmicrogels');
[ergebnis,hist_binmids,hist_binwidth,Ind_crossecn_mean] = auswertung3d(baseNamecsv, folder);
fullFileNamecsv = fullfile(folder, baseNamecsv);
[pathstr,filename,ext] = fileparts(fullFileNamecsv);
end





