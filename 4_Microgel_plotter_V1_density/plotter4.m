%plotter4
function [bla] = plotter4(xvec,yvec,zvec,all_localiz,filt_radius,climits,folder,filename,xmean,ymean,zmean,Rh)
interval = 257;
% offsetx = Rh - xmean + 5; %Ashvini und Eric
% offsety = Rh - ymean + 5;  %Ashvini und Eric
cmax0 = 1.2;
cmin0 = 0.2;
figure
set(gcf, 'Color', [1,1,1]);
view_xyz = subplot(2,2,1);
%view_xyz = plot;
%scatter3(xvec,yvec,zvec,16.5,all_localiz(:,6),'filled')%,'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
% scatter3(xvec+offsetx,yvec+offsety,zvec,16.5,all_localiz(:,4),'filled')%,'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6) %Ashvini plot wrt I_ratio
scatter3(xvec-(xmean),yvec-(ymean),(zvec-(zmean)),10,all_localiz(:,4)/100,'filled')% For shifting the center of the plot to center
hold on
%scatter3(xmean+offsetx,ymean+offsety,zmean,16.50,'x','MarkerEdgeColor','k','LineWidth',2)%Ashvini plot wrt I_ratio
scatter3(0,0,0,100,'x','MarkerEdgeColor','k','LineWidth',2)
hold on
az = -45;
el = 30;
lw = 1.5;
scaling_x = 0.98;
scaling_y = 0.96;
scaling_z = 0.96;
% scaling_x = 0.5;
% scaling_y = 0.5;
% scaling_z = 0.5;
view(az, el);
set(gca,'linewidth',lw)
axis square
axis equal
max_exp = max([xmean ymean zmean]);
% xlim([0 2*max_exp*scaling_x])
% ylim([0 2*max_exp*scaling_y])
% zlim([0 2*max_exp*scaling_z])
%%  Ashvini test
xlim([-300 300])
ylim([-300 300])
zlim([-300 300])
% Ashvini test end
%maketics(300);
maketics(150);%Ashvini 

%colormap jet
% caxis(climits)
[map] = Colormapblue2red(colormap);% Use this for Radial solvato values from blue to red
   colormap(map);
   cmax = cmax0;
cmin = -(cmax/(interval))+cmin0; % correction for colorbar to start the blue color at almost 0 
caxis([cmin cmax]);


set(view_xyz, 'FontName', 'Arial')
set(view_xyz,'FontSize',10,'FontWeight','normal') 
%txt = ('\fontname{Times} \it x\rm\fontname{Arial} / nm');
txt = ('x / nm');
xlabel(txt)
%txt = ('\fontname{Times} \it y\rm\fontname{Arial} / nm');
txt = ('y / nm');
ylabel(txt)
%txt = ('\fontname{Times} \it z\rm\fontname{Arial} / nm');
txt = ('z / nm');
zlabel(txt)


view_xy = subplot(2,2,2);
%scatter3(xvec,yvec,zvec,16.5,all_localiz(:,6),'filled')%,'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)


% scatter3(xvec+offsetx,yvec+offsety,zvec,16.5,all_localiz(:,4),'filled')%,'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)%Ashvini plot wrt I_ratio, for the points
scatter3(xvec-(xmean),yvec-(ymean),(zvec-(zmean)),10,all_localiz(:,4)/100,'filled')
%plot(xvec,yvec,'.')
%scatter(xvec,yvec,1,'MarkerFaceColor',all_localiz(:,6));
hold on
% scatter3(xmean+offsetx,ymean+offsety,600,16.50,'x','MarkerEdgeColor','k','LineWidth',2)% for the cross Ashvini
scatter3(0,0,0,100,'x','MarkerEdgeColor','k','LineWidth',2)
%c1 = min(all_localiz(:,6));
%c1 = min(all_localiz(:,4)); %Ashvini 0
c1 = 0;
%c2 = max(all_localiz(:,6));
%c2 = max(all_localiz(:,4));%Ashvini
c2 = 40;
%c2=(all_localiz(:,4)==40);%Ashvini
cm0 = linspace(c1,c2,100)'/c2;  
cm = [cm0,ones(100,1),ones(100,1)]; %die eigentliche Colormap
%colormap(cm); uncomment later test
[map] = Colormapblue2red(colormap);% Use this for Radial solvato values from blue to red
   colormap(map);
   cmax = cmax0;
cmin = -(cmax/(interval+50))+cmin0; % correction for colorbar to start the blue color at almost 0 
caxis([cmin cmax]);
    %create_rh(xmean+offsetx,ymean+offsety,Rh) 
    create_rh(0,0,Rh) % To create circle around the microgel
az = 0;
el = 90;
view(az, el);
set(gca,'linewidth',lw)
axis square
axis equal

% xlim([0 2*max_exp*scaling_x])
% ylim([0 2*max_exp*scaling_y])
% zlim([0 2*max_exp*scaling_z])
%%  Ashvini test
%%  Ashvini test
xlim([-300 300])
ylim([-300 300])
zlim([-300 300])
% Ashvini test end
maketics(150);
%colormap jet

[map] = Colormapblue2red(colormap);% Use this for Radial solvato values from blue to red
   colormap(map);
cmax = cmax0;
cmin = -(cmax/(interval))+cmin0; % correction for colorbar to start the blue color at almost 0 
caxis([cmin cmax]);

%caxis(climits)
%set(view_xy, 'Position', [0.575 0.65 0.30 0.30])
%set(view_xy, 'Position', [0.500 0.65 0.30 0.30])
set(view_xy, 'Position', [0.500 0.65 0.30 0.30])

set(gca, 'FontName', 'Arial')
set(gca,'FontSize',10,'FontWeight','normal') 
%txt = ('\fontname{Times} \it x\rm\fontname{Arial} / nm');
txt = ('x / nm');
xlabel(txt)
%txt = ('\fontname{Times} \it y\rm\fontname{Arial} / nm');
txt = ('y / nm');
ylabel(txt)
%txt = ('\fontname{Times} \it z\rm\fontname{Arial} / nm');
txt = ('z / nm');
zlabel(txt)
box on
grid off

view_zx = subplot(2,2,3);
%scatter3(xvec,yvec,zvec,16.5,all_localiz(:,6),'filled')%,'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
scatter3(xvec-(xmean),yvec-(ymean),(zvec-(zmean)),10,all_localiz(:,4)/100,'filled')%,'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)%%Ashvini plot wrt I_ratio
hold on
scatter3(0,0,0,100,'x','MarkerEdgeColor','k','LineWidth',2)
az = 0;
el = 0;
view(az, el);
set(gca,'linewidth',lw)
axis square
axis equal
% xlim([0 2*max_exp*scaling_x])
% ylim([0 2*max_exp*scaling_y])
% zlim([0 2*max_exp*scaling_z])
%%  Ashvini test
%%  Ashvini test
xlim([-300 300])
ylim([-300 300])
zlim([-300 300])
% Ashvini test end
% Ashvini test end
maketics(150);

%colormap jet

[map] = Colormapblue2red(colormap);% Use this for Radial solvato values from blue to red
   colormap(map);
  cmax = cmax0;
cmin = -(cmax/(interval))+cmin0; % correction for colorbar to start the blue color at almost 0 
caxis([cmin cmax]); 
%caxis(climits)
set(view_zx, 'Position', [0.15 0.15 0.30 0.30])
set(gca, 'FontName', 'Arial')
set(gca,'FontSize',10,'FontWeight','normal') 
%txt = ('\fontname{Times} \it x\rm\fontname{Arial} / nm');
txt = ('x / nm');
xlabel(txt)
%txt = ('\fontname{Times} \it y\rm\fontname{Arial} / nm');
txt = ('y / nm');
ylabel(txt)
%txt = ('\fontname{Times} \it z\rm\fontname{Arial} / nm');
txt = ('z / nm');
zlabel(txt)
grid off
box on

view_zy = subplot(2,2,4);
%plot(yvec,zvec,'.')
%scatter3(xvec,yvec,zvec,16.5,all_localiz(:,6),'filled')%,'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
%scatter3(xvec+offsetx,yvec+offsety,zvec,16.5,all_localiz(:,4),'filled')%,'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)%Ashvini to plot colormap wrt I_ratio
scatter3(xvec-(xmean),yvec-(ymean),(zvec-(zmean)),16.5,all_localiz(:,4)/100,'filled')%,'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)%Ashvini to plot colormap wrt I_ratio

hold on
%scatter3(600,ymean+offsety,zmean,16.50,'x','MarkerEdgeColor','k','LineWidth',2)%Ashvini plot wrt I_ratio
scatter3(0,0,0,16.50,'x','MarkerEdgeColor','k','LineWidth',2)%Ashvini plot wrt I_ratio
hold on
az = 90;
el = 0;
view(az, el);
set(gca,'linewidth',lw)
axis square
axis equal
% xlim([0 2*max_exp*scaling_x])
% ylim([0 2*max_exp*scaling_y])
% zlim([0 2*max_exp*scaling_z])
%%  Ashvini test
xlim([-300 300])
ylim([-300 300])
zlim([-300 300])
% Ashvini test end

maketics(150);
%colormap jet

[map] = Colormapblue2red(colormap);% Use this for Radial solvato values from blue to red
   colormap(map);
  cmax = cmax0;
cmin = -(cmax/(interval+50))+cmin0; % correction for colorbar to start the blue color at almost 0 
caxis([cmin cmax]); 
%caxis(climits)
%set(view_zy, 'Position', [0.575 0.15 0.30 0.30])

set(view_zy, 'Position', [0.500 0.15 0.30 0.30])
set(gca, 'FontName', 'Arial')
set(gca,'FontSize',10,'FontWeight','normal') 
%txt = ('\fontname{Times} \it x\rm\fontname{Arial} / nm');
txt = ('x / nm');
xlabel(txt)
%txt = ('\fontname{Times} \it y\rm\fontname{Arial} / nm');
txt = ('y / nm');
ylabel(txt)
%txt = ('\fontname{Times} \it z\rm\fontname{Arial} / nm');
txt = ('z / nm');
zlabel(txt)

grid off
box on
%colormap(map); %added by Eric and Ashvini

[map] = Colormapblue2red(colormap);% Use this for Radial solvato values from blue to red
   colormap(map);
   cmax = cmax0;
cmin = -(cmax/(interval+50))+cmin0; % correction for colorbar to start the blue color at almost 0 
caxis([cmin cmax]);
c = colorbar('vert');
%set(c,'Position',[0.91 0.1 0.1 0.80])
set(c,'Position',[0.85 0.1 0.1 0.80])
set(c,'LineWidth',lw)
ax = gca;
axpos = ax.Position;
cpos = c.Position;
cpos(3) = 0.2*cpos(3);
c.Position = cpos;
ax.Position = axpos;
set(gca, 'FontName', 'Arial','FontSize',11,'FontWeight','normal')
%txt = (['\fontname{Times} \it N \rm\fontname{Arial}_{loc}(',num2str(filt_radius),' nm)']);
txt = (['I_{ratio}(L/S)']);%Ashvini I_ratio plot
title(c,txt)
%colormap jet

[map] = Colormapblue2red(colormap);% Use this for Radial solvato values from blue to red
   colormap(map);
   cmax = cmax0;
cmin = -(cmax/(interval))+cmin0; % correction for colorbar to start the blue color at almost 0 
caxis([cmin cmax]);
% print(filename,'-dpng','-r600')
%txt = ([folder, filename,'_4plot.png']);

% txt = ([folder, filename,'_Iratio.tif']);%Ashvini plot wrt I_ratio
% screen2tif(txt)
 set(gcf, 'Units', 'Inches', 'Position', [0, 0, 6.9, 5.10], 'PaperUnits', 'Inches', 'PaperSize', [6.9, 5.10])
 saveas(figure(5),[folder, filename,'_Iratio.fig']);% Ashvini
 print(figure(5),[folder, filename,'_Iratio.png'],'-dpng','-r600')
 saveas(figure(5),[folder, filename,'_Iratio.png']);

end
