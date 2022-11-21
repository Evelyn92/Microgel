function density_vector = VisualizeDensityPlot(density_vector, save_add_path, processed_path, index, plot_config)
        Visible = false;
        h_bias = plot_config{1};
        d_max = plot_config{2};
        h_max = plot_config{3};
        className = plot_config{4};
        
        add = strcat(save_add_path,string(index),'.ply');

        test_local = pcread(add);
        all_localiz  = test_local.Location();

        zvec=all_localiz(:,3);
        xcoord=all_localiz(:,1);
        deltax = xcoord-median(xcoord);
        ycoord=all_localiz(:,2);
        deltay = ycoord-median(ycoord);
        %data = all_localiz;
        data = [deltax, deltay, zvec];

        %Visualization the pd
         %data(:,3) = 1-data(:,3); %flip z dimension
         cloud = pointCloud(data);
         g=figure("Visible", 'off');
         pcshow(cloud);
         if strcmp(className, 'FOCTS')
           % title('\fontsize{30}\fontname{Calibri} FOCTS','position',[-7.424 82.9050 0])
           title('\fontsize{30}\fontname{Calibri} FOCTS-'+string(index))
        elseif strcmp(className, 'ODS')
            title('\fontsize{30}\fontname{Calibri} ODS'+string(index))
        else
            title('\fontsize{30}\fontname{Calibri} PEG'+string(index))
        end
         
         fname_pc = strcat(processed_path,className,'/',className,'-',string(index),'-pc.png');
         saveas(g,fname_pc)

         %calculations and preparation of histogram ranges
        d_axis = sqrt(deltax.^2+deltay.^2);
        h = (zvec-median(zvec))+h_bias;

        clear count V d_range h_range;
        d_dist=10;
        h_dist=10;

        labelx_dist=100;
        labely_dist=100;

        d_range=0:d_dist:d_max;
        h_range=0:h_dist:h_max;

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


        if strcmp(className, 'FOCTS')
           % title('\fontsize{30}\fontname{Calibri} FOCTS','position',[-7.424 82.9050 0])
           title('\fontsize{30}\fontname{Calibri} FOCTS-'+string(index))
        elseif strcmp(className, 'ODS')
            title('\fontsize{30}\fontname{Calibri} ODS'+string(index))
        else
            title('\fontsize{30}\fontname{Calibri} PEG'+string(index))
        end

        
        xlabel('Dist. from symmetry axis / nm','FontSize',22);
        ylabel('Relative z / nm','FontSize',22);
        set(gca,'XTick',d_range/d_dist*labelx_dist/d_dist+0.5);
        set(gca,'XTickLabel',0:labelx_dist:max(d_range));
        set(gca,'YTick',h_range/h_dist*labely_dist/h_dist+0.5);
        set(gca,'YTickLabel',-max(h_range)/2:labely_dist:max(h_range)/2);

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

        LowerIntensity = 0;
        UpperIntensity = 1;
        if LowerIntensity == 0 && UpperIntensity == 2000      
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

        fname_2d = strcat(processed_path,className,'/',className,'-',string(index),'-density.png');
        saveas(gcf,fname_2d)
  
end