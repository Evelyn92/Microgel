im1 = imread('D:\Ashvini\PAINT\PAINT-set-2\2018-12-04-21C\Region2\Movie1\for_ugel_plotter\40to80_I_ratio_mean_whole range_imp\1.3d_Average_2D-density-distr.tif');
im2 = imread('D:\Ashvini\PAINT\PAINT-set-2\2018-12-05-33C\Region1\Movie1\for-microgel-plotter\40to80_whole_range_imp\1.3d_Average_2D-density-distr.tif');
figure
montage({im1, im2, 'cameraman.tif'})