clear all;
clc;
cd D:\forClone\Hiwi\Microgel\Jia_sample_selection
className{1} = 'FOCTS';
className{2} = 'ODS';
className{3} = 'PEG';
save_add{1} = './dataset/focts/tr/'; 
save_add{2} = './dataset/ods/tr/';
save_add{3} = './dataset/peg/tr/';

className = className{3};
index = 7;

if strcmp(className, 'FOCTS')
    max_index = 7;
    save_add_path = save_add{1};
    h_bias = 200;
    d_max = 600;
    h_max = 400;
elseif strcmp(className, 'ODS')
    max_index = 10;
    save_add_path = save_add{2};
    h_bias = 400;
    d_max = 600;
    h_max = 800;
else 
    max_index = 8;
    save_add_path = save_add{3};
    h_bias = 400;
    d_max = 420;
    h_max = 800;
end


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

save_path = strcat(save_add_path, 'loc-', int2str(index),'.txt');
writematrix(data, save_path)
