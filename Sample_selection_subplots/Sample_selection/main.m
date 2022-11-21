%%
clear all;
clc;
close all;
%% Configurations 
%Change the path to your folder
path = "C:\Users\evyji\Desktop\tmp\Sample_selection_subplots\Sample_selection";
cd (path) %The path to the main.m folder

%The max num of samples that will be extracted from the slides for each class.
max_samples = 15; %you can set 78 for ODS
%process the candidate samples from min_inx=1 to max_inx in each class 
max_inx = 10;

% The number of selected samples
selected_sample_num = 4;
selected_sample_num = min(max_inx, selected_sample_num);% less than the max number of processed samples

%% 1) read the slide files(***.3d) from the raw_data folder
addpath("./raw_dataset/")% The folder which contains LHA_(CLASSNAME)_NIPMAM_ALL.3d files
mkdir("./samples/")% Generated for the extracted point clouds
addpath("./samples/")

slides = dir('./raw_dataset/*.3d');
class_num = length(slides);
if class_num == 0
warning('There is no slide files in the folder.')
end
disp("The classes of microgel we have in the folder: ")
disp(slides(1).folder)

classes_nm = cell(class_num,1);
cloud_add = cell(class_num,1);
class_name = cell(class_num,1);
save_add = cell(class_num,1);

for idx = 1:class_num
    classes_nm{idx} = slides(idx).name;
    cloud_add{idx} = slides(idx).folder;
    cloud_add{idx} = strcat(cloud_add{idx},'\',classes_nm{idx});
    class_nm = split(classes_nm{idx}, "_");
    class_name{idx} = string(class_nm(2));
    disp(class_name{idx});
    save_add{idx} = strcat('./samples/', class_name{idx}, '/');
end
%% 2) read the entire slide: (the three .3d files contain the data from the Nano Letters paper)

ExtractSamples(class_name, cloud_add, save_add, max_samples)

disp("The extracted samples are saved in the path:")
for idx = 1:class_num
disp(save_add{idx})
end
%% 3）Sample selection
% the path to the processed_samples
processed_path = 'processed_samples/';
min_inx = 1;
SampleSelection(save_add, processed_path,class_name,min_inx, max_inx, selected_sample_num)


close all;

%% 4) Display the selected samples and discareded ones
SelectedVsDiscarded(path);
close all;
