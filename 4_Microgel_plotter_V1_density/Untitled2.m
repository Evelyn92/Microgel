
clc
clear all
close all

add = "./results/results_mtemp/go1.ply";
all_localiz  = dlmread(add);
pcshow(all_localiz(:,1:3));


figure()
add = "./results/results_mtemp/gg1.ply";
all_localiz  = dlmread(add);
pcshow(all_localiz(:,1:3));



