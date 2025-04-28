clear; clc;
load('data\Mask.mat');
Mask_all = flip(Mask_all,3);
File_Path ='data\';
dataname = [File_Path,'\images.mat'];
slice_choose = [7];
result_path = [File_Path,'result_',num2str(slice_choose),'\']; 
    
[Signal_Cr, Signal_PCr, Signal_glycoNOE] = Fitting(File_Path,dataname,slice_choose,Mask_all,result_path);
