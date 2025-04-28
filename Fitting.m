%% demo: Rewrite the code for CEST analysis to make it suitable for data processing of 2D, 3D, and multiple B1 amplitudes
% including registration, ROI definition, B0 correction, B1 correction,
% partial-side and full z-spectrum fitting.
% All the result are saved in result_path

% Xi.Xu updated on 2023.8.4@siat, ShenZhen
% Thanks to Chongxue.Bie for the source code
% 2024.5.21 Xi.Xu test wuhao data and change the load information

% clear all; clc; close all;
function [Signal_Cr, Signal_PCr, Signal_glycoNOE] = Fitting(File_Path,dataname,slice_choose,Mask_all,result_path)
addpath('Tool_Fitting');
addpath('unity');

Para.Freq_ppm = -[-4:0.2:3]; % [4~-3] % wuhao
Para.Freq_ppm_use =  -[-4:0.1:3];% interpolate into 0.1ppm
Para.Side_ppm='all'; % 'pos_+ppm', 'neg_-ppm', 'all_+-ppm'
Para.B1_input = [0.3];
Para.B1_output = Para.B1_input;
Para.B1_correction = 'n'; % y: B1 correction, n: no correction
Para.Mask = 'overall'; % 'Multismall'£¬ revising
Para.ROI_number = 2;
Para.SelectedSlice = slice_choose; %1: 2D data, slice_index 2024.5.21


index_cest = [1];
index_b1map = []; % if there isn't b1map, keep it as []
index_T2map = []; % if there isn't T2map, keep it as []

% result_path = [File_Path,'result_nufft','\'];
if ~exist(result_path, 'dir')
    mkdir(result_path);
    disp('build "result" dir successfully');
else
    disp('dir "result" has exsited');
end

%% Data Preparation
% all the data are sorted in ascended order in slice dim
B1map_registraiton = [];
T2w_registraiton = [];
% CEST:x,y,z, ppm, b1
images = [];
filename_img=[result_path,'images.mat'];
if exist(filename_img)
    load(filename_img);
else
    
    %% nufft
    %     datapath = [File_Path,'\recon_cs.mat'];
    images = load(dataname);
%     images = images.recon;
    [x,y,nSlice,nFreq] = size(images);
    save(filename_img,'images');
    
end

% Registration the T2w/B1map images to CEST images
filename_T2=[result_path,'T2w_registraiton.mat'];
if exist(filename_T2)
    load(filename_T2);
end

filename_B1map=[result_path,'B1map_registraiton.mat'];
if exist(filename_B1map)
    load(filename_B1map);
end

filename=[result_path,'images_registraiton.mat'];
if exist(filename)
    load(filename);
else
    [images_registraiton, images]   = Registration_img(Para,images);  % not perfect£¬slow and inaccurate£¬ selected-slice
    %     images_registraiton = images;
    filename=[result_path,'images_registraiton.mat'];
    save(filename,'images_registraiton');
    
    if  ~isempty(index_T2map)
        img_T2w = readDICOM_2D([File_Path,FileNameList{index_T2map}]);
        T2w_registraiton = Registration_T2(Para,images,img_T2w);
        save([result_path,'T2w_registraiton.mat'],'T2w_registraiton');
    end
    
    % B1 map
    if ~isempty(index_b1map)
        img_B1map = readDICOM_2D([File_Path,FileNameList{index_b1map}]);
        B1map_registraiton = Registration_B1(Para,images,img_B1map);% need certificating
        save([result_path,'B1map_registraiton.mat'],'B1map_registraiton');
    end
end

%% Mask: x,y,z,tube
% filename=[result_path,'Mask.mat'];
% % if exist(filename)
% %     load('E:\2023\CEST\Code_use\Calf_3D_150spokes\result_LplusS_csm150_40_wuhao\Mask.mat');%wuhao
%       load([result_path,'\Mask.mat']);
% % slice_choose = 21-slice_choose;
% %       Mask_all(:,:,1) = Mask_all(:,:,slice_choose);
% % else
%     S0 = squeeze(images_registraiton(:,:,:,1,1)); % x,y,slice
% % %
%     for slice=Para.SelectedSlice
% %         Mask = DrawMask(S0(:,:,slice),Para.ROI_number);
%         Mask = DrawMask(S0(:,:,slice),1);
%         if Para.Mask == 'overall'
%             Mask_all(:,:,slice)  = get_mask(Mask);
% %         elseif Para.Mask == 'Multismall'
% %             Mask_all(:,:,slice,:) = Mask;
%         end
%     end
%     save(filename,'Mask_all');
% % end

%% B0 correction

filename=[result_path,'images_B0correction.mat'];
if exist(filename)
    load(filename);
else
    % B0map: x,y,slice,1,b1
    % delta_freq_map: x,y,slice,1,b1
    % images_B0correction: x,y,slice,freq, b1
    [images_B0correction, delta_freq_map,B0map] = B0correction(Para, images,Mask_all);
    
    save(filename,'images_B0correction','B0map','delta_freq_map');
end
Para.Freq_ppm = Para.Freq_ppm_use; % 20230925
%% B1 correction
filename=[result_path,'images_B1correction.mat'];
if exist(filename)
    load(filename);
else
    %     x,y,slice,freq, b1
    if ~isempty(index_b1map)&&Para.B1_correction == 'y'
        images_B1correction = B1correction(images_B0correction,Mask_all,Para,B1map_registraiton);
        save(filename,'images_B1correction');
    else
        images_B1correction = images_B0correction;
    end
end
%% CEST fitting
% initialize varables and reshape the data
[~,index]=sort(Para.Freq_ppm,'ascend');
freq_ppm = Para.Freq_ppm(:,index);
Para.Freq_ppm = freq_ppm;
images_B1correction = images_B1correction(:,:,:,index,:);  % x,y,slice,freq, b1
images_Fit = permute(images_B1correction,[1 2 5 4 3]);% x,y,b1,freq, slice


% +ppm
if strcmp(Para.Side_ppm, 'pos')
    [map_Cr,map_PCr,map_glycoNOE,res_B0map]=Fitting_PixelWise_down(images_Fit,Mask_all,B0map,Para,result_path);
elseif strcmp(Para.Side_ppm, 'neg')
    [map_Cr,map_PCr,map_glycoNOE,res_B0map]=Fitting_PixelWise_up(images_Fit,Mask_all,B0map,Para,result_path);
else %strcmp(Para.Side_ppm, 'all')
    [map_Cr,map_PCr,map_glycoNOE,res_B0map,Signal_Cr, Signal_PCr, Signal_glycoNOE]=Fitting_PixelWise_use(images_Fit,Mask_all,B0map,Para,result_path);
end

%% imshow and save maps
%ShowSave_map(map_Cr,map_PCr,map_glycoNOE,res_B0map,B1map_registraiton,T2w_registraiton,images_Fit, Mask_all,result_path,Para);
ShowSave_map_ROI(map_Cr,map_PCr,map_glycoNOE,res_B0map,B1map_registraiton,T2w_registraiton,images_Fit, Mask_all,result_path,Para,images_registraiton);
end
