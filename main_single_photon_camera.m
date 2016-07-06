%% main_single_photon_camera.m

% Implements framework presented in
% "PHOTON-EFFICIENT IMAGING USING A SINGLE-PHOTON CAMERA"
% D.Shin, F.Xu, D.Venkatraman, R.Lussana, F.Villa
% F.Zappa, V.Goyal, F.N.C.Wong, J.H.Shapiro
% in Nature Communications

% Uses variants of MATLAB packages available freely online:
% SPIRAL-TAP by Z.H.Harmany et al.
% OMP implementation by S.Becker

%% * Load single-photon camera data

clc; clear; close all;

spad_file_name = 'data/scene_man_flower_N100';

% frame parameters 
num_frames = 20;
load(spad_file_name);
[nr,nc] = size(D);
bin_ranges = 0:127;
% load functions in dir
addpath(genpath([pwd '/functions']));
C = zeros(nr,nc);
T = cell(nr,nc);
for i=1:nr
    for j=1:nc
        frames = F{i,j};
        dats = D{i,j};
        inds = find(frames<num_frames);
        if(~isempty(inds))
            dats_new = dats(inds);
            C(i,j) = length(dats_new);
            T{i,j} = [T{i,j} dats_new'];
        end
    end
end
sc = 1.5;
I_raw = C/(max(C(:)/sc));
I_raw(I_raw<0) = 0;
I_raw(I_raw>1) = 1;
load('data/data_supp'); % load supp data
t2cm = 100*((389e-12)*(3e8)/2);

fprintf('* Finished loading SPAD camera data\n');

%% * Perform image reconstruction from photon-sparse SPAD camera data

run_conventional() % run conventional filtered histogram LIDAR approach
run_proposed_step1() % run step 1 of our proposed framework
run_proposed_step2() % run step 2 of our proposed framework
run_proposed_step3() % run step 3 of our proposed framework
fprintf('* Finished image reconstruction from photon-sparse SPAD camera data!\n') 

%% * Plot results

range_C = [0,1];
range_D = [74,80];

vizsc = 12;
figure; 
subplot(2,3,1);
imshow(uint8(repmat(C_ML/vizsc*255,[1,1,3]))); axis image;
title({'intensity','(conventional)'})
subplot(2,3,4);
imshow(uint8(repmat(C_map/vizsc*255,[1,1,3]))); axis image;
title({'intensity','(proposed)'})
subplot(2,3,2);
imagesc(D_ML,range_D); axis image;
set(gca,'xtick',[],'ytick',[])
title({'depth','(conventional)'})
subplot(2,3,5);
imagesc(D_map,range_D); axis image; 
set(gca,'xtick',[],'ytick',[])
title({'depth','(proposed)'})
colormap('hot');

load('data/data_truth');
max_cm = 6;

D_ML_fin = D_ML; D_MAP_fin = D_map; 
D_truth_fin(~M_fin) = 0; D_ML_fin(~M_fin) = 0; D_MAP_fin(~M_fin) = 0;
D_ML_fin(isnan(D_ML_fin)) = 0;
E_ML = abs(D_truth_fin-D_ML_fin)*t2cm;
E_map = abs(D_truth_fin-D_MAP_fin)*t2cm;
err_ML = mean(abs(E_ML(E_ML>0)));
err_map = mean(abs(E_map(E_map>0)));

subplot(2,3,3);
imagesc(E_ML,[0,max_cm]); axis image;
set(gca,'xtick',[],'ytick',[])
title({'MAE(old depth)',['=' num2str(err_ML) ' cm']})
subplot(2,3,6);
imagesc(E_map,[0,max_cm]); axis image; 
set(gca,'xtick',[],'ytick',[])
title({'MAE(our depth)',['=' num2str(err_map) ' cm']})

fprintf('* Finished plotting reconstruction results!\n') 

