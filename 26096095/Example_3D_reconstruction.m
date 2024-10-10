% This script explains how to use the RMM software package.
% This code was developed and tested by Matlab R2021a.
% 
% Code developed by Sungsam Kang (Korea University Superdepth imaging group)
%                                                              2024. 6. 25.
%% set function and data path
% Run this section to add CLASS functions and sample data folder in Matlab 
% search paths. 
path = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(path));
%% load test data
% Load numerical simulation data included in the RMM software package.
% sample: 'cameraman.tif' image with optinal input/output aberrations
load('numerical_data_3D_beads.mat');


% Rk: reflection matrix in spatial frequency domain
% dx: pixel resolution (um)
% Nx, Ny: number of pixels 
% NA: numerical aperture
% lambda: wavelength
% k0: vacuum wave number
% n0: refractive index of immersion medium
% ground_truth_data: data struct for ground-truth images used in the simulation
%       dx: pixel resolution of ground-truth target image
%       dz: axial interval of ground-truth target image
%       z_val: z values of ground-truth target image
%       Target: ground-truth target images [y, x, z]

%% generate logical index map for dk-space mapping
indx_dk = RM_get_dk_logical_index(Nx, Ny);
%% generate RMM image 
RMMimg = RM_get_image(Rk, indx_dk, Nx, Ny);
%
figure(1), imagesc(abs(RMMimg).^2);
axis image;colormap(gray);axis off;colorbar;
title('RMM image');
%% 3D reconstruction
zscan = -7.0:0.5:7.0;
imgStack3D = RM_3D_recon(Rk, zscan, Nx, Ny, dx, n0, k0, NA);

%% display section images
z = 3;

zindx_gt = find(ground_truth_data.z_val >= z, 1);
zindx = find(zscan >= z, 1);

figure(1), imagesc(ground_truth_data.Target(:,:, zindx_gt));
axis image;colormap(hot);axis off;colorbar;
title(sprintf('Ground-truth image at z = %.1fum', z));

figure(2), imagesc(abs(imgStack3D(:,:, zindx )).^2, [0 0.02]);
axis image;colormap(hot);axis off;colorbar;
title(sprintf('Reconstructed image at z = %.1fum', zscan(zindx)));
%% display MIP image of 3D stack
%% ground-truth image
[gt_MIPimg, gt_MIP_z_indx] = max(ground_truth_data.Target, [], 3);

gt_zmap = ground_truth_data.z_val(gt_MIP_z_indx);

H = ((gt_zmap)+8)/16;
S = gt_MIPimg;
V = gt_MIPimg;

MIP_disp_gt = hsv2rgb(cat(3, H, S, V));
figure(1)
imshow(MIP_disp_gt, InitialMagnification=300,Interpolation="bilinear");
%% 3D reconstruction by RMM
[MIPimg, MIP_z_indx] = max(abs(imgStack3D).^2, [], 3);
zmap = zscan(MIP_z_indx);
sat_lv = 0.8;

H = ((zmap)+8)/16;
S = (MIPimg/max(MIPimg(:))/sat_lv);
V = (MIPimg/max(MIPimg(:))/sat_lv);
S(S>1) = 1;
V(V>1) = 1;
MIP_disp = hsv2rgb(cat(3, H, S, V));
figure(2)
imshow(MIP_disp, InitialMagnification=300,Interpolation="bilinear");
