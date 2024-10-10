% This script explains how to use the RMM software package.
% This code was developed and tested by Matlab R2021a.
% 
% Code developed by Sungsam Kang (Korea University Superdepth imaging group)
%                                                              2024. 6. 25.
%% set function path
% Run this section to add CLASS functions and sample data folder in Matlab 
% search paths. 
path = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(path));
%% load test data
% Load numerical simulation data included in the RMM software package.
% sample: 'cameraman.tif' image with optinal input/output aberrations
load('numerical_data_RMM_simulation.mat');


% E_in: Incident electric field 
% E_out: Reflected electric field by the sample
% E_in_ab: Incident electric field under input aberration
% E_out_ab: Reflected electric field under input aberration
% L: ROI size (um)
% dx: pixel resolution (um)
% Nx, Ny: number of pixels 
% NA: numerical aperture
% lambda: wavelength
% nb: number of images measured in the experiment while scanning angle
% ground_truth_data: data struct for ground-truth images used in the simulation
%       dx: pixel resolution of ground-truth target image
%       Target: ground-truth target image
%       input_aberration, output_aberration: ground-truth aberration maps

%% construct reflection matrix without aberration
% convert individual 2D images into column vectors
E_in = reshape(E_in, Nx*Ny, nb);
E_out = reshape(E_out, Nx*Ny, nb);

% construct reflection matrix
RR = E_out * E_in';

% construct reflection matrix in k-space
Rk = RM_fft(RR, Nx, Ny);

%% construct reflection matrix with aberration
% convert individual 2D images into column vectors
E_in_ab = reshape(E_in_ab, Nx*Ny, nb);
E_out_ab = reshape(E_out_ab, Nx*Ny, nb);

% construct reflection matrix
RR_ab = E_out_ab * E_in_ab';

% construct reflection matrix in k-space
Rk_ab = RM_fft(RR_ab, Nx, Ny);

%% generate logical index map for dk-space mapping
indx_dk = RM_get_dk_logical_index(Nx, Ny);

%% generate RMM image 
[RMMimg, RMMimg_deconv] = RM_get_image(Rk, indx_dk, Nx, Ny);
[RMMimg_ab, RMMimg_ab_deconv] = RM_get_image(Rk_ab, indx_dk, Nx, Ny);

% 
figure(1), imagesc(ground_truth_data.Target);
axis image;colormap(gray);axis off;colorbar;
title('ground-truth target image');

%
figure(2), imagesc(abs(RMMimg));
axis image;colormap(gray);axis off;colorbar;
title('RMM image: no aberration case');

%
figure(3), imagesc(abs(RMMimg_deconv));
axis image;colormap(gray);axis off;colorbar;
title('deconvolved RMM image: no aberration case');
%
figure(4), imagesc(abs(RMMimg_ab));
axis image;colormap(gray);axis off;colorbar;
title('RMM image: with aberration case');

%
figure(5), imagesc(abs(RMMimg_ab_deconv));
axis image;colormap(gray);axis off;colorbar;
title('deconvolved RMM image: with aberration case');

%% CLASS application
% Apply CLASS algorithm for correcting input and output aberrations
[CLASS_img, Rk_cor, CLASS_result]=CLASS(Rk_ab, indx_dk, Nx, Ny);

%% Display CLASS result
% display CLASS image
figure(6), imagesc(abs(CLASS_result.Final_image_deconv));
axis image;colormap(gray);axis off;colorbar;
title('CLASS image');

% display input and output aberration maps obtained by CLASS algorithm
figure(7), imagesc(angle(CLASS_result.Ain_tot));
axis image;colormap(jet);axis off;colorbar;
title('input aberration map');

figure(8), imagesc(angle(CLASS_result.Aout_tot));
axis image;colormap(jet);axis off;colorbar;
title('outpout aberration map');
