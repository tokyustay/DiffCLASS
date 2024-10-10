% This script demonstrates an example of CLASS algorithm. 
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
% Load sample data included in the RMM software package.
% sample: USAF target under 1.5mm thick glass plate 
load('Experimental_data_USAF.mat');


% E_in: incident electric field measured by a mirror
% E_out: Reflected electric field by the sample
% L: ROI size (um)
% dx: pixel resolution (um)
% Nx, Ny: number of pixels 
% NA: numerical aperture
% lambda: wavelength
% nb: number of images measured in the experiment while scanning angle
%% construct reflection matrix
% convert individual 2D images into column vectors
E_in = reshape(E_in, Nx*Ny, nb);
E_out = reshape(E_out, Nx*Ny, nb);

% construct reflection matrix
RR = E_out * E_in';

% construct reflection matrix in k-space
Rk = RM_fft(RR, Nx, Ny);

%% generate logical index map for dk-space mapping
indx_dk = RM_get_dk_logical_index(Nx, Ny);

%% generate RMM image 
[RMMimg, RMMimg_deconv] = RM_get_image(Rk, indx_dk, Nx, Ny);
%
figure(1), imagesc(abs(RMMimg));
axis image;colormap(gray);axis off;colorbar;
title('RMM image');
set(gca, 'FontSize', 17);

%% CLASS application
% Apply CLASS algorithm for correcting input and output aberrations
[CLASS_img, Rk_cor, CLASS_result]=CLASS(Rk, indx_dk, Nx, Ny);

%% Display CLASS result
% display CLASS image
figure(2), imagesc(abs(CLASS_result.Final_image));
axis image;colormap(gray);axis off;colorbar;
title('CLASS image');
set(gca, 'FontSize', 17);

% display input and output aberration maps obtained by CLASS algorithm
figure(3), imagesc(angle(CLASS_result.Ain_tot));
axis image;colormap(jet);axis off;colorbar;
title('input aberration map');
set(gca, 'FontSize', 17);
figure(4), imagesc(angle(CLASS_result.Aout_tot));
axis image;colormap(jet);axis off;colorbar;
title('outpout aberration map');
set(gca, 'FontSize', 17);
