% Matlab script for generating reflection matrix simulation data.
% This script simulates reflection matrix data under known ground-truth target
%
% This script require the CLASS functions and sample data folder. 
%
% This code was developed and tested by Matlab R2021a.
% 
% Code developed by Sungsam Kang (Korea University Superdepth imaging group)
%                                                             2024. 06. 25
%% [0] set function path
% Run this section to add CLASS functions and sample data folder in Matlab 
% search paths. 
path = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(path));

%% [1] set simulation parameters
% Sampling period and number of sampling points are determined by the
% Nyquist sampling limit for a given maximum bandwidth of the optical
% system (by numerical aperture of the objective lens) and size of the ROI.
% All the units in position domain are set to micro-meter scale.

lambda = 1.0; % wavelength [um]
n0 = 1.33; % refractive index of immersion medium, assume water
NA = 1.2; % numerical aperture of objective lens
Lx = 25; % ROI size along x axis [um]
Ly = 25; % ROI size along y axis [um]

k0 = 2*pi/lambda;  % wavenumber, [um^-1]
kmax = k0*NA; % maximum wavenumber supported by the numerical aperture. [um^-1]
dx = pi/kmax; % optimum sampling period by Nyquist sampling limit. [um]
Nx = round(Lx/dx/2)*2; % number of sampling pixels along x-direction. We set it as even number for CLASS applications
Ny = round(Ly/dx/2)*2; % number of sampling pixels along y-direction. We set it as even number for CLASS applications

dkx = 2*kmax/Nx; % sampling period in spatial frequency domain for x-axis. [um^-1]
dky = 2*kmax/Ny; % sampling period in spatial frequency domain for y-axis. [um^-1]
Nkx = round(kmax/dkx); % number of sampling points for the numerical aperture
Nky = round(kmax/dky);

% Generate aperture mask of objective lens according to the sampling
% condition

x = -Lx/2:dx:(Lx/2-dx);
y = -Ly/2:dx:(Ly/2-dx);
kx = -kmax:dkx:(kmax-dkx);
ky = -kmax:dky:(kmax-dky);

mask = kx.^2 + (ky').^2 < kmax^2;

%% [2] set ground-truth images
% load ground-truth target image (target reflectivity) 
% We set the size of the ground-truth target image to [Ny*2, Nx*2] with
% pixel resolution dx/2.
% It will result the maximum spatial frequency bandwith of the ground-truth 
% target image as 2*k0*NA which coincides with the bandwidth limit of 
% reflection matrix microscopy.
% For the sake of memory, CLASS functions and this simulation code use
% single precision for the floating points.

ground_truth_target_image = single(imread('cameraman.tif'));
ground_truth_target_image = ground_truth_target_image/max(ground_truth_target_image(:));
ground_truth_target_image = imresize(ground_truth_target_image, [Ny*2, Nx*2]);

figure(1), imagesc(ground_truth_target_image);
axis image;colormap(gray);axis off;
title('Ground-truth target image');

% load sample aberration maps stored in the 'RMM_sample_data' folder. One can use
% arbitrary aberration patterns instead of ths sample files.

load('numerical_data_ground_truth_aberrations.mat');
% set size of aberration maps [Nx, Nx]
ground_truth_input_aberration = imresize(ground_truth_ab_in, [Ny, Nx]);
ground_truth_output_aberration = imresize(ground_truth_ab_out, [Ny, Nx]);

ground_truth_input_aberration = ground_truth_input_aberration./abs(ground_truth_input_aberration);
ground_truth_output_aberration = ground_truth_output_aberration./abs(ground_truth_output_aberration );

figure(2), imagesc(angle(ground_truth_input_aberration));
axis image;colormap(jet);axis off;colorbar;
caxis([-pi pi]);
title('Ground-truth input aberration map')

figure(3), imagesc(angle(ground_truth_output_aberration));
axis image;colormap(jet);axis off;colorbar;
caxis([-pi pi]);
title('Ground-truth output aberration map')

%% [3] Simulate reflection matrix with assuming plane wave illumination
% Generate plane wave illumination with dkx, and dky sampling step
input_index = find(mask(:) == 1); % indices for input spatial frequency within NA.
kr = sqrt((kx.^2 + (ky').^2).*mask); % radial spatial frequency (sqrt(kx^2 + ky^2))
[~, ascending_kr_index] = sort(kr(input_index));
input_index = input_index(ascending_kr_index);

Inc = diag(mask(:)); % diagonal matrix of incident angles in spatial frequency domain
Inc = Inc(:, input_index); % illumination angle matrix within NA
Inc = single(reshape(Inc, Ny, Nx, []));
% applay padding with considering the spatial frequency bandwidth
Inc = padarray(Inc, [Nky, Nkx, 0]);
% Generate a stack of incident electric field maps in position domain
% Individual complex maps of Inc corresponds to the incident electric field
% in the focal plane.
E_in = fftshift(fftshift(ifft2(ifftshift(ifftshift(Inc, 1), 2)), 1), 2); 
% Generate output electric field after target plane
E_out = E_in .*ground_truth_target_image;

% Construct reflection matrix from image stack
% Vectorize individual images of 3D image stack and mapping into 2D matrix
E_in = reshape(E_in, Nx*Ny*4, []);
E_out = reshape(E_out, Nx*Ny*4, []);

% Generate reflection matrix under plane wave basis
% In principle, R should be given by E_ref*(E_Inc)^(-1). However, since
% E_Inc is semi-unitary matrix, E_Inc*E_Inc' = I (identity matrix). 
R = E_out*E_in';

% Matrix Fourier transform to convert R into spatial frequency domain
Rk = RM_fft(R, Nx*2, Ny*2);
% Crop matrix in spatial frequency domain with [NxNy, NxNy] dimension
ROI = zeros(Ny*2, Nx*2);
ROI(Ny/2+1:Ny/2*3, Nx/2+1:Nx/2*3) = 1;
ROI = find(ROI==1);

Rk = Rk(ROI, ROI);

% apply aperture mask of objective lens
Rk = Rk.*mask(:).*mask(:)';
% Normalize the reflection matrix by its Frobinius norm
Rk = Rk/norm(Rk(:));
% Generate aberrated reflection matrix in spatial frequency domain
Rk_ab = ground_truth_output_aberration(:).*Rk.*(ground_truth_input_aberration(:).');


%% [3-1] Simulate multiple scattering (white noise) in spatial frequency domain
% generate random numbers 
Nk = single(normrnd(0, 1, Ny*Nx, Ny*Nx)+1i*normrnd(0, 1, Ny*Nx, Ny*Nx));
% apply NA mask
Nk = Nk.*mask(:).*mask(:)';
% normalize Nk by its Frobinius norm
Nk = Nk/norm(Nk(:));

%% [3-2] Construct final reflection matrix with including aberration and multiple scattering
gamma = 0.0; % energy ratio between the reflection signal and multiple scattering noise
Rk_ab_ms = Rk_ab + Nk*sqrt(gamma); 
%% Reconstruct target image from reflection matrix
% generate logical index numbers for dk-space mapping
indx_dk = RM_get_dk_logical_index(Nx, Ny);

[RMMimg, RMMimg_deconv] = RM_get_image(Rk_ab_ms, indx_dk, Nx, Ny);
% [~, image_with_ab_ms] = RM_get_image(Rk_ab_ms, indx_dk, Nx, Ny);

figure(4), imagesc(abs(RMMimg_deconv));
axis image;colormap(gray);axis off;colorbar;
title('Reconstructed target image without aberration');

%% CLASS application
[CLASS_img, Rk_cor, CLASS_result]=CLASS(Rk_ab_ms, indx_dk, Nx, Ny);


figure(6), imagesc(abs(CLASS_result.Final_image_deconv));
axis image;colormap(gray);axis off;
title('Reconstructed target image after CLASS');

