function [imgStack_3D, Hs] = RM_3D_recon(M_kin_kout, z_val, Nx, Ny, dx, n0, k0, NA)
% Reconstruct 3D image from reflection matrix.
%
% input    
%       M_kin_kout: reflection matrix in spatial frequency domain
%       z_val: axial location to be reconstructed. It can be [N, 1] array containing z values. [um]
%       Nx: Number of sampling points of ROI along x-direction
%       Ny: Number of sampling points of ROI along y-direction
%       dx: pixel resolution along x, y direction. [um]
%       n0: refractive index of immersion medium or surrounding material.
%       k0: vaccum wavenumber (2pi/wavelength) [rad/um]
%       NA: numerical aperture
%       
% output
%       imgStack_3D: 3D image stack according to axial location z_val.
%       Hs: Transfer function of numerical propagation with propagation distance 1um
%
% Code developed by Sungsam Kang (Korea University Superdepth imaging group)
%                                                              2024. 6. 25
%
% See also RM_get_dk_logical_index.m, RM_get_image.m


% Use CPU multi-cores for logical indexing (This option requires 'indexing_array_mt.mex') 
multi_thread_logical_indexing = true; 

if (~rem(Nx, 2)) && (~rem(Ny, 2))
    Nky=Ny/2;  % spatial frequency bandwidth for y axis
    Nkx=Nx/2;  % spatial frequency bandwidth for x axis
else
    error('Nx, and Ny should be even integer!!');
end
dkx = 2*pi/(Nx*dx);
dky = 2*pi/(Ny*dx);

% Set aperture mask of objective lens
kx = (-Nkx:(Nkx-1))*dkx;
ky = (-Nky:(Nky-1))*dky;
ky = ky';

APmask = kx.^2 + ky.^2 < k0^2*NA^2;
%%
indx_dk = RM_get_dk_logical_index(Nx, Ny);

imgStack_3D=zeros(Nx*2, Nx*2, length(z_val), 'single');
Hs = [];
for ii=1:length(z_val)
    phs=(exp(1i*z_val(ii)*sqrt(k0^2 * n0^2 - kx.^2 -ky.^2).*APmask));
    kkc = M_kin_kout.*(phs(:).*phs(:).'); % reflection
    Hs(:,:, ii) = phs;
    imgStack_3D(:,:, ii)= RM_get_image(kkc, indx_dk, Nx, Ny);
end

end