function [img, img_deconv]=RM_get_image(M_kin_kout, indx_dk, Nx, Ny)
% Reconstruct image from a spatial frequency domain reflection matrix M_kin_kout.
% M_kin_kout should be a square matrix whose dimension is [Nx^2, Nx^2] with
% assuming imaging .
% logical indexing for delta-k space mapping can be accelerated by using
% 'indexing_array_mt.mex' function which can use multi-thread CPU cores.
% 
% input variable:
%       M_kin_kout: spatial frequency domain reflection matrix whose
%           dimension is [Nx^2, Nx^2] with sampling ROI [Nx, Nx].
%       indx_dk: logical indices for delta-k space mapping. Number of
%           elements should coincides with M_kin_kout (Nx^4). Data type should
%           be int32 for 'indexing_array_mt.mex' function.
%       Nx: Number of pixels of ROI along x-axis. Should be even number.
%       Ny: Number of pixels of ROI along y-axis. Should be even number.
%
%
% Output variables:
%       img: reconstructed image. Its size is doubled by [Ny*2, Nx*2] due 
%           to the increased bandwidth of image reconstruction.
%       img_deconv: Deconvolved reconstructed image with assuming ideal
%           transfer function given by the self-convolution of the aperture
%           function. Default threshold level th = 0.01 is applied to prevent
%           high spatial frequency noise amplification. In case of strong
%           multiple scattering, increase th value for better contrast.
% 
%
% Code developed by Sungsam Kang (Korea University Superdepth imaging group)
%                                                              2024. 06. 25
%
% References
%   CASS microscopy: S. Kang, et. al., Nature Photonics 9 (4), 253-258 (2015)
%                       https://www.nature.com/articles/nphoton.2015.24
%
%   See also RM_get_dk_logical_index, CLASS

sz = size(M_kin_kout, 1); % matrix dimension
th = 0.01; % threshold level for deconvolution


if ~sz==Nx*Ny
    error('input matrix dimension should be [Nx*Ny, Nx*Ny] !!');
end
    
if rem(Nx, 2) || rem(Ny, 2)
    error('Nx, and Ny should be even integer!!');
end

if ~numel(M_kin_kout)==numel(indx_dk)
    error('check dk mapping index size');
end

% check whether using mex file for logical indexing
if exist('indexing_array_mt', 'file') == 3
    dkk = indexing_array_mt(M_kin_kout, indx_dk, sz, sz*4);
else
    dkk = zeros(sz*4, sz, 'single');
    dkk(indx_dk) = M_kin_kout;
end

imgk=(sum(dkk, 2));
imgk = reshape(imgk, Ny*2, Nx*2);
img=ifftshift(ifft2(ifftshift(imgk)));
img = img*4; % factor 4 was applied to match graylevel values with the diagonal elements of M_kin_kout.

% deconvolution
kx = -Nx:(Nx-1);
ky = -Ny:(Ny-1);

% mask = kx.^2 + (ky').^2 < Nx*Ny/4;
mask = kx.^2/(Nx/2)^2 + (ky').^2 /(Ny/2)^2< 1;
H = conv2(mask, mask, 'same');
H = H/max(H(:));
H(H<th) = 1;

img_deconv = ifftshift(ifft2(ifftshift(imgk./H)));

end