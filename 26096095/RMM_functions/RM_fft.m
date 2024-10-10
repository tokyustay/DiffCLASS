function M_kin_kout=RM_fft(M_xin_xout, Nx, Ny)
% This function calculate the Fourier transform of position domain 
% reflection matrix. 
%
% Input: 
%       M_xin_xout: Position domain square shaped reflection matrix whose 
%           column and raw dimension coincides with Nx^2 x Nx^2.
%       Nx: Number of sampling points of ROI along x-direction
%       Ny: Number of sampling points of ROI along y-direction
%
% Output:
%       M_kin_kout: Spatial frequency domain reflection matrix
%
%
%   See also RM_ifft, RM_get_image
%

if ~ (size(M_xin_xout, 1) == Nx*Ny)
    error('Check matrix dimension!!');
end

M_xin_xout=reshape(M_xin_xout, Ny, Nx, []);
M_kout_xin=fftshift(fftshift(M_xin_xout, 1), 2);
clear M_xin_xout;
M_kout_xin=fft2(M_kout_xin);
M_kout_xin=fftshift(fftshift(M_kout_xin, 1), 2);

M_xin_kout=reshape(M_kout_xin, Ny*Nx, []).';
clear M_kout_xin

M_xin_kout= reshape(M_xin_kout, Ny, Nx, []);
M_xin_kout = fftshift(fftshift(M_xin_kout,1),2);
M_kin_kout = ifft2(M_xin_kout); 

clear M_xin_kout;

M_kin_kout = fftshift(fftshift(M_kin_kout,1),2);
M_kin_kout = reshape(M_kin_kout, Ny*Nx, []);

M_kin_kout=M_kin_kout.';

end