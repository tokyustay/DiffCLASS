function M_xin_xout=RM_ifft(M_kin_kout, Nx, Ny)
% This function calculate the Fourier transform of position domain 
% reflection matrix. 
%
% Input: 
%       M_kin_kout: Spatial frequency domain square shaped reflection 
%           matrix whose column and raw dimension coincides with Nx^2 x Nx^2.
%       Nx: Number of sampling points of ROI along x-direction
%       Ny: Number of sampling points of ROI along y-direction
%
% Output:
%       M_xin_xout: Position domain reflection matrix
%
%
%   See also RM_fft, RM_get_image
%

if ~ (size(M_kin_kout, 1) == Nx*Ny)
    error('Check matrix dimension!!');
end

M_kin_kout=reshape(M_kin_kout, Ny, Nx, []);
M_xout_kin=ifftshift(ifftshift(M_kin_kout, 1), 2);
clear M_kin_kout;
M_xout_kin=ifft2(M_xout_kin);
M_xout_kin=fftshift(fftshift(M_xout_kin, 1), 2);

M_kin_xout=reshape(M_xout_kin, Ny*Nx, []).';
clear M_xout_kin

M_kin_xout= reshape(M_kin_xout, Ny, Nx, []);
M_kin_xout = ifftshift(ifftshift(M_kin_xout,1),2);
M_xin_xout = fft2(M_kin_xout); 

clear M_kin_xout;

M_xin_xout = fftshift(fftshift(M_xin_xout,1),2);
M_xin_xout = reshape(M_xin_xout, Ny*Nx, []);

M_xin_xout=M_xin_xout.';

end