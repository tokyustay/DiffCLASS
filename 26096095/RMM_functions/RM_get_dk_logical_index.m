function indx_dk = RM_get_dk_logical_index(Nx, Ny)
% Get indices of delta-k space mapping for logical indexing. 
% input:
%       Nx: Number of pixels of ROI along x-axis. Should be even number.
%       Ny: Number of pixels of ROI along y-axis. Should be even number.
%       
% output:
%       indx_dk: indices for logical indexing of delta-k mapping. [int32]
%
% Example:
%       For spatial frequency domain reflection matrix M_kin_kout whose
%       elements ~ [Nx*Ny, Nx*Ny] with sampling ROI [Nx, Ny]
%       mapping delta-k space:
%           M_delta_k_space = zeros(4*NxNy, Nx*Ny, 'single');
%           M_delta_k_space(indx_dk) = M_kin_kout;
%       inverse mapping (delta-k space matrix -> spatial frequency domain
%       matrix)
%           M_kin_kout = M_delta_k_space(indx_dk);
%           M_kin_kout = reshape(M_kin_kout, Nx*Ny, Nx*Ny);
% 
% Code developed by Sungsam Kang (Korea University Superdepth imaging group)
%                                                              2024. 06. 25
%       

if (~rem(Nx, 2)) && (~rem(Ny, 2))
    Nky=Ny/2;
    Nkx=Nx/2;
else
    error('Nx, and Ny should be even integer!!');
end


[kx, ky]=meshgrid(-Nkx:Nkx-1, -Nky:Nky-1);
kx=kx(:);
ky=ky(:);
n1=(2*Nkx)*(2*Nky); % number of elements in reflection matrix
dk_matrix=zeros(n1*4, n1);
kmap=padarray(ones(Nky*2, Nkx*2), [Nky, Nkx]);
for ii=1:n1
    tmp_in=circshift(kmap, [-ky(ii), -kx(ii)]);
    dk_matrix(:, ii)=tmp_in(:);
end
indx_dk = uint32(find(dk_matrix==1));
end