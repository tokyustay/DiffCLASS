function Rkc = RM_ab_correction(Rk, ab_in, ab_out)
% Correct input and output aberrations from CLASS result.
%
% input    
%       Rk: reflection matrix in spatial frequency domain
%       ab_in, ab_out: input and output aberration obtained by CLASS
%       algorithm
%       
% output
%       Rkc: aberration corrected reflection matrix
%
% Code developed by Sungsam Kang (Korea University Superdepth imaging group)
%                                                              2024. 6. 25
%
% See also CLASS

Rkc=ab_out(:).*Rk.*ab_in(:).';

end