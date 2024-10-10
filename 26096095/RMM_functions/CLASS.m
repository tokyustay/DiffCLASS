function [CLASS_img, M_kin_kout, CLASS_result]= CLASS(M_kin_kout, indx_dk, Nx, Ny)
% CLASS main code
% If you are running this code for the first time, run 'initiize_CLASS_functions.m'.
%   input: 
%   M_kin_kout: kin - kout matrix with number of input and output pixels (2NAsz)^2
%              where NAsz is the number of pixels corresponding to the radius of
%              the numerical aperture (or input and output mask).
%   indx_dk: indices for logical indexing of delta-k mapping. [int32]
%   Nx: Number of sampling points along x-direction. Should be even number
%       (twice Nkx_max).
%   Ny: Number of sampling points along y-direction. Should be even number
%       (twice Nkx_max).
%
%   output:
%   CLASS_img: Final CLASS image
%   M_kin_kout: aberration corrected reflection matrix
%   CLASS_result: struct of CLASS result.
%           CLASS_result.Ain_tot: input aberration map
%           CLASS_result.Aout_tot: output aberration map
%           CLASS_result.Ain: input aberration maps while iteration
%           CLASS_result.Aout: output aberration maps while iteration
%           CLASS_result.images: CLASS image while iteration
%           CLASS_result.Final_image: Final CLASS image after deconvolution
%           CLASS_result.rms: rms phase of updated input and output phase
%                             maps while iteration.
%           CLASS_result.enh: intensity enhancement while iteration.
%
% Code developed by Sungsam Kang (Korea University Superdepth imaging group)
%                                                              2024. 06. 25
%
% References
%   CASS microscopy: 
%       S. Kang, et. al., Nature Photonics 9 (4), 253-258 (2015)
%               https://www.nature.com/articles/nphoton.2015.24
%   CLASS microscopy:
%       S. Kang, et. al., Nature Communications 8, 2157 (2017)
%               https://www.nature.com/articles/s41467-017-02117-8
%
%
%   See also RM_fft, RM_ifft, RM_get_image, RM_get_dk_logical_index

t_start = tic;
%% Set CLASS internal parameters
% iteration parameters
itN = 6; % iteration number
pow_itN = 5; % power iteration number for phase correlation

% Low spatial frequency blocking for fast convergence
kfilter = 6; % block size in spatial frequency domain. Set 0 if do not want this option.

% Use CPU multi-cores for logical indexing (This option requires 'indexing_array_mt.mex') 
multi_thread_logical_indexing = true; 

% progress display option
display_progress_fig = true;
display_progress_fig_num = 10001;

% Internal parameters
if (~rem(Nx, 2)) && (~rem(Ny, 2))
    Nky=Ny/2;  % spatial frequency bandwidth for y axis
    Nkx=Nx/2;  % spatial frequency bandwidth for x axis
else
    error('Nx, and Ny should be even integer!!');
end

%% Set aperture mask of objective lens
kx = -Nkx:(Nkx-1);
ky = -Nky:(Nky-1);

APmask = kx.^2/Nkx^2 + (ky').^2 /Nky^2< 1;

kx2 = -2*Nkx:(2*Nkx-1);
ky2 = -2*Nky:(2*Nky-1);
HPFmask = kx2.^2 + (ky2').^2 > kfilter^2;
HPFindex = find(HPFmask(:) == 1);

%% Prepare output struct
CLASS_result.Ain_tot = ones(Ny, Nx, 1, 'single');
CLASS_result.Aout_tot = ones(Ny, Nx, 1, 'single');
CLASS_result.Ain = ones(Ny, Nx, itN, 'single');
CLASS_result.Aout = ones(Ny, Nx, itN, 'single');
CLASS_result.images = zeros(Ny*2, Nx*2, itN+1, 'single');
CLASS_result.rms = zeros(itN, 2, 'single');
CLASS_result.enh = zeros(itN+1, 1, 'single');
%% prepare figure window for displaying iteration step
if display_progress_fig
    hfig = figure(display_progress_fig_num);
    set(hfig, 'pos', [50 749 1160 240 ]);
    tt = tiledlayout(1, 5, 'TileSpacing', 'compact', 'Padding', 'none');
end
%% CLASS main iteration
for ii = 1:itN
    t_itn = tic;
    % input correction
    % mapping on delta k space
    if multi_thread_logical_indexing
        dkk = indexing_array_mt(M_kin_kout, indx_dk, Nx*Ny, 4*Nx*Ny);
    else
        dkk = zeros(4*Nx*Ny, Nx*Ny, 'single');
        dkk(indx_dk) = M_kin_kout;
    end
    ref_vec = sum(dkk, 2);
    
    if ii == 1
        initial_norm = norm(ref_vec);
    end

    CLASS_result.enh(ii) = norm(ref_vec)/initial_norm;
    % set output image while iteration
    img_while_iteration = ifftshift(ifft2(ifftshift(reshape(ref_vec, Ny*2, Nx*2))));

    % calculate enhancement of the previous iteration step
  
    % get input phase correlation with power iteration
    if kfilter == 0
        ph = PowIt(dkk, pow_itN, APmask);
    else
        ph = PowIt(dkk, pow_itN, APmask, HPFindex);
    end       
    
    abi = exp(1i*angle(reshape(ph, Ny, Nx)).*APmask);
    CLASS_result.Ain(:,:, ii) = abi;
    CLASS_result.Ain_tot=CLASS_result.Ain_tot.*CLASS_result.Ain(:,:, ii);
    
    clear dkk;
    
    % matrix correction
    M_kin_kout=M_kin_kout.*abi(:).';
    
    
    % output correction
    % mapping delta-k space for the phase conjugation process
    if multi_thread_logical_indexing
        dkk = indexing_array_mt(M_kin_kout.', indx_dk, Nx*Ny, Nx*Ny*4);
    else
        dkk = zeros(Nx*Ny*4, Nx*Ny, 'single');
        dkk(indx_dk) = M_kin_kout.';
    end
    
    % get output phase correlation with power iteration   
    if kfilter == 0
        ph = PowIt(dkk, pow_itN, APmask);
    else
        ph = PowIt(dkk, pow_itN, APmask, HPFindex);
    end 
    abo = exp(1i*angle(reshape(ph, Ny, Nx)).*APmask);
    CLASS_result.Aout(:,:, ii) = abo;
    CLASS_result.Aout_tot=CLASS_result.Aout_tot.*CLASS_result.Aout(:,:, ii);
    clear dkk;

    % matrix correction
    M_kin_kout=abo(:).*M_kin_kout;

    % set output struct
    CLASS_result.images(:,:, ii) = img_while_iteration;
    CLASS_result.rms(ii,1)=rms(angle(abi(APmask == 1)/mean2(abi(APmask == 1))));
    CLASS_result.rms(ii,2)=rms(angle(abo(APmask == 1)/mean2(abo(APmask == 1))));

    % display image
    if display_progress_fig
        CLASS_display(tt, CLASS_result, itN, sprintf('CLASS iteration %d...', ii), ii)
    end
    
    % display progress
    if ii == 1
        txtbyte = fprintf('CLASS iteration ... %d  [%.2fsec]\n', ii, toc(t_itn));
    elseif ii>1
        for tb = 1:(txtbyte-20)
            fprintf('\b');
        end
%         txtbyte = fprintf('CLASS iteration ... %d  [%.2fsec]', ii, toc(t_itn));
        txtbyte = fprintf('%d  [%.2fsec]\n', ii, toc(t_itn)) + 20;
    end
    
end
%% 


[CLASS_img, CLASS_img_deconv] = RM_get_image(M_kin_kout, indx_dk, Nx, Ny);
CLASS_result.images(:,:, end) = CLASS_img;
CLASS_result.Final_image = CLASS_img;
CLASS_result.Final_image_deconv = CLASS_img_deconv;
image_k = fft2(CLASS_img);
CLASS_result.enh(end) = norm(image_k(:))/initial_norm/4;

t_stop = toc(t_start);
fprintf('CLASS finished in [%.2fsec]\n', t_stop);

if display_progress_fig
    CLASS_display(tt, CLASS_result, itN, sprintf('CLASS finished in [%.2fsec]\n', t_stop), itN);
end

end

function ph = PowIt(T, n, APmask, varargin)
% Calculate first singular vector of input matrix T using power method. 
% input variables 
%   T: input matrix
%   n: number of iteration
if gpuDeviceCount("available")>0
    g = gpuDevice;
    s = whos('T');
    if g.AvailableMemory > s.bytes*0.8
        T = gpuArray(T);
        use_gpu = true;
    else
        use_gpu = false;        
    end
end
    
sz0 = size(T); % initial matrix size
m = size(T, 2); % number of input basis
input_basis_indx = find(APmask == 1);
v0 = ones(numel(input_basis_indx), 1)/sqrt(m); % initial vector
v = v0;

T = T(:, input_basis_indx);

if nargin>3
    HPFindex = varargin{1};
    T = T(HPFindex, :);
end

for ii = 1:n    
    % Be care about the position of brackets. Without these brackets, Matlab
    % automatically calculate the matrix multiplication from left-hand-side. 
    % With brackets, it cacluates matrix-vector product in sequence.
    v = T'*(T*v);  
    v = v./abs(v); % for finding phase only solution
end

ph = ones(m, 1);
ph(input_basis_indx) = v;

if use_gpu
    ph = gather(ph);
end
end

function CLASS_display(tt, CLASS_result, itN, img_title, ii)
title(tt,img_title);
% image before correction
ax1 = nexttile(tt, 1);
imagesc(abs(CLASS_result.images(:,:, 1)));colormap(ax1, gray);axis image;axis off;colorbar;
title(ax1, 'Initial image');

% image after correction
ax2 = nexttile(tt, 2);
imagesc(abs(CLASS_result.images(:,:, ii)));colormap(ax2, gray);axis image;axis off;colorbar;
title(ax2, 'CLASS image');


% Aberration maps found while CLASS iteration
ax3 = nexttile(tt, 3);
imagesc(angle(CLASS_result.Ain_tot));colormap(ax3, jet);axis image;axis off;
title(ax3, 'input aberration map');
ax4 = nexttile(tt, 4);
imagesc(angle(CLASS_result.Aout_tot));colormap(ax4, jet);axis image;axis off;
title(ax4, 'output aberration map');

% CLASS performance curves
ax5 = nexttile(tt, 5);
yyaxis left
plot(1:itN, CLASS_result.rms);
ylabel('RMS phase');
yyaxis right
plot(0:(itN), CLASS_result.enh);
xlim([0 itN]);
ylabel('CLASS enhancement');

drawnow;

end
