clear all;
clc;

N = 512^2; %signal length
scal = 1/4;
n = N*(scal^2);
%s = 2000;
%load image
image = 'lovett';
if strcmp(image, 'lovett')
    im = imread('test_images/Lovett_Hall.jpg');
% im = imread('test_images/Catt_Hall3.jpg');
%     im = imresize(im,[512 512]);
    im = im2double(im);
    im = im(:,(400-250):(400+250),:);
    %im = im(:,(500-250):(500+250),:);
    nn = 512;
    im_sq = imresize(im,[nn nn],'bicubic');
    imgray = rgb2gray(im_sq);
else
    error('Unknown image')
end           
imgray = imresize(imgray,scal);
nn = scal*nn;
lev = floor(log(nn)/log(2));

%important: 'per' ensures you get a non-redundant transform.
wavmode = 'per';
wavname = 'haar';
dwtmode(wavmode,'nodisp')
[c,c_ind] = wavedec2(imgray,lev,wavname);


pr = struct;
%Fixed parameters
%pr.n = 1000; %length of the input signal
pr.b = 1; %number of blocks if signal is block-sparse; otherwise keep 1
pr.tol1 = 1e-5; %error tolerance for measurements
pr.tol2 = 1e-7;
pr.max_iter = 15;
pr.R = 4.25; %period of the modulo function
pr.rho = 3;%spread of the true measurements, y =A*z
pr.del = 1; %truncation factor for supp estimation
pr.spgl_opts = spgSetParms('verbosity',0);
%Tuned parameters
%pr.mspan1 = [100:100:2000];
%pr.mspan2 = [600:100:1000];
%pr.mspan=[pr.mspan1,pr.mspan2];
%pr.mspan=100:100:1000;
pr.mspan = 6000;
pr.num_trials = 1;
pr.s_span = 800:800:800; % sparsity
pr.amp = 1; %amplification factor 
pr.del_p = 0.005; % ps = del*m (sparsity pertaining to error in p)
pr.method = 'basis-pursuit';
pr.init_method = 'true-rcm';
pr.svd_opt = 'svd';
pr.plot_method = 'mean-error';

s = pr.s_span(1);
m=pr.mspan(1);
abs_c = abs(c');
[sorted_c, indx_c] = sort(abs_c,'descend');

z = zeros(size(c',1),1);
z(indx_c(1:s)) = c(indx_c(1:s)); 
supp_c = indx_c(1:s);
norm_z = norm(z);
z = z/norm(z);
[y_clip, y_p, A] = modulo_measure_signal(m,z,pr.R);
%[y_clip, y, A] = clip_measure_signal(m,z,pr.R);
y_clip_nz = nonzeros(y_clip);
nzind=find(y_clip);
mnz=nnz(y_clip);
switch pr.method
    case 'cosamp'
        x = cosamp((y_clip)/sqrt(m), A/sqrt(m),s,100); %Its = 100
        disp('Error with CoSaMP')
        norm(x-z)/norm(z)
        
    case 'cosamp-nz'
        x = cosamp((y_clip_nz)/sqrt(mnz), A(nzind,:)/sqrt(mnz),s,100); %Its = 100
        disp('Error with CoSaMP-NZ')
        norm(x-z)/norm(z)

    case 'basis-pursuit'
        x = spg_bp(A/sqrt(m),(y_clip)/sqrt(m),pr.spgl_opts);
         disp('Error with BP')
        norm(x-z)/norm(z)
        
    case 'basis-pursuit-nz'
        x = spg_bp(A(nzind,:)/sqrt(mnz),(y_clip_nz)/sqrt(mnz),pr.spgl_opts);
         disp('Error with BP-NZ')
        norm(x-z)/norm(z)
        
        %x = l1eq_pd(x,A/sqrt(m), [], (y_mod-pr.R*p)/sqrt(m)); % l1-magic implementation -- slow
 end

z = norm_z*z;
x = norm_z*x;
im1 = waverec2(c,c_ind,wavname);
im2 = waverec2(z, c_ind, wavname);
im3 = waverec2(x,c_ind, wavname);
figure, imshow([im1 im2 im3]), axis image;
title('Actual v/s Sparse v/s Reconstructed');
psnr_err = psnr(im3,im2);
save(['./lovett/rconst_clip_',pr.init_method,'_amp_',num2str(pr.amp),'_r_',num2str(pr.R),'_s_',...
num2str(pr.s_span(1)),'_',num2str(pr.s_span(end)),'_m_',num2str(pr.mspan(1)),...
'_',num2str(pr.mspan(end)),'_',pr.method,'_num_trials_',num2str(pr.num_trials),num2str(psnr_err)],'im2','im3','psnr_err');



% p_err = y_p - p;
% p_err_idx = find(p_err~=0);
% y_true = A*z;
% disp('value of y corresponding to the errors are:')
% disp(y_true(p_err_idx))
% 
% 
% y_sorted = sort(y_true, 'ComparisonMethod','abs');
% y_sorted(1:length(p_err_idx))



% construct_subplots(reconst_err,pr,['rconst_',pr.init_method,'_amp_',num2str(pr.amp),'_r_',num2str(pr.R),'_s_',...
%     num2str(pr.s_span(1)),'_',num2str(pr.s_span(end)),'_m_',num2str(pr.mspan(1)),...
%     '_',num2str(pr.mspan(end)),'_',pr.method,'_num_trials_',num2str(pr.num_trials)],pr.plot_method,1);

% construct_plot(init_err,pr,['init_','r_',num2str(pr.R),'_s_',...
%     num2str(pr.s_span(1)),'_',num2str(pr.s_span(end)),'_m_',num2str(pr.mspan(1)),...
%     '_',num2str(pr.mspan(end)),'_',pr.method],pr.method);