clear all;clc;
% addpath(genpath('../CTV_code/')) % linux/MacOS platform
addpath(genpath('..\CTV_code\')) % windows platform
%% load data
hsi_name = 'pure_DCmall_small';
load([hsi_name,'.mat'])
clean_data       = Ori_H;
clean_data       = Normalize(clean_data);
[M,N,p]        = size(clean_data);
tic
gaussian_level = 0.2;
sparse_level   = 0.0;
noise_data       = GetNoise(clean_data,gaussian_level,sparse_level);
D = reshape(noise_data,[M*N,p]);
mpsnr = zeros(3,1);
mssim = zeros(3,1);
ergas = zeros(3,1);
[mpsnr(1),mssim(1),ergas(1)]=msqia(clean_data, noise_data);
%% CTV-RPCA
it =2;
fprintf('======== CTV-RPCA  ========\n')
opts.rho = 1.5;
ctv_out = ctv_rpca(noise_data,opts);
[mpsnr(it),mssim(it),ergas(it)]=msqia(clean_data, ctv_out);
%% RPCA
it =3;
D       = zeros(M*N,p) ;
for i=1:p
    bandp = noise_data(:,:,i);
    D(:,i)= bandp(:);
end
fprintf('========   RPCA  ========\n')
A_hat = rpca_m(D);
rpca_out = reshape(A_hat,[M,N,p]);
[mpsnr(it),mssim(it),ergas(it)]=msqia(clean_data, rpca_out);

showband = 103;
figure;
Y = clean_data(:,:,showband);
subplot(2,2,1);imshow(Y,[]);title('original band')
Y = noise_data(:,:,showband);
subplot(2,2,2);imshow(Y,[]);title(['noise, psnr:',num2str(mpsnr(1))])
Y = ctv_out(:,:,showband);
subplot(2,2,3);imshow(Y,[]);title(['ctv-rpca, psnr:',num2str(mpsnr(2))])
Y = rpca_out(:,:,showband);
subplot(2,2,4);imshow(Y,[]);title(['rpca, psnr:',num2str(mpsnr(3))])