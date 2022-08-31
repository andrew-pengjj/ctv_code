clear all;clc
addpath(genpath('../../CTV_code/')) % linux/MacOS platform
%addpath(genpath('..\..\CTV_code\')) % windows platform
h = 20;
w = 20;
band = 200;
r_s  = 0.1;
rho_s  = 0.3;
smooth_flag = 1;

[D,RLmat,RSmat]=generate_M(h,w,band,r_s,rho_s,smooth_flag);
diff_v = diff3(RLmat,[h,w,band]);
normD = norm(D,'fro');
noise_data = reshape(D,[h,w,band]);
%% CTV-RPCA
out = ctv_rpca(noise_data);
A_hat = reshape(out,[h*w,band]);
it = 1;
err(it) = norm(A_hat-RLmat,'fro')/normD;
%% RPCA
A_hat = rpca_m(D);
it = 2;
err(it) = norm(A_hat-RLmat,'fro')/normD;
%% CTV with varing lambda
ratio = max(2-rho_s*2,1+r_s*2);
clear opts;
opts.lambda = 3*ratio/sqrt(h*w);
out = ctv_rpca(noise_data,opts);
A_hat = reshape(out,[h*w,band]);
it = 3;
err(it) = norm(A_hat-RLmat,'fro')/normD;
%% RPCA with varing lambda
lambda = ratio/sqrt(h*w);
A_hat = rpca_m(D,lambda);
it = 4;
err(it) = norm(A_hat-RLmat,'fro')/normD;
err_ctv = min(err(1),err(3));
err_rpca = min(err(2),err(4));
fprintf('======== Result ============\n')
fprintf('error of ctv-rpca :  %.6f\n ',err_ctv);
fprintf('error of rpca :  %.6f\n ',err_rpca);

