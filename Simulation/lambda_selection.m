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
normD = norm(D,'fro');
noise_data = reshape(D,[h,w,band]);
constant = [0.1,0.2,0.5,0.8,1,1.2,1.5,2,3,5,8,10];
len = length(constant);
err = zeros(len,1);
%% CTV with varing lambda
for it = 1:len
    clear opts;
    opts.lambda = 3*constant(it)/sqrt(h*w);
    out = ctv_rpca(noise_data,opts);
    A_hat = reshape(out,[h*w,band]);
    err(it) = norm(A_hat-RLmat,'fro')/normD;
end
figure;plot(constant,err,'r-*');
xlabel('constant of lambda')
ylabel('relative error')