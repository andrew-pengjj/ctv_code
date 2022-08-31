function [Omat,RLmat,RSmat] = generate_M(h,w,band,r_s,rho_s,smooth_flag)

% Oct 2021
% written by Jiangjun Peng

%% parameter 
    % Input 
        % h,w:  the size of each image
        % band: band number
        % r_s: low-rank ratio of L
        % rho_s: sparsity of S
        % sigma: noise level
        % smooth_flag: generate_L mode, 1: local_smooth
        % noise_mode: G,B,U,P, more detail can see generate_N.m
        % lambda: the parameter of noise_mode = P
    % Output
        % RLmat: real low-rank matrix
        % RSmat: real sparse matrix
        % Omat:  observation matrix Omat = RLmat+RSmat+Noise
rk = round(band*r_s);
k_num = round(rho_s*h*w*band);
thresK = round(band*0.5);
RLmat = generate_L(h,w,band,rk,thresK,smooth_flag);
RSmat = generate_S(h*w,band,k_num);
Omat  = RLmat + RSmat;