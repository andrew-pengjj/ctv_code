function [noise_data,E]=GetNoise(original_data,gaussian_level,sparse_level)

% Oct 2021
% written by Jiangjun Peng

if nargin<2
    gaussian_level = 0;
    sparse_level   = 0;
end
if nargin<3
    sparse_level   = 0;
end
noise_data = original_data;
[M,N,p]    = size(original_data);
if length(gaussian_level)==1
    gaussian_level = gaussian_level*ones(1,p);
end
if length(sparse_level)==1
    sparse_level = sparse_level*ones(1,p);
end

%---------------------------------------------------------------
%                       add noise
%---------------------------------------------------------------
% S&P noise
if mean(sparse_level) ~=0
    for i =1:p
        noise_data(:,:,i)= imnoise(noise_data(:,:,i),'salt & pepper',sparse_level(i));
    end
    fprintf('Sparse noise had been added to data, the level is %.4f\n',mean(sparse_level));
    E = reshape(noise_data - original_data,[M*N,p]);
end 
% Gaussian noise
if mean(gaussian_level) ~=0
    for i =1:p
        noise_data(:,:,i)=noise_data(:,:,i)  + gaussian_level(i)*randn(M,N);
    end
    fprintf('Gaussian noise had been added to data, the level is %.4f\n',mean(gaussian_level));
end
