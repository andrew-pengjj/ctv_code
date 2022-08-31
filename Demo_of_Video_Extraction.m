clear all;close all;clc
addpath(genpath('../CTV_code/')) % linux/MacOS platform
%addpath(genpath('..\CTV_code\')) % windows platform
data_name="airport";
[original_data,gt_fore]=GetVideoMask(data_name);
[M,N,p]=size(original_data);
InputTensor = original_data;
InputMatrix       = reshape(InputTensor,[M*N,p]);
%% CTV-RPCA
it =1;
fprintf('======== CTV-RPCA  ========\n')
opts.rho    = 1.5;
opts.weight = 10;
[ctv_out,E] = ctv_rpca(InputTensor,opts);
E_ctv = reshape(E,[M,N,p]);
auc(it) = MAUC(gt_fore,abs(E_ctv));

%% RPCA
it =2;
fprintf('======== RPCA  ========\n')
[A_hat,E_hat,iter] = rpca_m(InputMatrix);
rec_tensor = reshape(A_hat,[M,N,p]);
E_rpca = reshape(E_hat,[M,N,p]);
auc(it) = MAUC(gt_fore,abs(E_rpca));

index = 4;
figure
subplot(2,2,1);imshow(InputTensor(:,:,index),[]);title('observed')
subplot(2,2,2);imshow(abs(gt_fore(:,:,index)),[]);title('groundtruth')
subplot(2,2,3);imshow(abs(E_ctv(:,:,index)),[]);title(['CTV-RPCA:',num2str(auc(1))])
subplot(2,2,4);imshow(abs(E_rpca(:,:,index)),[]);title(['RPCA:',num2str(auc(2))])