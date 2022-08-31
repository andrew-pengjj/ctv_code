function [tensor_data_r,original_data,InputTensor]=GetVideoMask(data_name)

% Oct 2021
% written by Jiangjun Peng
fprintf('======== data_name:%s ========\n',data_name);
load(strcat(data_name,'_mask_r.mat'));
GroundTruth       = tensor_data_r;
original_data     = GroundTruth;
GroundTruth(find(GroundTruth~=0))=1;
load(strcat(data_name,'_r.mat'));
InputTensor       = tensor_data_r;
InputTensor       = Normalize(InputTensor);
[M,N,p]    = size(InputTensor);
ForeTensor = zeros([M,N,p]);
for i = 1:p
    tmp           = GroundTruth(:,:,i);
    nonzero_index = find(tmp~=0);
    tmp_input     = InputTensor(:,:,i);
    tmp_fore      = ForeTensor(:,:,i);
    tmp_fore(nonzero_index)    = tmp_input(nonzero_index);
    ForeTensor(:,:,i) = tmp_fore;
end