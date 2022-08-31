function [result]=MAUC(image_data1,image_data2)  
%����AUCֵ,image_data1Ϊԭʼ������ǩ,image_data2Ϊ�������õ�����Ϊ����ĸ���
[~,~,p] = size(image_data1);
auc_result = zeros(1,p);
for i=1:p
    auc_result(i)=AUC(image_data1(:,:,i),image_data2(:,:,i));
end
result = mean(auc_result);