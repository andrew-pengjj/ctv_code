function [result]=MAUC(image_data1,image_data2)  
%计算AUC值,image_data1为原始样本标签,image_data2为分类器得到的判为正类的概率
[~,~,p] = size(image_data1);
auc_result = zeros(1,p);
for i=1:p
    auc_result(i)=AUC(image_data1(:,:,i),image_data2(:,:,i));
end
result = mean(auc_result);