function [result]=AUC(test_targets,output)  
%计算AUC值,test_targets为原始样本标签,output为分类器得到的判为正类的概率
% 均为行或列向量  
test_targets = test_targets(:);
output       = output(:);
[~,I]=sort(output);  
M=0;N=0;  
for i=1:length(output)  
    if(test_targets(i)==1)  
        M=M+1;%正类样本数
    else  
        N=N+1;  %负类样本数
    end  
end  
sigma=0;  
for i=M+N:-1:1  
    if(test_targets(I(i))==1)  
        sigma=sigma+i;%正类样本rank相加  
    end  
end
result=(sigma-(M+1)*M/2)/(M*N);  