function [result]=AUC(test_targets,output)  
%����AUCֵ,test_targetsΪԭʼ������ǩ,outputΪ�������õ�����Ϊ����ĸ���
% ��Ϊ�л�������  
test_targets = test_targets(:);
output       = output(:);
[~,I]=sort(output);  
M=0;N=0;  
for i=1:length(output)  
    if(test_targets(i)==1)  
        M=M+1;%����������
    else  
        N=N+1;  %����������
    end  
end  
sigma=0;  
for i=M+N:-1:1  
    if(test_targets(I(i))==1)  
        sigma=sigma+i;%��������rank���  
    end  
end
result=(sigma-(M+1)*M/2)/(M*N);  