%% ���ͼƬ�е����岿��
function P_Image=remove_edge(D)
%  input: 
%         D��ԭʼͼ�񣬴��б߿�
%  output
%         P_Image:�����������߿��ͼƬ

ind1=mean(D,1);
column=find(ind1~=255);
tmp=D(:,column);
ind2=mean(D,2);
row=find(ind2~=255);
P_Image=tmp(row,:);