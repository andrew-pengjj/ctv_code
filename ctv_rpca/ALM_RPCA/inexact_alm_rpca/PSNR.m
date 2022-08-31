%% PSNR 
function [mse,psnr]=PSNR(O,C)
%    input 
%          O:ԭͼ
%          C�����Ƚϵ�ͼƬ
%    output
%          mse:������
%          psnr:��ֵ�����
mse=sum(sum(abs(O-C).*abs(O-C)))/(size(O,1)*size(O,2));
psnr=20*log10(1/sqrt(mse));