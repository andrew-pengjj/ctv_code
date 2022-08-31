function [data,U,V]=generate_L(h,w,band,rk,thresK,smooth_flag)

% Oct 2021
% written by Jiangjun Peng

if smooth_flag ==1
    [~,U]=generate_U(h,w,rk,thresK);
    flag = 0;
    if rk>thresK
        rk = thresK;
        flag = 1;
    end
    V = normrnd(0,1,band,rk);
    for i=1:rk         
        V(:,i) = smooth(V(:,i));
    end
else
    if rk>thresK
        rk = thresK;
    end
    U = normrnd(0,1,h*w,rk);
    V = normrnd(0,1,band,rk);
end
data = U*V';
data = (data -min(data(:)))/(max(data(:))-min(data(:)));
