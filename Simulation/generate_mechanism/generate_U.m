function [mask,coef_matrix]=generate_U(n1,n2,rk,thresK)

% Oct 2021
% written by Jiangjun Peng

mask = zeros(n1,n2);
class = rk;
if rk > thresK
    rk = thresK;
end
coef_tensor = zeros(n1,n2,rk);
support = normrnd(0,1,class,rk);
init_center = randperm(n1*n2,class);
center_axis = zeros(class,2);
for j = 1:class
    row_id = mod(init_center(j),n1);
    col_id = (init_center(j)-row_id)/n1+1;
    center_axis(j,1) = row_id;
    center_axis(j,2) = col_id;
end
for i = 1:n1
    for j = 1:n2
        a = repmat([i,j],class,1);
        dist = sum((a-center_axis).*(a-center_axis),2);
        [c,c_id]=min(dist);
        mask(i,j)=c_id;
        coef_tensor(i,j,:)=support(c_id,:);
    end
end
coef_matrix = reshape(coef_tensor,[n1*n2,rk]);


            