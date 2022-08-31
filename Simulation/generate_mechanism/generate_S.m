function [outer_matrix] =generate_S(row,col,k_num)

% Oct 2021
% written by Jiangjun Peng

outer_matrix = zeros(row,col);
rand_index   = randperm(row*col);
choose_index = rand_index(1:k_num);
for i =1:length(choose_index)
    if rand(1,1)>=0.5
        outer_matrix(choose_index(i))=1;
    else
        outer_matrix(choose_index(i))=-1;
    end
end
    
