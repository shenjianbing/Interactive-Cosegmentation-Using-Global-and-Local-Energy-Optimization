%    version 1.0, 
%    rewritten and updated by Shiming Xiang,  Oct. 6, 2010  

function MR   = regularized_matrix(M, ind_labeled, lambda)

% M:                the main matrix, is a sparse matrix
% ind_labeled:      the indices of the data points which are labeled by the
%                   user, including the labeled possitive and negative
%                   samples. It is a row vector.
% lambda:           the regularization parameter. Usually, it is a large number

% return:
% MR:               the regularized matrix

% N = size(M,1);
% MR = M +  lambda * speye(N);
% 
% return;
% 





N = size(M,1);
num_labeled = length(ind_labeled);

MR = M;
for i = 1:num_labeled
    ind_temp = ind_labeled(i);
    MR(ind_temp, ind_temp) = MR(ind_temp, ind_temp) + lambda;   
end





