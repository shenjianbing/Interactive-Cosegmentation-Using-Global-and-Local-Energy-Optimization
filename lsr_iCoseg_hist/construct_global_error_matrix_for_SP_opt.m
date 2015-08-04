%    This code is written by Shiming Xiang
% 
%    Address:    National Laboratory of Pattern Recognition (NLPR), 
%                      Institute of Automation, Chinese Academy of Sciences, Beijing, 100190, China
%    Email:        smxiang@gmail.com

%    This software can be used freely for research purposes.
%    Published reports of research  using this code (or a modified version) should cite the 
%    article that describes the algorithm: 
%
%    Shiming Xiang, Feiping Nie, and Changshui Zhang. 
%    Semi-Supervised Classification via Local Spline Regression. 
%    IEEE Transactions on Pattern Analysis and Machine Intelligence, vol.32, no.11, pp.2039-2053, 2010. 

%    This code has been compiled and tested using matlab 6.5  and matlab     7.0

%    version 1.0, 
%    updated by Shiming Xiang,  Oct. 6, 2010  

% =========================================================================


function M = construct_global_error_matrix_for_SP_opt(X, nb, lambda_local)

% X:        the input matrix, each column is a data point;
% nb:       the neighbor indces, recored in each column. In each column, the
%           first is itself, and then its their neighbors
%

% lambda_local:                    the local regularization parameter, lambda in Eq.(9)

% reutrn:
% M:  a sparse matrix the global reconstructtion error matrix


% [K, N] = size(nb);                                     % K:        the k-nn parameter; N: number of data points 
N = length(nb);

% --------------------------------------------------------------------------------------------------------
% parameter check
% --------------------------------------------------------------------------------------------------------
lambda_here = 0.0001;
if nargin >  3
   lambda_here = lambda_local;
end


% --------------------------------------------------------------------------------------------------------
% assemble the local error matrices
% --------------------------------------------------------------------------------------------------------
% size_window = floor(sqrt(K));
% nzmax = calculate_how_many_nonzeros(height, width,  size_window);
% [pIR, pJC] = calculate_index_for_sparse(height, width,  size_window, nzmax);
% values       = zeros(nzmax, 1);
% 
% size_matrix = height * width;
M= sparse(N,N);
for i = 1 : N
    ind = nb{i};
    V = X(:, ind);                        % get the local patch    
    M_local = construct_local_error_matrix_for_image(V, lambda_here);
    K = length(ind);
    for ii = 1 : K
%          nLow  = pJC(ind(ii));
%     	nHigh   = pJC(ind(ii) + 1) - 1;  
        for jj = 1 : K
% 				nRowIndex = my_binary_search(pIR,  nLow,  nHigh, ind(jj));
%                 values(nRowIndex) = values(nRowIndex) + M_local(ii, jj); 
            M(ind(ii),ind(jj)) = M(ind(ii),ind(jj))+ M_local(ii, jj); 
        end
    end

end

% pIC = zeros(nzmax, 1);
% for i = 1: size_matrix
%     pIC(pJC(i) :  pJC(i + 1) -1) = i; 
% end
% 
% clear pJC;
% 
% M = sparse(pIR, pIC, values, size_matrix, size_matrix, nzmax);


return



