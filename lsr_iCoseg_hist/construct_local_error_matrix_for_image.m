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


function M = construct_local_error_matrix_for_image(X, lambda)
%X:             all the points in a nieghborhood, here thay are expressed in low-dimensional space
%                 In X, each column is a data point, 
%lambda:     the regularization parameter in each neighborhoods


% return:
%    M:    return the inverse matrix


[d, num] = size(X);

% --------------------------------------------------------------------------------------------------------
% calculate the squared distance
% --------------------------------------------------------------------------------------------------------
X2 = sum(X.^2, 1);
dist2 = repmat(X2, num, 1) + repmat(X2', 1, num) - 2 * X'* X;   % squared distance

% --------------------------------------------------------------------------------------------------------
% constrcu the coefficient matrix 
% --------------------------------------------------------------------------------------------------------
distance = sqrt(dist2 +  0.0000000001)  +   lambda * eye(num);     % to avoid imaginary number by adding a very small number

%construct the coefficient matrix
P = ones(num, 1);
P = [P, X'];                 % add the column

% % the following four sentences are written according to the algorithm
C = zeros(d + 1, d + 1);

A = [distance, P; P', C];   % We should point out that  only the linear polynomials are considered here. 
                                        % This may degrade the accuracy of the spline.
                                        % But this treatemnt will facilate the computation when d is greater than 3 (d > 3)

% --------------------------------------------------------------------------------------------------------
% evaluate the local error matrix 
% --------------------------------------------------------------------------------------------------------
M = pinv(A);                       
M = M(1 : num, 1 : num);  

return





