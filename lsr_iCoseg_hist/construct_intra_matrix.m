function [ M_intra ,Y0] = construct_intra_matrix( color_SP, labelInd, edges_s,options )
%CONSTRUCT_INTRA_MATRIX Summary of this function goes here
%   Detailed explanation goes here
if nargin <  4
   options.type_interaction = 1; 
   options.gamma_global     = 10000.0; % lambda0 in paper
   options.lambda_local       = 0.0001; % lambda in paper
   options.is_labeled = 1;
end



% check again as some fields can  also be missed even in the case that "options" is transferred
type_interaction = 1;
if isfield(options, 'type_interaction')
    type_interaction = options.type_interaction;
end

gamma_global = 10000.0;
if isfield(options, 'gamma_global')
    gamma_global = options.gamma_global;
end

lambda_local = 0.0001;
if isfield(options, 'lambda_local')
    lambda_local = options.lambda_local;
end

% --------------------------------------------------------------------------------------------------------
% convert image to matrix
% --------------------------------------------------------------------------------------------------------
X = color_SP;                    %In X, each column is a data point(SP histogram)
[RR, CC] = size(X);
X = X + 0.00000001 * randn(RR, CC);      % this sentence is actually unnecessary although I maintain it here

% --------------------------------------------------------------------------------------------------------
% extract the user labeled superpixels
% -------------------------------------
N = size(color_SP,2);
if  options.is_labeled == 1
    label_foreground = labelInd.fg;
    label_background = labelInd.bg;
    index_labeled = [label_foreground; label_background];
    Y0 = zeros(N,1);

    if type_interaction == 1                                         % interactive image segmentation
        Y0(label_foreground) = 1;
        Y0(label_background) = 0;
    else
        Y0(label_foreground) = 1;                               % interactive image matting
        Y0(label_background) = 0;
    end

    clear label_foreground;
    clear label_background;
else
    Y0 = zeros(N,1);
end
% --------------------------------------------------------------------------------------------------------
% construct the k nearest neighborhood
% --------------------------------------------------------------------------------------------------------
k = 8;
nb = zeros(k,N); % including itself
nb(1,:) = [1:N];
connect_nb = cell(1,N);%the neighbours for each superpixels  每个point邻接的节点的坐标    
temp1 = edges_s(:,1);
temp2 = edges_s(:,2);   
for i = 1:N
    temp =  temp1(temp2==i);
    temp =  [temp;temp2(temp1==i)];
    temp = unique(temp);
    dist = sqdist(X(:,temp),X(:,i));
    [~,index] = sort(dist);
    temp =[i; temp(index)];% the nearest sort according the distance, including itself 根据距离，得到最近的排序 包含自己
    if length(temp)>k
        connect_nb{i} = temp(1:k+1);
    else
        connect_nb{i} = temp;
    end
end


% --------------------------------------------------------------------------------------------------------
% construct the global error matrix 
% --------------------------------------------------------------------------------------------------------
 
% M = construct_global_error_matrix_for_image(X, nb, lambda_local);     % very slow!!!!!!!, commented by Xiang, Oct. 6, 2010


M = construct_global_error_matrix_for_SP_opt(X, connect_nb, lambda_local);  % maybe possible error? when constructing indices for sparse matrix.
                                                                                                                                      % in this case, it is nice to use: 
                                                                                                                                      % M = construct_global_error_matrix_for_image(X, nb, lambda_local); 


clear X;
clear nb;

% --------------------------------------------------------------------------------------------------------
% solve the linear system
% --------------------------------------------------------------------------------------------------------
if  options.is_labeled == 1
    Y0 = gamma_global .* Y0;
    M = regularized_matrix(M, index_labeled, gamma_global);
    M_intra = sparse(M);
else
    M_intra = sparse(M);
end


return;
end

