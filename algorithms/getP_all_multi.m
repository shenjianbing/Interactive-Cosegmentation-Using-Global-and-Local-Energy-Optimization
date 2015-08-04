function [P11 P12 P21 P22] = getP_all_multi(W,N1,N2)
n_Y_all = sum(N2);

W11 = W(1:N1,1:N1);
W12 = W(1:N1,N1+1:N1+n_Y_all);
W21 = W(N1+1:N1+n_Y_all,1:N1);
W22 = W(N1+1:N1+n_Y_all,N1+1:N1+n_Y_all);

iD = diag(sparse(1./sum(W11')));    P11 = iD*W11;   clear iD;
iD = diag(sparse(1./sum(W12')));    P12 = iD*W12;   clear iD;
iD = diag(sparse(1./sum(W21')));    P21 = iD*W21;   clear iD;
iD = diag(sparse(1./sum(W22')));    P22 = iD*W22;   clear iD;

% P11 = P11 + P11' - P11'*P11;
% P22 = P22 + P22' - P22'*P22;

% row_sum = sum(W11'); row_sum = row_sum.^0.5; iD = diag(sparse(1./row_sum));    P11 = iD*W11*iD;   clear iD row_sum;
% row_sum = sum(W22'); row_sum = row_sum.^0.5; iD = diag(sparse(1./row_sum));    P22 = iD*W22*iD;   clear iD row_sum;