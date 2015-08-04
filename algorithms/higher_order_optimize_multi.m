function [Posteriors n_X n_Y] = higher_order_optimize_multi(linesX, lambda, mu, epsilon, full_connect, hs ,hr, name, img_path, out_path, bbb)

beta = 60;
multi = size(hs,2);

%% information of pixels
[img lab_img points_p edges_p colors_p lab_colors_p] = getPropertiesForPixels([img_path name '.bmp']);

%% over-segmentation
M  = 30;
data_path = [out_path 'regions/']; mkdir(data_path);
for i=1:multi
    [segs{i} labels{i} seg{i} colors_s{i} lab_colors_s{i} edges_s{i}] = msseg(double(img),reshape(lab_img, size(img,1)*size(img,2), size(img,3)),hs{i},hr{i},M,full_connect);
    [imgMasks,segOutline,imgMarkup]=segoutput(im2double(img),double(labels{i})); 
    imwrite(imgMarkup, [out_path 'regions/' name '_' int2str(hs{i}) '_' int2str(hr{i}) '_' int2str(M) '.bmp']);  
    clear imgMasks segOutline imgMarkup;
end;

st=clock;
%% make affinity matrix
[W weights_x n_X n_Y] = getW_multi(lab_colors_p,lab_colors_s,edges_p,edges_s,seg,beta);
[PX PXY PYX PY] = getP_all_multi(W,n_X,n_Y);

%% u vectors
nlabels = size(linesX,2);
linesY = getScribbelsSeg_multi(linesX,labels,n_Y);
u_vec = [linesX ; linesY];  sum_initial_likeli = sum(u_vec,2);
idx1(1:n_X) = lambda/(1+mu+lambda); idx1(find(sum_initial_likeli(1:n_X)==0)) = 0;
idx2(1:sum(n_Y)) = lambda/(1+epsilon+lambda); idx2(find(sum_initial_likeli(n_X+1:end)==0)) = 0;
Omega=diag(sparse([idx1 idx2])); clear idx1 idx2 sum_initial_likeli;

%% joint matries
mu_ba = 1/(1+mu);   epsilon_ba = 1/(1+epsilon);
nP1 = horzcat(mu_ba*PX,(1-mu_ba)*PXY);
nP2 = horzcat((1-epsilon_ba)*PYX,epsilon_ba*PY);
Pi = vertcat(nP1,nP2); clear nP1 nP2;
idx(1,1:size(Pi,1))=1;  E=diag(sparse(idx));  clear idx; data = Omega*u_vec;

st1=clock;
B = E -(E-Omega)*Pi;
Likelihoods = B\data;
fprintf(' took %.2f second\n',etime(clock,st1));

iD = diag(sparse(1./sum(Likelihoods')));
Posteriors = iD*Likelihoods;        clear iD;
fprintf(' took %.2f second\n',etime(clock,st));
