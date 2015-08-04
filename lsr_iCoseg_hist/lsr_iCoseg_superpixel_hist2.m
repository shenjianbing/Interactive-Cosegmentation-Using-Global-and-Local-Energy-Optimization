function [YS,testParam] = lsr_iCoseg_superpixel_hist2(histSP,superpixel_labelInd,edges_s,colors_s,param)
%LSR_ICOSEG_SUPERPIXEL Summary of this function goes here
%  histSP{i};%d x n
%  colors{i};%n x d
%  superpixel_labelInd{i}.fg;%n x 1
%  superpixel_labelInd{i}.bg;%n x 1
%  edges_s{i};% the index of adjacence superpixel
%  YS: the segmentation results for superpixels
   options.type_interaction = 1; 
   options.gamma_global     = param.lambda0;%100000.0;% lambda0 in paper
   options.lambda_local     = 0.0001; % lambda in paper
   options.is_labeled = 1;
   options.n_unlabeled_sp = 0;
lambda2 = param.lambda2;
lambda1 = param.lambda1;
lambda0 = options.gamma_global;
niter = param.niter;
n_img = length(histSP);
n_scri = length(superpixel_labelInd);


%% scribbles images are the supervised images
n_sp = zeros(n_img,1);%the number of superpixel in each image
superpixel_labelInd_all.fg = [];
superpixel_labelInd_all.bg = [];
M = cell(n_img,1);
empty_labelInd.hg=[1];
empty_labelInd.bg=[1];
S=0;

%get all scribbled superpixels
scribbled_sp_colors_all_fg = [];
scribbled_sp_colors_all_bg = [];
scribbled_sp_hist_all_fg = [];
scribbled_sp_hist_all_bg = [];
for i = 1:n_scri
    scribbled_sp_colors_all_fg = [scribbled_sp_colors_all_fg;colors_s{i}(superpixel_labelInd{i}.fg,:)];
    scribbled_sp_colors_all_bg = [scribbled_sp_colors_all_bg;colors_s{i}(superpixel_labelInd{i}.bg,:)];
    scribbled_sp_hist_all_fg = [scribbled_sp_hist_all_fg,histSP{i}(:,superpixel_labelInd{i}.fg)];
    scribbled_sp_hist_all_bg = [scribbled_sp_hist_all_bg,histSP{i}(:,superpixel_labelInd{i}.bg)];
end
scribbled_sp_hist_all = [scribbled_sp_hist_all_fg,scribbled_sp_hist_all_bg];%d x n
scribbled_sp_colors_all = [scribbled_sp_colors_all_fg;scribbled_sp_colors_all_bg]';%d x n
scribbled_sp_all_ind.fg = [1:size(scribbled_sp_colors_all_fg,1)]';
scribbled_sp_all_ind.bg = [size(scribbled_sp_colors_all_fg,1)+1:size(scribbled_sp_colors_all,2)]';

% labelParam.fea_style = 'histogram';

%makeGMM colors
if ~exist('EM.m','file')
    addpath('GMM-GMR-v2.0');
end
A = scribbled_sp_colors_all;
bgPix = scribbled_sp_all_ind.bg ;
fgPix = scribbled_sp_all_ind.fg ;
ncenters = param.ncenter;
GMMparam.ncenters = ncenters;
X = A(:,bgPix);
[Unaries, Mu, Sigma, ncentersB] = EM_init_kmeans(X, ncenters);
[GMMparam.unarypotB, GMMparam.muB, GMMparam.sigmaB] = EM(X, Unaries, Mu, Sigma);
GMMparam.ncentersB = ncentersB;

X = A(:,fgPix);
[Unaries, Mu, Sigma, ncentersF] = EM_init_kmeans(X, ncenters);
[GMMparam.unarypotF, GMMparam.muF, GMMparam.sigmaF] = EM(X, Unaries, Mu, Sigma);
GMMparam.ncentersF = ncentersF;
labelParam = GMMparam;
labelParam.fea_style = 'colors';

% initialization
S=0;
% scale = ones(n_img,1);
start = zeros(n_img,1);
start(1)=1;
d = size(histSP{1},1);
h_bar = zeros(d,1);

E0 = 0;
e = cell(n_img,1);
YS = cell(n_img,1);
M_HtH = cell(n_img,1);
test_label_global = cell(n_img,1);
testGMMprob = cell(n_img,1);
label_index = cell(n_img,1);
M_intra = cell(n_img,1);
Y0 = cell(n_img,1);
stop = zeros(n_img,1);
% main routine
for i = 1:n_img
    tempC_i = colors_s{i};
    tempC_i = tempC_i';%d x n
    temp_i = histSP{i};%d x n
    n_sp(i)=size(temp_i,2);
    e{i} = ones(n_sp(i),1);
    YS{i} = e{i}*0.6;

%     M_inter = (n_img-1)*temp_i'*temp_i;  
    M_HtH{i} = sparse(temp_i'*temp_i);
    stop(i) = start(i) + n_sp(i) - 1;
    if i == n_img
        n = stop(i);
    else
        start(i+1) = stop(i) + 1;
    end  
    
    
    
    if i<=n_scri
%         [ label_global_ind ] = get_labels_from_scribbles( temp_i,scribbled_sp_hist_all, scribbled_sp_all_ind,labelParam );
        [ label_global_ind GMMprob] = get_labels_from_scribbles( tempC_i,scribbled_sp_colors_all, scribbled_sp_all_ind,labelParam );
        %test
        test_label_global{i} = label_global_ind;%just for test
        testGMMprob{i} = GMMprob;%just for test
        
        label_local_ind = superpixel_labelInd{i};
        % keep global labels and local labels consistent
        logical_fg = zeros(n_sp(i),1);              logical_bg = zeros(n_sp(i),1);
        logical_fg(label_global_ind.fg) =1;         logical_bg(label_local_ind.bg) =1;
        temp_log = logical_fg+logical_bg;
        logical_fg(temp_log==2) = 0; 
        label_global_ind.fg = find(logical_fg);
        
        logical_fg = zeros(n_sp(i),1);              logical_bg = zeros(n_sp(i),1);
        logical_fg(label_local_ind.fg) = 1;         logical_bg(label_global_ind.bg) =1;
        temp_log = logical_fg + logical_bg;
        logical_bg(temp_log==2) = 0; 
        label_global_ind.bg = find(logical_bg);  
        
        %combine global labels and local labels 
        label_ind.fg = unique([label_local_ind.fg;label_global_ind.fg]);
        label_ind.bg = unique([label_local_ind.bg;label_global_ind.bg]);        

        [ M_intra{i} ,Y0{i}] = construct_intra_matrix( colors_s{i}', label_ind, edges_s{i} );
        
        Y1 = zeros(size(Y0{i})); Y1(label_global_ind.fg) = lambda0-lambda1;
        M_intra{i} = M_intra{i} - diag(Y1);  Y0{i} = Y0{i} - Y1;

         %GMM constraint
        label_index{i} = label_ind;
        
        YS{i} = M_intra{i}\Y0{i};
        YS{i} = double(YS{i}>=0.5);
        h_bar = h_bar+temp_i*YS{i};
        %%%compute energy
        Mi = M_intra{i}-diag(Y0{i})+lambda2*M_HtH{i};
        ys_new = YS{i};
        Vi = lambda2*temp_i'*temp_i*ys_new;
        E0 = E0+ys_new'*Mi*ys_new-2*ys_new'*Vi;

    else

        [ label_global_ind  GMMprob] = get_labels_from_scribbles( tempC_i,scribbled_sp_colors_all, scribbled_sp_all_ind,labelParam );
        %test
        test_label_global{i} = label_global_ind;%just for test
        testGMMprob{i} = GMMprob;%just for test
    
        [ M_intra{i} ,Y0{i}] = construct_intra_matrix( colors_s{i}', label_global_ind, edges_s{i} );
        Y1 = zeros(size(Y0{i})); Y1(label_global_ind.fg) = lambda0-lambda1;
        M_intra{i} = M_intra{i} - diag(Y1);  Y0{i} = Y0{i} - Y1;
         %GMM constraint
        label_index{i} = label_global_ind;    
        
        YS{i} = (M_intra{i} )\(Y0{i});           
        YS{i} = double(YS{i}>=0.5);     
        h_bar = h_bar+temp_i*YS{i};
        %%%compute energy
        Mi = M_intra{i}-diag(Y0{i})+lambda2*M_HtH{i};
        ys_new = YS{i};
        Vi = lambda2*temp_i'*temp_i*ys_new;
        E0 = E0+ys_new'*Mi*ys_new-2*ys_new'*Vi;
 
    end
  
end

%% iterative approach
h_bar = h_bar/n_img;
count = 1;
options_c=optimset('Algorithm','interior-point','display','off');

while count<=niter 
    
    E1 = 0;
    for i=1:n_img
        temp_i = histSP{i};%d x n
        Mi = M_intra{i}-diag(Y0{i})+lambda2*M_HtH{i};
        Vi = lambda2*temp_i'*h_bar;
        ys_old = YS{i};        

        label_ind = label_index{i};
        lb = zeros(n_sp(i),1);lb(label_ind.fg)=1;
        ub = ones(n_sp(i),1);ub(label_ind.bg)=0;
        [ys_new,Ei] = fmincon(@(ye_old)myfun2(ys_old,Mi,Vi),ys_old,[],[],[],[],lb,ub,[],options_c);

        Yp{i} = ys_new;
        ys_new = double(ys_new>=0.5);
        YS{i} =ys_new;
        h_bar = h_bar+temp_i*(ys_new-ys_old)/n_img;
        Vi = lambda2*temp_i'*h_bar;
        E1 = E1+ys_new'*Mi*ys_new-2*ys_new'*Vi;
    end
    EE{count} = E1;

    fprintf('enger=%d\n',E1);
    count  = count+1;
    if abs(E1-E0)<0.001
        break
    else
        E0 = E1;
    end
end
fprintf('niter=%d\n',count-1);

testParam.test_label_global = test_label_global;
testParam.testGMMprob = testGMMprob;
testParam.EE = EE;
testParam.Yp = Yp;
for i = 1:n_img
    YS{i} = double(YS{i}>=0.5);
end

end





