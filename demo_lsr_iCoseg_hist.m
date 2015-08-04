clear all;
close all;

addpath 'msseg'
addpath 'others'
addpath 'algorithms'
addpath 'lsr_iCoseg_hist'
%% settings
type_of_seed = 0; % 0: scribbles or 1: trimatp
full_connect = 0; % 0: adjacency graph on region layer or 1: full connection
reset(RandStream.getGlobalStream);

% parameters of over segmentation
lambda = 100000; mu = 0.002; epsilon = 0.2;
hs = 10; hr = 7;  % hs{2} = 10; hr{2} = 10;   hs{3} = 10; hr{3} = 15;hs{1} = 10; hr{1} = 7;
M = 30;
nbins = 20;%histogram parameters
% parameters of interaction cosegmentation with local spline regression
param.lambda0 = 1e5;
param.lambda1 = 1e3;
param.lambda2 = 1e-7;%0;%
param.gamma = 1e4;
param.niter =5;
param.ncenter=1;
%% paths

dataset = 'cow\';%'iCoseg\skate2\';%'scaleimage\Horse200\';%
img_path = ['Datasets\images\',dataset];
scribbles_path = ['Datasets\scribbles\',dataset];
out_path = ['./results/',dataset];


%% read image names
imgstyle = 'bmp' ;
img_dir = dir([img_path '*.' imgstyle]);
scribbles_dir = dir([scribbles_path '*.bmp']);
% img_names = {};
% scribbles_names = {};
n_img = length(img_dir);
n_scri = length(scribbles_dir);
img_names_t = cell(n_img,1);
for i =1:n_img
    img_names_t{i} = strtok( img_dir(i).name,'.');  
end
scribbles_names = cell(n_scri,1);
scri_img_idx = zeros(n_scri,1);%index of scribbles image
for i =1:n_scri
    temp_name =  strtok(scribbles_dir(i).name,'.');  
    scribbles_names{i} = temp_name;
    j=1;
    for j = 1:n_img
        if strcmp(temp_name,img_names_t{j})
            scri_img_idx(i) = j;
            break
        end
    end
    
end

%% over segmentation
fprintf('over segmentation\n');
overseg_data = [out_path 'regions/overseg_data_' int2str(hs) '_' int2str(hr) '_' int2str(M) '_' int2str(n_img) '.mat'];
if exist(overseg_data,'file')
    load(overseg_data);
else
    histSP_t = cell(n_img,1);
    labels_t = cell(n_img,1); colors_s_t = cell(n_img,1);
    lab_colors_s_t = cell(n_img,1); edges_s_t = cell(n_img,1);
    seg_t = cell(n_img,1); d_edges_t = cell(n_img,1);
    tic
    for i = 1:n_img
        [ histSP_t{i} labels_t{i} colors_s_t{i}  lab_colors_s_t{i} edges_s_t{i} seg_t{i} d_edges_t{i}] = ...
          over_segmentation( img_path, out_path, img_names_t{i}, nbins, hs, hr, M,full_connect,imgstyle);
    end
    toc
    save(overseg_data,'histSP_t', 'labels_t', 'colors_s_t',  'lab_colors_s_t', 'edges_s_t', 'seg_t','d_edges_t');
end
%% resort image name ,scribbled image in the front
un_scri_img_idx = ones(n_img,1);%index of un-scribbles image
un_scri_img_idx(scri_img_idx) = 0;
un_scri_img_idx = find(un_scri_img_idx);
re_img_idx = [scri_img_idx;un_scri_img_idx];
img_names = cell(n_img,1); histSP = cell(n_img,1);
labels = cell(n_img,1); colors_s = cell(n_img,1);
lab_colors_s = cell(n_img,1); edges_s = cell(n_img,1);
seg = cell(n_img,1); d_edges = cell(n_img,1);
for i  = 1:n_img
    img_names{i} = img_names_t{re_img_idx(i)};
    histSP{i} = histSP_t{re_img_idx(i)};
    labels{i} = labels_t{re_img_idx(i)};
    colors_s{i} = colors_s_t{re_img_idx(i)};
    lab_colors_s{i} = lab_colors_s_t{re_img_idx(i)};
    edges_s{i} = edges_s_t{re_img_idx(i)};
    seg{i} = seg_t{re_img_idx(i)};
    d_edges{i} = d_edges_t{re_img_idx(i)};
end
clear un_scri_img_idx  re_img_idx histSP_t labels_t colors_s_t lab_colors_s_t edges_s_t seg_t

%% get scribbles label of superpixel
superpixel_labelInd = cell(n_scri,1);
for i = 1:n_scri
    scribs_img_name = [scribbles_path scribbles_names{i} '.bmp'];
    [lines] = seed_generation(scribs_img_name,type_of_seed);%lines label index: one line represent one label index
    fg = unique(labels{i}(find(lines(:,1))));
    bg = unique(labels{i}(find(lines(:,2))));
    tmp1 = [fg;bg]; nf = length(fg);
    [b1 m1 n1] = unique(tmp1,'first');
    temp_label.fg = b1(m1<=nf); temp_label.bg = b1(m1>nf);
    superpixel_labelInd{i} = temp_label;
end

fprintf('lsr_iCoseg_superpixel\n');

tic
st=clock;
% get the segmentation results for superpixels
[YS,testParam] = lsr_iCoseg_superpixel_hist2(histSP,superpixel_labelInd,edges_s,colors_s,param);

T = etime(clock,st);
toc
% save results for superpixels
    if ~exist([out_path, 'lsr_hist'],'file')
        mkdir([out_path, 'lsr_hist']);
    end

save([out_path, 'lsr_hist/T.mat'],'T');
test_label_global = testParam.test_label_global;
testGMMprob = testParam.testGMMprob;

% transform the superpixels results to pixels results 
for i = 1:n_img
    YS_temp = YS{i};
    [h,w,d] = size(labels{i});
    ind = find(YS_temp);
    Y = zeros(h,w);
    for j = 1:length(ind)
        Y(labels{i} == ind(j))=1;
    end
    Y_mask = uint8(repmat(Y,[1,1,3]));
    im = imread([img_path ,img_names{i},'.' imgstyle]);

   [imgMasks,segOutline,imgMarkup]=segoutput(im2double(im),double(Y+1));

%     figure;
%     imshow(Y, []);
    % save pixels results
    file_save = [out_path, 'lsr_hist/',img_names{i}, '_segmentation.bmp'];
    file_save_img = [out_path, 'lsr_hist/',img_names{i}, '_mask_img.bmp'];
    if ~exist([out_path, 'lsr_hist'],'file')
        mkdir([out_path, 'lsr_hist']);
    end
    imwrite(Y,  file_save);
%     imwrite(im,  file_save_img);
    imwrite(imgMarkup,  file_save_img);
    
    %%%result of png
%     [token,remain] = strtok(dataset,'\');
%     out_path_png = ['ObjectDiscovery\Results\',token,'\Ours',remain];
%         if ~exist(out_path_png,'file')
%             mkdir(out_path_png);
%         end
%     file_save = [out_path_png,img_names{i},'.png'];
%     imwrite(Y,  file_save);
       
    
    %% iCoseg
%     groundtruth_path = [img_path,'GroundTruth\'];
%     gtImage = imread([groundtruth_path,img_names{i},'.png']);
%     groundtruth = double(gtImage(:,:,1)>0);
%     P(i) =sum(groundtruth(:)==Y(:)) ./ prod(size(groundtruth));
%     Jar(i) =sum( (Y(:)==1) & (groundtruth(:)==1) ) ./ sum( (Y(:) | groundtruth(:))==1 );
    
    %% MSRC
    groundtruth_path = [img_path,'GroundTruth\'];
    gtImage = imread([groundtruth_path,img_names{i},'.bmp']);
    groundtruth = double(gtImage(:,:,1))./255;
    P(i) =sum(groundtruth(:)==Y(:)) ./ prod(size(groundtruth));
    Jar(i) =sum( (Y(:)==1) & (groundtruth(:)==1) ) ./ sum( (Y(:) | groundtruth(:))==1 );
end
MP = mean(P); MJ = mean(Jar);
save ([out_path,'lsr_hist/PJ'] ,'P','Jar','MP','MJ','img_names');
fprintf('P=%f\nJ=%f\n',MP,MJ);



