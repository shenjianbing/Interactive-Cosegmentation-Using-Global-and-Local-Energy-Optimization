function [ histSP labels colors_s  lab_colors_s edges_s seg d_edges ] = over_segmentation( img_path, out_path, img_name, nbins, hs, hr, M, full_connect,imgstyle)
%OVER_SEGMENTATION Summary of this function goes here
%   Detailed explanation goes here

% get superpixels and their histogram feature


%% information of pixels
[img lab_img points_p edges_p colors_p lab_colors_p] = getPropertiesForPixels([img_path img_name '.' imgstyle]);
img_v = double(reshape(img,[],3));

%% over-segmentation
data_path = [out_path 'regions/']; 
if ~exist(data_path)
    mkdir(data_path);
end

    [segs labels seg  colors_s  lab_colors_s edges_s d_edges ] = msseg(double(img),reshape(lab_img, size(img,1)*size(img,2), size(img,3)),hs,hr,M,full_connect);
    [imgMasks,segOutline,imgMarkup]=segoutput(im2double(img),double(labels)); 
    imwrite(imgMarkup, [out_path 'regions/' img_name '_' int2str(hs) '_' int2str(hr) '_' int2str(M) '.bmp']);  
    [h,w] = size(labels);
    imgColors = zeros(w*h,3);
    temp_lables = reshape(labels,w*h,1);
    for i = 1:size(colors_s,1)
        imgColors(temp_lables ==i,1) = colors_s(i,1);
        imgColors(temp_lables ==i,2) = colors_s(i,2);
        imgColors(temp_lables ==i,3) = colors_s(i,3);
    end
    imgColors = reshape(imgColors,h,w,3)./255;
    imwrite(imgColors, [out_path 'regions/' img_name '_' int2str(hs) '_' int2str(hr) '_' int2str(M) 'colors.bmp']);  
    clear imgMasks segOutline imgMarkup;
    %% histogram of superpixel
    n_superpixel = length(seg);
    temp_hist = zeros(nbins*3,n_superpixel);
    for j =1:n_superpixel
        temp_index = seg{j};
        SP = img_v(temp_index,:);
        temp_hist(:,j)  = makehist( nbins,SP );
    end
    histSP = temp_hist;%d x n
    clear temp_hist
end

