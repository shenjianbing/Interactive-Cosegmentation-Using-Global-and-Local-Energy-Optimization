function save_segmentation_results(Posteriors, n_X, name, img_path, out_path, bg_sub)

img = im2double(imread([img_path name '.bmp']));    h = size(img,1);    w = size(img,2);
nlabels = size(Posteriors,2);

if nargin < 6, bg_sub = 1; end;

%% Display the result of the segmentation for pixels
data_path = [out_path 'segments/']; mkdir(data_path);
[vals,inds] = max(Posteriors(1:n_X,:)');
mask = reshape(inds,h,w);

[imgMasks,segOutline,imgMarkup]=segoutput_t(img,mask,bg_sub);
figure; imshow(imgMarkup);  

if nargin < 7
    imwrite(imgMarkup, [out_path 'segments/' name '_ours.bmp']);
    save([out_path 'segments/' name '_ours.mat'],'mask');
else
    imwrite(imgMarkup, [out_path 'segments/' name '_ours.bmp']);
    save([out_path 'segments/' name '_ours.mat'],'mask');
end;
clear imgMasks segOutline imgMarkup;

disp_colors(1,1:3) = [1,0,0];
disp_colors(2,1:3) = [0,1,0];
disp_colors(3,1:3) = [0,0,1];
disp_colors(4,1:3) = [1,1,0];
disp_colors(5,1:3) = [1,0,1];
disp_colors(6,1:3) = [0,1,1];
disp_colors(7,1:3) = [1,1,1];
disp_colors(8,1:3) = [0,0,0];

data_path = [out_path 'labels/']; mkdir(data_path);
if nlabels <= 8,
    label_img = zeros(h,w,3);
    for nc=1:3,
        tmp = label_img(:,:,nc);
        for i=1:nlabels, tmp(find(mask==i)) = disp_colors(i,nc); end;
        label_img(:,:,nc) = tmp; clear tmp;
    end;
    imwrite(label_img, [out_path 'labels/' name '_seg.bmp']);
end;

