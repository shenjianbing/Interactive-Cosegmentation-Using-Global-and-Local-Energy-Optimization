function [img lab_img points edges colors lab_colors] = getPropertiesForPixels(img_name)

img = imread(img_name); lab_img = colorspace('Lab<-', img);
h = size(img,1); w = size(img,2);
[points edges] = lattice(h,w,1);

if size(img,3) > 1
    tmp = img(:,:,1); colors(:,1) = tmp(:); clear tmp;
    tmp = img(:,:,2); colors(:,2) = tmp(:); clear tmp;
    tmp = img(:,:,3); colors(:,3) = tmp(:); clear tmp;
    
    tmp = lab_img(:,:,1); lab_colors(:,1) = tmp(:); clear tmp;
    tmp = lab_img(:,:,2); lab_colors(:,2) = tmp(:); clear tmp;
    tmp = lab_img(:,:,3); lab_colors(:,3) = tmp(:); clear tmp;
else
    colors = img(:);
    lab_colors = lab_img(:);
end;