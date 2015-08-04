function [W weights_x n_X n_Y] = getW_multi(imgValsX,imgValsY,edgesX,edgesY,seg,beta)

n_X = size(imgValsX,1);   % n_X: the number of pixels
n_layers = size(seg,2);

n_stack = n_X;
imgVals = imgValsX;
total_edges = edgesX;

for k=1:n_layers
    n_Y(k) = size(imgValsY{k},1);   % n_X: the number of pixels, n_Y: the number of regions
    imgVals     = [imgVals ; imgValsY{k}];
    % links between pixels and regions
    edgesXY = [];
    for i=1:n_Y(k)
        col1 = (i+n_stack)*ones(size(seg{k}{i},1),1);
        col2 = seg{k}{i};
        edgesXY = [edgesXY ; [col1 col2]];
        edgesXY = [edgesXY ; [col2 col1]];
        clear col1 col2;
    end;
    total_edges = [total_edges ; edgesY{k}+n_stack ; edgesXY];
    n_stack = n_stack + n_Y(k);
    clear edgesXY;
end;

weights=makeweights(double(total_edges),imgVals,beta);
weights_x=weights(1:size(edgesX,1),:);
W=adjacency(double(total_edges),weights,n_stack);

