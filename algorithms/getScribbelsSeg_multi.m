function [linesY] = getScribbelsSeg_multi(linesX,labels,n_Y)

n_layers = size(labels,2);

linesY = [];
for k=1:n_layers
    local_lines = zeros(n_Y(k),size(linesX,2));
    for i=1:size(linesX,2)
        idx = labels{k}(find(linesX(:,i)==1));
        local_lines(idx,i) = 1;
        clear idx;
    end;
    linesY = [linesY ; local_lines];
    clear local_lines;
end;
idx = find(sum(linesY,2) > 1);
linesY(idx,:) = 0;
