function [ histSP ] = makehist( nbins,SP )
%
% SP: n x 3 vector
% nbins: number of  bins (1D)
% pixel range: 0 to 255
% histSP: nbins*3 x 1 vector


[n,d] = size(SP);
histSP = zeros(d,nbins);
for i = 1:d
    histSP(i,:) = hist(SP(:,i),nbins);
end
histSP = reshape(histSP,[],1);

end

