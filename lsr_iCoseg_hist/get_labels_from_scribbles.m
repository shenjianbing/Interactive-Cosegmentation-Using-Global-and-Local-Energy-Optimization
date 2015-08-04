function [ label_ind  GMMprob] = get_labels_from_scribbles( feature,scribbled_sp_fea_all, scribbled_sp_all_ind,labelParam )
%GET_LABELS_FROM_SCRIBBLES Summary of this function goes here
%   Detailed explanation goes here
switch labelParam.fea_style
    case 'histogram'
        n_labels = 6;
        histSP = feature;
        scribbled_sp_hist_all = scribbled_sp_fea_all;
        hist_foregroud = sum(scribbled_sp_hist_all(:,scribbled_sp_all_ind.fg),2);%d x 1;
        hist_backgroud = sum(scribbled_sp_hist_all(:,scribbled_sp_all_ind.fg),2);%d x 1;
        %normalize 
        hist_foregroud =  hist_foregroud/sum(hist_foregroud);
        hist_backgroud =  hist_backgroud/sum(hist_backgroud);
        [d,n_sp] = size(histSP);
        temp = sum(histSP,1);
        temp = repmat(temp,d,1);
        histSP = histSP./temp;

        dist = sqdist(histSP,hist_foregroud);
        [~,index] = sort(dist,1);
        label_ind.fg = index(1:n_labels);
        label_ind.bg = index(end-n_labels+1:end);
        GMMprob =[];


    case 'colors'
        
        %%%% GMM
        n_labels = 6;
        colorsSP = feature;
        scribbled_sp_colors_all = scribbled_sp_fea_all;
        GMMparam = labelParam;
        A = colorsSP;
        unarypotb = GMMparam.unarypotB;        muB = GMMparam.muB;        sigmaB = GMMparam.sigmaB;
        unarypotF = GMMparam.unarypotF;        muF = GMMparam.muF;        sigmaF = GMMparam.sigmaF;    
        ncenters = GMMparam.ncenters;     ncentersF = GMMparam.ncentersF;   ncentersB = GMMparam.ncentersB;
        n = size(colorsSP,2);
        stabF = 0.01;
        stabB = 0.01;

        wo = zeros(n,ncentersB);
        wb = zeros(n,ncentersF);

        for i=1:ncentersB
            L = chol(inv(sigmaB(:,:,i) + stabB*eye(3)));
            tmp = L*(A - repmat(muB(:,i),1,n));
        %     mean( sum( (A - repmat(muB(:,i),1,n)).^2, 1 ))
            wo(:,i) = -log(unarypotb(i)) + 0.5*log(det(sigmaB(:,:,i) + stabB*eye(3))) + 0.5* sum( tmp.^2, 1)'; 
        end
        
        for i=1:ncentersF
            L = chol(inv(sigmaF(:,:,i)+stabF*eye(3)));
            tmp = L*(A - repmat(muF(:,i),1,n));
            wb(:,i) = -log(unarypotF(i)) + 0.5*log(det(sigmaF(:,:,i) + stabF*eye(3))) + 0.5* sum( tmp.^2, 1)'; 
        end

        wo2 = min(wo, [], 2) + 3/2 * log(2*pi);
        wb2 = min(wb, [], 2) + 3/2 * log(2*pi);

        while min(wo2) < 0
            stabB = stabB + 0.1;
            for i=1:ncentersB
                L = chol(inv(sigmaB(:,:,i) + stabB*eye(3)));
                tmp = L*(A - repmat(muB(:,i),1,n));
                wo(:,i) = -log(unarypotb(i)) + 0.5*log(det(sigmaB(:,:,i) + stabB*eye(3))) + 0.5* sum( tmp.^2, 1)';
            end
            wo2 = min(wo, [], 2) + 3/2 * log(2*pi);
        end
        while min(wb2) < 0
            stabF = stabF + 0.1;
            for i=1:ncentersF
                L = chol(inv(sigmaF(:,:,i)+stabF*eye(3)));
                tmp = L*(A - repmat(muF(:,i),1,n));
                wb(:,i) = -log(unarypotF(i)) + 0.5*log(det(sigmaF(:,:,i) + stabF*eye(3))) + 0.5* sum( tmp.^2, 1)';
            end
            wb2 = min(wb, [], 2) + 3/2 * log(2*pi);        
        end

        nk =6;%ceil(0.01*n);%
        posterior = wb2./(wb2+wo2);
        [~,index] = sort(posterior,1);
        GMMprob.posteriorF = 1-posterior;
        GMMprob.likeliF = wb2;
        nt = sum((1-posterior)>0.5);
        if nt<nk
            nk=nt;
        end
        label_ind.fg = reshape(index(1:nk,:),[],1);
        

        posterior = wo2./(wb2+wo2);
        [~,index] = sort(posterior,1);
        GMMprob.posteriorB = 1-posterior;
        GMMprob.likeliB = wb2;
        nt = sum((1-posterior)>0.5);
        if nt<nk
            nk=nt;
        end
        label_ind.bg = reshape(index(1:nk,:),[],1);    
    otherwise
        
end
end

