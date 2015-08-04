function D = spblkdiag(C)
    % spblkdiag(C)
    %
    %   Efficiently concatenate a cell array of sparse matrices into
    %   one sparse block diagonal matrix.
    %
    % C: a cell array of diagonal images
    
    n = numel(C);
    
    % Count total nonzeros + allocate COO arrays
    Dnz = 0;
    
    for i=1:n
        Dnz = Dnz + nnz(C{i});
    end
    
    Di = zeros(Dnz, 1);
    Dj = Di;
    Dv = Di;
    
    % Fill COO arrays
    start = 0;
    nrows = 0;
    ncols = 0;
    
    for i=1:n
        % Unpack C matrix
        Cnz = nnz(C{i});
        [Ci Cj Cv] = find(C{i});
        
        Di(start+1:start+Cnz) = Ci + nrows;
        Dj(start+1:start+Cnz) = Cj + ncols;
        Dv(start+1:start+Cnz) = Cv;
        
        % Update indices
        start = start + Cnz;
        nrows = nrows + size(C{i}, 1);
        ncols = ncols + size(C{i}, 2);
    end
    
    assert(start == Dnz);
    
    % Construct block diagonal matrix
    D = sparse(Di, Dj, Dv, nrows, ncols);
end