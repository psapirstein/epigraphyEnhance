function [nA] = normalizeCells(A)
    nA = A;
    for i=1:size(A,1)
        colA = squeeze(A(i,:,:));
        nTmp = sqrt( sum( colA.^2, 2 ) );
        nTmp(nTmp == 0) = 1;
        nA(i,:,:) = bsxfun(@rdivide, colA, nTmp);
    end
end