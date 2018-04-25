function [SP] = slic2_5D(imD, imN, SPwidth, m)
    % Implementation based on Achanta, Shaji, Smith, Lucchi, Fua & Susstrunk
    % SLIC Superpixels, with modifications depth images
    
    %Requires:
    % imD: an MxN raster comprising depth values rescaled to x/y pixel scale
    % imN: an MxNx3 raster comprising normals. Background should be zeros
    % SPwidth: the width in pixels of an expected incision (user-supplied)
    % m: optional parameter tuning between distance and normals weighting
    fprintf('SLIC ');
    delchar = fprintf('superpixels: Smoothing normals... ');
    maxIter = 25; %Most superpixels will be fully resolved within 5 iterations
    imgSz = size(imD);
    imN = squeeze(imN(:,:,3));   %Using just the Z-component of the normal
    imBG = imN<eps;                %Zero out negative normals (facing away)
    imN = wlsfilter(imN,0.25,1); %Gentle smoothing by weighted least squares
    imN(imBG) = 0;               
    
    %Seeds grid of superpixels according to pixel width of incisions
    step = ceil(0.75*SPwidth);
    if step < 5, step =5; end
    strips = ceil(imgSz/step);
    clXCs = linspace(step/2,imgSz(1)-step/2,strips(1));
    clXCs = repmat(ceil(clXCs'),1,strips(2));
    clYCs = linspace(step/2,imgSz(2)-step/2, strips(2));
    clYCs = repmat(ceil(clYCs),strips(1),1);
    
    clK = prod(strips);             %Number of clusters (k)
    Zscale = SPwidth/3;             %Increase impact of depth differences
    Zscale = max(min(Zscale,4),2);  %Constrained within range of 2-4
    imD = Zscale*imD;
    if nargin < 4, m = 0.25; end    %Default m at 0.25
    Dscale = m^2*clK/prod(imgSz)/Zscale; %Tuning distance function to 1*m
    
    %Initialize clusters, labels, and distances
    clMDs = zeros(strips);        %median cluster depth
    clMNs = zeros([strips,3]);    %median normal
    imLbl = -ones(imgSz);         %label by cluster number
    imDists = inf(imgSz);         %pixel distances from cluster centers
    for clX = 1:strips(1) %Initialize cluster values from 3x3 median values
        for clY = 1:strips(2)
            xR = clXCs(clX,clY)-1:clXCs(clX,clY)+1;
            yR = clYCs(clX,clY)-1:clYCs(clX,clY)+1;
            clMDs(clX,clY) = mean(reshape(imD(xR,yR),[],1));
            clMNs(clX,clY) = mean(reshape(imN(xR,yR),[],1));
        end
    end
        
    fprintf(repmat(char(8),delchar,1));
    delchar = fprintf('superpixels: Processing clusters... ');
    progT = progressTimer(1,0);
    iter = 0;
    updL = true(strips); %Marks which clusters to update in XY grid
    while iter < maxIter && sum(updL(:))
        iter = iter+1; oldLbl = imLbl;
        for clX = 1:strips(1)
            clXmin = max(min(round(clXCs(clX,:)))-step,1);
            clXmax = min(max(round(clXCs(clX,:)))+step,imgSz(1));            
            clXsubs = transpose((clXmin:clXmax));
            clXlen = clXmax-clXmin+1;
            
            for clY = find(updL(clX,:)) %Only process altered clusters
                clYmin = max(round(clYCs(clX,clY))-step,1);
                clYmax = min(round(clYCs(clX,clY))+step,imgSz(2));
                clYlen = clYmax-clYmin+1; %Subimage around cluster:
                clImD = reshape(imD(clXmin:clXmax,clYmin:clYmax),[],1);
                clImN = reshape(imN(clXmin:clXmax,clYmin:clYmax),[],1);
                
                %3D distance function
                clYlocal = repmat(clYmin:clYmax,clXlen,1);
                clXlocal = reshape(repmat(clXsubs,1,clYlen),[],1);
                cldiffXYZ = (reshape(clYlocal,[],1)-clYCs(clX,clY)).^2 + ...
                    (clXlocal-clXCs(clX,clY)).^2 + (clImD-clMDs(clX,clY)).^2;
                %Distance function based on the Z component of the normal
                cldiffN = abs(clImN-clMNs(clX,clY));
                
                clminDist = reshape(cldiffN + cldiffXYZ*Dscale,clXlen,clYlen);
                %Save all pixels whose distance is closer to this cluster
                %and assign the cluster number
                prevD = imDists(clXmin:clXmax,clYmin:clYmax);
                clUpdateIDX = find(clminDist < prevD);
                if ~isempty(clUpdateIDX)
                    prevL = imLbl(clXmin:clXmax,clYmin:clYmax);
                    prevL(clUpdateIDX) = clX + (clY-1)*strips(1);
                    prevD(clUpdateIDX) = clminDist(clUpdateIDX);

                    imDists(clXmin:clXmax,clYmin:clYmax) = prevD;
                    imLbl(clXmin:clXmax,clYmin:clYmax) = prevL;
                end
            end
        end
        [clMDs,clMNs,clXCs,clYCs,clIDX] = updateSP(imLbl,imD,imN,strips);
        changedL = unique(imLbl(imLbl~=oldLbl));
        changePct = length(changedL)/prod(strips);
        if changePct > 0.01
            updL = false(strips); updL(changedL) = true;
            prcCh = 2^(1-changePct) - 1;
            if prcCh<0.1, progT.update(0.1); else, progT.update(prcCh); end
        else, updL = false; %Conclude when only 1% of clusters may change
        end
    end
    progT.done;
    fprintf(repmat(char(8),delchar,1));
    
    delchar = fprintf('superpixels: Enforcing connectivity');
    %Split up non-continuous regions
    Lblnew = unique(imLbl(:))'; maxlabel = max(Lblnew);
    for cl = Lblnew
        if length(clIDX{cl}) < 2, continue; end
        clsubs = [];
        [clsubs(:,1),clsubs(:,2)] = ind2sub(imgSz,clIDX{cl});
        clmin = min(clsubs); clmax = max(clsubs);
        climlbl = imLbl(clmin(1):clmax(1),clmin(2):clmax(2));
        [bl,blnum] = bwlabel(climlbl==cl, 4);
        if blnum > 1  %More than one region has the same label
            %Keep label of the largest component, relabelling the smaller
            blct = zeros(blnum,1);
            for bn = 1:blnum, blct(bn) = sum(bl(:)==bn); end
            [~,blasc] = sort(blct);
            for bn = blasc(1:end-1)'
                maxlabel = maxlabel+1; %Generate a new label
                [blsubX,blsubY] = ind2sub(size(bl),find(bl==bn));
                blsubX = blsubX + clmin(1) - 1;
                blsubY = blsubY + clmin(2) - 1;
                imLbl(sub2ind(imgSz,blsubX,blsubY)) = maxlabel;
            end
        end
    end
    [clMDs,clMNs,clXCs,clYCs,clIDX] = updateSP(imLbl,imD,imN);
    
    %Merge the small clusters (many just split apart) with closest neighbor
    smallclustersExist = true;
    while smallclustersExist
        imNeighbors = getNeighbors(imLbl); %Identify neighboring regions
        mergecl = cellfun(@length,clIDX);
        mergecl = find(mergecl > 0 & mergecl < SPwidth); %Ignore empty clusters
        if isempty(mergecl), smallclustersExist = false; end
        for i = 1:length(mergecl)
            cl = mergecl(i);
            adjcl = find(imNeighbors(:,1)==cl);
            adjcl = imNeighbors(adjcl,2)'; %#ok<FNDSB>
            if length(adjcl) > 1 %Merge with the most similar cluster
                cldiffXYZ = (clXCs(cl)-clXCs(adjcl)).^2 + ...
                    (clYCs(cl)-clYCs(adjcl)).^2 + (clMDs(cl)-clMDs(adjcl)).^2;
                cldiffN = (clMNs(cl)-clMNs(adjcl)).^2;
                [~,clminDist] = min(cldiffN+cldiffXYZ*Dscale);
                adjcl = adjcl(clminDist);
            end
            imLbl(clIDX{cl}) = adjcl;
            clIDX{adjcl} = [clIDX{adjcl};clIDX{cl}];
            clIDX{cl} = [];
        end
    end
    
    %Merge all the background clusters, whose normal has been set to [0,0,0]
    bgLabels = find(clMNs==0)';
    if ~isempty(bgLabels)
        bgnum = bgLabels(1);
        clIDX{bgnum} = cell2mat(clIDX(bgLabels));
        for cl = bgLabels(2:end)
            clIDX{cl} = [];
            clXCs(cl) = Inf;
        end
        imLbl(clIDX{bgnum}) = bgnum;
        clXCs(bgnum) = prod(imgSz); %The background will be moved to the end
        
        [bl,blnum] = bwlabel(imLbl==bgnum, 4);
        if blnum > 1  %Split up non-continguous background areas
            maxlabel = max(imLbl(:));
            for bn = 1:blnum
                maxlabel = maxlabel+1;
                imLbl(bl==bn) = maxlabel;
                clXCs(maxlabel) = prod(imgSz)+maxlabel;
                clYCs(maxlabel) = prod(imgSz)+maxlabel;
            end
        end
    end
    
    %Resort labels by distances from origin, removing deleted values
    Lblold = unique(imLbl(:));
    [~,LIDXnew] = sort(clXCs(Lblold).^2+clYCs(Lblold).^2);
    Lblnew = zeros(length(Lblold),1);
    for i=1:length(Lblnew), Lblnew(i) = find(LIDXnew==i); end
    Lblremap = zeros(max(imLbl(:)),1);
    for i=1:length(Lblremap)
        remap = find(Lblold==i);
        if remap>0, Lblremap(i) = Lblnew(remap); end
    end    
    for x=1:imgSz(1)
        for y=1:imgSz(2), imLbl(x,y) = Lblremap(imLbl(x,y)); end
    end
    [~,clMNs,~,~,clIDX] = updateSP(imLbl,imD,imN);
    fprintf(repmat(char(8),delchar,1));
    
    SP = struct;
    SP.background = find(clMNs==0); %Save indices of background segments
    SP.labels = imLbl;
    SP.pixelCoords = clIDX;
    fprintf('segmentation complete. ');
end

function [mDs,mNs,cXs,cYs,clidx] = updateSP(imgl,ds,ns,strp)
    %Update superpixels with median values of depth / centroid for x,y
    clmax = max(imgl(:));
    sz = size(imgl);
    clidx = accumarray(imgl(:),1:prod(sz),[clmax 1],@(n) {sort(n)});
    
    mDs = zeros(clmax,1); mNs = mDs; cXs = mDs; cYs = mDs;
    for cl = 1:clmax
        imidx = clidx{cl};
        if ~isempty(imidx)
            mDs(cl) = median(ds(imidx));
            mNs(cl) = median(ns(imidx));
            [xidx,yidx] = ind2sub(sz,imidx);
            cXs(cl) = mean(xidx);
            cYs(cl) = mean(yidx);
        end
    end
    if nargin > 3 %Reshape for processing by seeds in a grid
        mDs = reshape(mDs,strp); mNs = reshape(mNs,strp);
        cXs = reshape(cXs,strp); cYs = reshape(cYs,strp);
    end
end

function lpairs = getNeighbors(iml)
    %Fast 4-connected neighbor function based on diff
    sz = size(iml);
    [dxy(:,1),dxy(:,2)] = find(diff(iml,1,1));
    dxy(dxy(:,2)>sz(2)-1,:) = [];
    dxyIDX = sub2ind(sz,dxy(:,1),dxy(:,2));
    dxyIDX2 = sub2ind(sz,dxy(:,1)+1,dxy(:,2));
    dxp1 = unique([iml(dxyIDX),iml(dxyIDX2)],'rows');
    
    dxy = [];
    [dxy(:,1),dxy(:,2)] = find(diff(iml,1,2));
    dxy(dxy(:,1)>sz(1)-1,:) = [];
    dxyIDX = sub2ind(sz,dxy(:,1),dxy(:,2));
    dxyIDX2 = sub2ind(sz,dxy(:,1),dxy(:,2)+1);
    dxp2 = unique([iml(dxyIDX),iml(dxyIDX2)],'rows');
    
    lpairs = unique([dxp1;dxp2;fliplr(dxp1);fliplr(dxp2)],'rows');
end

function imgOut = wlsfilter(imgIn, lambda, alpha)
%WLSFILTER Edge-preserving smoothing based on the weighted least squares(WLS) 
%   optimization framework, as described in Farbman, Fattal, Lischinski, and
%   Szeliski, "Edge-Preserving Decompositions for Multi-Scale Tone and Detail
%   Manipulation", ACM Transactions on Graphics, 27(3), August 2008.    
    [imgrows,imgcols] = size(imgIn);
    smallNum = 0.0001;
    np = imgrows*imgcols;
    L = log(imgIn+eps);
    % Compute affinities between adjacent pixels based on gradients of L
    dx = diff(L, 1, 2); dy = diff(L, 1, 1);
    dx = -lambda./(abs(dx).^alpha + smallNum);
    dy = -lambda./(abs(dy).^alpha + smallNum);
    dx = padarray(dx, [0 1], 'post'); dy = padarray(dy, [1 0], 'post'); 
    dx = dx(:); dy = dy(:);
    % Construct a five-point spatially inhomogeneous Laplacian matrix
    B(:,1) = dx; B(:,2) = dy;
    d = [-imgrows,-1];
    A = spdiags(B, d, np,np);
    w = padarray(dx, imgrows, 'pre'); w = w(1:end-imgrows);
    n = padarray(dy, 1, 'pre'); n = n(1:end-1);
    A = A + A' + spdiags(1-(dx+w+dy+n), 0, np,np);
    % Solve
    imgOut = A\reshape(imgIn,[],1);
    imgOut = reshape(imgOut, imgrows, imgcols);
end