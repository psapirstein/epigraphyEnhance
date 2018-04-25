function [RGout] = slicMerge(RG, imD, imN)
    % Superpixel merging algorithm, related to the approach in:
    % M.S. Chaibou, P.-H. Conze, K. Kalti, B. Solaiman, & M.A. Mahjoub
    %   "Adaptive strategy for superpixel-based region-growing image
    %    segmentation" Journal of Electronic Imaging 26.6, 061605 (2017)
    fprintf('Merging superpixels');
    delchars = fprintf(': initializing regions ');
    
    %Superpixel linkages (matrix including a to b, and b to a):
    maxlbl = max(RG.labels(:));
    LS = getLinks(RG.labels,RG.background); %Link matrix A,B,score
    RG.regions = unique(LS(:,1)); %Superpixels that seed initial regions
    RG.SPlist = num2cell(transpose(1:maxlbl)); %Superpixels in each region
    RG.silhouettes = cell(maxlbl,1); %List of pixels encircling each region
    RG.silhouettes(RG.regions) = getPerimeter(RG, RG.regions);
    %Initialize scores for adjacent superpixels
    [RG.planes,RG.points,RG.normals] = initializeRegions(RG,imD,imN);
    LS = scoreLinks(RG,LS);
    %Save initial similarity scores for lookup during later border tests
    RG.SPborders = accumarray(LS(:,1),LS(:,2),[maxlbl 1],@(n){n});
    RG.SPvalues = sparse(LS(:,1),LS(:,2),LS(:,3),maxlbl,maxlbl);
    
    fprintf(repmat(char(8),delchars,1));
    delchars = fprintf(': processing ');
    
    adaptSim = 2;    %Initial threshold
    Ait = 1.05;      %Tunes the threshold for merging adjacent regions
    stopSim = 0.9;   %Stop merging below this threshold score
    simRange = adaptSim-stopSim;
    mergePct = 0.01; %Minimum percentage required before merging
    its = 0;
    
    nRegions = length(RG.regions);
    [mergeList,mergeScore] = compareSimilarity(RG,LS);
    
    progT = progressTimer(1,0);
    progMin = 0.01;
    while adaptSim > stopSim
        nMerged = sum(mergeScore > adaptSim);
        if nMerged < nRegions*mergePct
            adaptSim = adaptSim/Ait; %Lower the threshold for merging
        else
            mergeList = mergeList(mergeScore > adaptSim,:);
            Ait = max(1+nMerged/nRegions, 1.01);
            adaptSim = adaptSim*Ait;
            [RG,LS,newLS] = mergeRegions(RG,mergeList,LS);
            LS = scoreLinks(RG,LS,newLS);
            
            its = its+1;
            nRegions = length(RG.regions);
            [mergeList,mergeScore] = compareSimilarity(RG,LS);
            progMin = max((1-(adaptSim-stopSim)/simRange)^1.5,progMin);
            progT.update(progMin);
        end
    end
    progT.done;
    fprintf(repmat(char(8),delchars+its,1));
    fprintf(' complete.\n');
    rmfnames = {'points','normals','silhouettes',...
        'SPborders','SPlist','SPvalues'};
    RGout = rmfield(RG,rmfnames);
    RGout.planes = RG.planes(:,1:6);
    RGout.segmentCurvature = RG.planes(:,8);
end

function [LS] = scoreLinks(RG, LSin, lnkidx)
    LS = LSin; initializing = false;
    fullA = LS(:,1); fullB = LS(:,2);
    if nargin < 3 %Process everything (initialization)
        procA = fullA; procB = fullB;
        lnkidx = 1:length(procA);
        initializing = true;
    else          %Process the supplied indices of the link table
        procA = LS(lnkidx,1); procB = LS(lnkidx,2);
    end
    if size(LS,2) < 3, LS = [LS, -ones(size(LS,1),5)]; end
    
    previdx = 0;
    for il = 1:length(procA)
        aidx = procA(il);
        if aidx > previdx
            Actr = RG.planes(aidx,1:3);
            Anorm = RG.planes(aidx,4:6)';
            previdx = aidx;
        end
        %Edge weights: median angle of B points relative to the plane of A
        Bpts = RG.points{procB(il)}-Actr; BptD = Bpts*Anorm;
        LS(lnkidx(il),4) = median(BptD./sqrt(sum(Bpts.^2,2)));
        LS(lnkidx(il),7) = std(BptD);
    end
    pointAngle = asind(min(max(LS(:,4),-1),1));
    planeAngle = acosd(sum(RG.planes(fullA,4:6).*RG.planes(fullB,4:6),2));
    %Rescale to -1 to 1, clipped to 45 degrees (maximum angle tested)
    pointAngle = max(min(pointAngle/45,1),-1);
    planeAngle = max(min(planeAngle/45,1),-1);
    %Reduce the angle between close scores
    rescale = @(x,n,y) x + x.*sqrt(1-(abs(x)-y)/(1-y))*(n-y)/y;
    planeAngle = rescale(planeAngle,0.1,0.2);
    
    cFac = 2; %Lowering matches across convex edges / raising concave
    convexEdge = pointAngle < 0;
    pointAngle(~convexEdge) = pointAngle(~convexEdge).^cFac;
    pointAngle(convexEdge) = abs(pointAngle(convexEdge)).^(1/cFac);
    angleScore = (1-pointAngle).*(1-planeAngle);
    %Effect on standard deviation of B by fitting to the plane for A
    %Anything above 1.5 receives the minimum score (~4.5 times worse)
    stdScore = max(log(LS(:,7)./max(RG.planes(fullB,7),0.001)),0)/1.5;
    stdScore = max(1-stdScore,0);
    %Very different curvatures have a significant impact
    curvScore = RG.planes(fullA,8)-RG.planes(fullB,8);
    curvScore = 1-min(abs(curvScore)/(3*1.4826*mad(curvScore)),1);
    
    %Region-to-region matching depends primarily on angle scores
    simRegion = (angleScore.*curvScore*3 + stdScore)/4;
    
    %Initialize the region boundary scores with the region scores
    if initializing, LS(:,5) = simRegion; end
    regSPnum = cellfun(@length,RG.SPlist);
    procRegion = find( regSPnum(procA) > 1 | regSPnum(procB) > 1 );
    numSP = length(RG.pixelCoords);
    for il = 1:length(procRegion)
        ilA = procA(procRegion(il)); ilB = procB(procRegion(il));
        %Superpixels from the 2nd region bordering the 1st, & vice-versa        
        ilspA = RG.SPborders{ilA}; ilspB = RG.SPborders{ilB};
        ilspA = ilspA(ismembc(ilspA,RG.SPlist{ilB}));
        ilspB = ilspB(ismembc(ilspB,RG.SPlist{ilA}));
        %Lookup original similarity scores of the border superpixels
        [ilspA,ilspB] = meshgrid(ilspA,ilspB);
        spvals = nonzeros(RG.SPvalues(ilspA(:) + (ilspB(:)-1)*numSP));
        if isempty(spvals), spvals = 0; end
        LS(lnkidx(procRegion(il)),5) = sum(spvals(:))/numel(ilspA);
    end
    %Weights for border scores from the length of the shared perimeter
    rgcircum = cellfun(@length,RG.silhouettes);
    crcmA = rgcircum(procA); crcmB = rgcircum(procB);
    wtBeta = ones(length(procA),1);
    previdx = 0;
    for il = 1:length(procA)
        aidx = procA(il);
        if aidx > previdx
            pxA = RG.silhouettes{aidx};
            previdx = aidx;
        end
        wtBeta(il) = sum(ismembc(pxA,RG.pixelCoords{procB(il)}));
    end
    wtBeta = min( sqrt(wtBeta .* (crcmA+crcmB)./(crcmA.*crcmB)), 1);
    LS(lnkidx,6) = wtBeta;
    LS(:,3) = max(simRegion + LS(:,6).*LS(:,5), 1.0e-6);
end

function [mergeList,mergeScores] = compareSimilarity(RG,LS)
    linkIDX = [1; 1+find(diff(LS(:,1))>0); size(LS,1)+1];
    rgnum = length(RG.pixelCoords);
    bestMatch = zeros(rgnum,1);
    bestScore = bestMatch;
    for il=1:length(RG.regions)
        linkRange = linkIDX(il):(linkIDX(il+1)-1);
        [bestScore(RG.regions(il)),nearidx] = max(LS(linkRange,3));
        bestMatch(RG.regions(il)) = LS(linkRange(nearidx),2);
    end
    
    hasMatch = find(bestMatch);
    penalty = 0.90; %Non-mutual match scores are penalized
    nonmutualMatch = bestMatch(bestMatch(hasMatch)) ~= hasMatch;
    bestScore(nonmutualMatch) = penalty * bestScore(nonmutualMatch);
    %Combine the links and scores
    mergeList = [transpose(1:rgnum),bestMatch,bestScore];
    mergeList = mergeList(hasMatch,:);
    %Duplicate the rows without a reciprocal match
    mergeList = [mergeList; mergeList(nonmutualMatch,[2 1 3])];
    mergeList(:,1:2) = sort(mergeList(:,1:2),2);
    mergeList = sortrows(mergeList);
    %Reduce to single links, A < B, with averaged scores
    mergeList(1:2:end-1,3) = (mergeList(1:2:end-1,3) + mergeList(2:2:end,3))/2;
    mergeList = mergeList(1:2:end-1,:);
    
    %Remove links that do not have the best score for their region
    bestScore = accumarray(mergeList(:,1),mergeList(:,3),[rgnum 1],@max);
    mergeList(mergeList(:,3) < bestScore(mergeList(:,2)),:) = [];
    mergeList = sortrows(mergeList, [1 3], {'ascend' 'descend'});
    mergeList = mergeList( [1; 1+find(diff(mergeList(:,1))>0)], :);
    %Remove duplicate second regions with lower scores
    mergeList = sortrows(mergeList, [2 3], {'ascend' 'descend'});
    mergeList = mergeList( [1; 1+find(diff(mergeList(:,2))>0)], :);
    %If a second region matches a first, remove the second set (because of
    %the sorting, the second will have a lower score than the first)
    delIDX = ismember(mergeList(:,1),unique(mergeList(:,2)));
    mergeList(delIDX,:) = [];
    mergeList = sortrows(mergeList);
    
    mergeScores = mergeList(:,3);
    mergeList = mergeList(:,1:2);
end

function [RG,newLS,newidx] = mergeRegions(RG,mL,oldLS)
    mLA = mL(:,1);
    nparents = length(mLA);
    newSPlist = cell(nparents,1);
    ncrd = newSPlist; nPts = newSPlist; nNs = newSPlist; nPrm = newSPlist;
    newSPborder = newSPlist;
    areaRatio = cellfun(@length,RG.pixelCoords(mL));
    areaRatio = areaRatio(:,1)./sum(areaRatio,2);
    for is=1:length(mLA)
        parentR = mLA(is); childR = mL(is,2);
        childPx = RG.pixelCoords{childR};
        %Add the child to the parent
        rgcontent = sort([RG.SPlist{parentR}; RG.SPlist{childR}]);
        rgint = sort([RG.pixelCoords{parentR}; childPx]);
        newSPlist{is} = rgcontent; ncrd{is} = rgint;
        nPts{is} = [RG.points{parentR}; RG.points{childR}];
        nNs{is} = [RG.normals{parentR}; RG.normals{childR}];
        RG.labels(childPx) = parentR;
        
        %Combine perimeter pixels, removing those from interior
        regperim = [RG.silhouettes{parentR}; RG.silhouettes{childR}];
        nPrm{is} = setdiff(regperim,rgint);
        %Concatenate border superpixels, removing any inside the new region
        bordRG = [RG.SPborders{parentR}; RG.SPborders{childR}];
        newSPborder{is} = setdiff(bordRG,rgcontent);
        %Update the plane data for each region
        newN = RG.planes([parentR;childR], 4:6);
        newN = areaRatio(is)*newN(1,:) + (1-areaRatio(is))*newN(2,:);
        RG.planes(mLA(is),:) = fitPlane(nPts{is},nNs{is}, newN/norm(newN));
        %Clear the discarded child region
        RG.SPlist{childR} = []; RG.pixelCoords{childR} = [];
        RG.silhouettes{childR} = []; RG.SPborders{childR} = [];
        RG.points{childR} = []; RG.normals{childR} = [];
    end
    RG.SPlist(mLA) = newSPlist; RG.pixelCoords(mLA) = ncrd;
    RG.normals(mLA) = nNs; RG.points(mLA) = nPts;
    RG.silhouettes(mLA) = nPrm; RG.SPborders(mLA) = newSPborder;
    RG.planes(mL(:,2),:) = zeros(nparents,8);
    %Regenerate the link table and active regions
    newLS = getLinks(RG.labels,RG.background);
    RG.regions = unique(newLS(:,1));
    %Copy the previous similarity scores for unchanged regions
    unchgidx = ismember(newLS(:,1:2),mLA); %All links with merged regions
    unchgidx = sum(unchgidx,2)<1; %The set with no merged regions
    unchgidxold = ismember(oldLS(:,1:2),mL(:));
    unchgidxold = sum(unchgidxold,2)<1;
    %Regions that must be recalculated are scored -1
    newLS = [ newLS, -ones(length(unchgidx),5) ]; %Padding, columns 3-7
    for ldim = 3:7, newLS(unchgidx,ldim) = oldLS(unchgidxold,ldim); end
    newidx = find(~unchgidx);
end

function [planedata] = fitPlane(pts,ptsN,plNin)
    %Weighted planar fitting to points, provided a seed normal
    np = size(pts,1);
    gctr = L1median(pts); %Geometric center approximation
    pctr = pts-repmat(gctr,np,1);
    %Initialize plane from unweighted points, then reweight iteratively
    [plNout,~] = eig(pctr'*pctr/(np-1));
    plNout = plNout(:,1);
    for its=1:3 %Weighted plane normal beginning with supplied vector
        wts = abs(pctr*plNout);
        wts = 1-wts/max(wts);
        if any(isnan(wts)), wts = ones(np,1); end
        ptw = repmat(wts,1,3).*pctr;
        %PCA / eigenvectors of covariance matrix from weighted points
        [plNout,~] = eig(ptw'*ptw/(np-1));
        plNout = plNout(:,1);
    end    
    if plNin*plNout<0, plNout = -plNout; end
    if plNout(3)<0, plNout = -plNout; end
    
    %Estimating curvature by the Chord-and-Normal method for each piont
    %See Zhang and Cheng 2008, Curvature Estimation of 3D Point Cloud
    %Surface Through the Fitting of the Normal Section Curvatures,
    %ASIAGRAPH 2008 Proceedings
    plPtvec = sqrt(sum(crossProdcols(plNout,pctr).^2,2));
    zerolen = find(plPtvec<eps);
    if ~isempty(zerolen)
        %A point in the scan coincides with the center of the fitted plane
        plPtvec(zerolen) = []; ptsN(zerolen,:) = []; pctr(zerolen,:) = [];
    end
    kGeomMed = sqrt(sum(crossProdcols(plNout,ptsN).^2,2))./plPtvec;
    %Returning the median signed curvature from the set of points
    kGeomMed = median(kGeomMed.*sign(sum(ptsN.*pctr,2)));
    planedata = [gctr,plNout',std(pctr*plNout),kGeomMed];
end

function [x] = crossProdcols(a,b)
    x = [a(2)*b(:,3)-a(3)*b(:,2), a(3)*b(:,1)-a(1)*b(:,3), ...
        a(1)*b(:,2)-a(2)*b(:,1)];
end

function [mctr] = L1median(pts)
    % Ref: Hossjer and Croux (1995) "Generalizing Univariate Signed Rank
    % Statistics for Testing and Estimating a Multivariate Location
    % Parameter", Non-parametric Statistics, 4, 293-308.
    % Translated from the Gauss code of Hossjer and Croux (1995) in Matlab
    % by Sabine Verboven, Antwerp University. Condensed by P. Sapirstein.
    mrobj = @(X,M) sum(sqrt(sum( (X-repmat(M,size(X,1),1)).^2, 2)));
    tol = 0.01; maxits = 20;
    np = size(pts,1);
    %initializing starting value for m
    mctr = median(pts,1); its = 1;
    while its <= maxits
        mold = mctr;
        ptctr = pts-repmat(mctr,np,1);
        wts = 1./(sqrt(sum(ptctr.^2,2))+1.0e-8);
        delta = sum(ptctr.*repmat(wts,1,3),1) / sum(wts);
        maxhalf = 0;
        nd = sqrt(sum(delta.^2,2));
        if nd >= tol, maxhalf = log2(nd/tol); end
        mctr = mold + delta;   %computation of a new estimate
        nstep = 0;
        while mrobj(pts,mctr) >= mrobj(pts,mold) && (nstep <= maxhalf)
            nstep = nstep+1;
            mctr = mold + delta./(2^nstep);
        end
        if (nstep >= maxhalf)
            mctr = mold;
            return;
        end
        its = its+1;
    end
end

function [lpairs] = getLinks(iml,lbg)
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
    %Delete all rows that contain an index from the background
    bgrows = [];
    for is=1:length(lbg)
        [ri,~] = find(lpairs == lbg(is));
        bgrows = unique([bgrows;ri]);
    end
    lpairs(bgrows,:) = [];
end

function [S] = getPerimeter(RG,indices)
    sz = size(RG.labels);
    S = cell(length(indices),1);
    for ir=1:length(indices)
        pIDX = RG.pixelCoords{indices(ir)};
        if ~isempty(pIDX)
            [bx,by] = ind2sub(sz,pIDX);
            xshift = min(bx)-2; yshift = min(by)-2;
            bx = bx-xshift; by=by-yshift;
            bimg = false(max(bx)+1,max(by)+1);
            bimg(sub2ind(size(bimg),bx,by)) = true;
            
            pimg = zeros(size(bimg));
            pimg(1:end-1,:) = diff(bimg,1,1)>0;
            pimg(2:end,:) = pimg(2:end,:) | diff(bimg,1,1)<0;
            pimg(:,1:end-1) = pimg(:,1:end-1) | diff(bimg,1,2)>0;
            pimg(:,2:end) = pimg(:,2:end) | diff(bimg,1,2)<0;
            [bx,by] = find(pimg);
            bx = bx+xshift; by = by+yshift;
            includexy = bx > 0 & by > 0 & bx <= sz(1) & by <= sz(2);
            
            S{ir} = sub2ind(sz, bx(includexy), by(includexy));
        end
    end
end

function [rpln,pts3D,ptsN] = initializeRegions(RG,imDs,imNs)
    imsz = size(RG.labels);
    improd = numel(RG.labels);
    rgidx = RG.regions';
    rpln = zeros(length(RG.SPlist),8); %1:3 Plane center, 4:6 normal
        %7 standard deviation from the fitted plane; 8: curvature estimate
    pts3D = cell(size(rpln,1),1);
    ptsN = pts3D;
    for ir=rgidx
        pIDX = RG.pixelCoords{ir}; pt3D = [];
        %3D points in pixel units for each superpixel
        [pt3D(:,1),pt3D(:,2)] = ind2sub(imsz,pIDX);
        pt3D(:,3) = reshape(imDs(pIDX),[],1);
        pt3D(:,2) = -pt3D(:,2); %Flip ys from image encoding
        
        rgN = reshape(imNs([pIDX, pIDX+improd, pIDX+2*improd]),[],3);
        ptMask = sum(rgN.^2,2) > eps; %Zeroed normals are masked out
        pt3D = pt3D(ptMask,:); rgN = rgN(ptMask,:);
        %Median normal from the data is used for weighting the plane normal
        rgNmed = median(rgN,1);
        rpln(ir,:) = fitPlane(pt3D,rgN,rgNmed/norm(rgNmed));
        pts3D{ir} = pt3D; ptsN{ir} = rgN;
    end
end