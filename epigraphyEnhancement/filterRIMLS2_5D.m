function [imgout] = filterRIMLS2_5D(imgin,ptcloud)
    %Set up the parallel computing
    numWorkers = 7; %Suited to an 8-core machine
    paralleljob = gcp('nocreate');
    if isempty(paralleljob), paralleljob = parpool('local8',numWorkers);
    elseif paralleljob.NumWorkers ~= numWorkers
        delete(paralleljob);
        paralleljob = parpool('local8',numWorkers);
    end
    fprintf('Refining surface model by RIMLS: ');
    %Initialize variables
    Pts = double(ptcloud.Location);
    Ns = double(ptcloud.Normal);
    [~,procIDX] = sort(Pts(:,1));
    Pts = Pts(procIDX,:); Ns = Ns(procIDX,:);
    
    imgout = imgin;
    clear ptcloud imgin;
    
    imgSz = size(imgout.Mask);
    pixW = imgout.PixD;
    minDepth = min(imgout.Depths(:));
    
    imgCs = [imgout.xCoords(:),imgout.yCoords(:),imgout.Depths(:)];
    imgNs = imgCs;
    for k=1:3
        imgNs(:,k) = reshape(imgout.Normals(:,:,k),[],1);
    end
    imgIDX = transpose(1:prod(imgSz));
    
    procIDX = find(imdilate(imgout.Mask,strel('disk',2)));
    imgCs = imgCs(procIDX,:);
    imgNs = imgNs(procIDX,:);
    imgIDX = imgIDX(procIDX);
    [~,procIDX] = sort(imgCs(:,1));
    imgCs = imgCs(procIDX,:);
    imgNs = imgNs(procIDX,:);
    imgIDX = imgIDX(procIDX);
    clear procIDX;
    
    projIter = 20;
    projAccu = 0.01;
    minDelta = sqrt(2) * projAccu * pixW;
    
    fSize = 5;
    fRad = pixW * fSize;
    minPtcutoff = 5; %Do not process a RIMLS surface with fewer points
    
    prcLin = 50;
    progT = progressTimer(ceil(imgSz(1)/prcLin),0);
    rStep = progT.reportStep;
    for i=1:prcLin:imgSz(1)
        xcolEnd = min(i+prcLin-1,imgSz(1));
        xLimits = [imgout.xCoords(i,1),imgout.xCoords(xcolEnd,1)];
        firstidx = find(Pts(:,1)>xLimits(1)-fRad,1,'first');
        lastidx = find(Pts(:,1)<xLimits(2)+fRad,1,'last');
        lPts = Pts(firstidx:lastidx,:);
        lNs = Ns(firstidx:lastidx,:);
        
        firstidx = find(imgCs(:,1)>=xLimits(1),1,'first');
        lastidx = find(imgCs(:,1)<=xLimits(2),1,'last');
        lCtrs = imgCs(firstidx:lastidx,:);
        
        KDobj = KDTreeSearcher(lPts(:,1:2)); %Search just XY coordinates
        [lIDX,lDs] = knnsearch(KDobj,lCtrs(:,1:2),'k',100);
        
        newDs = lCtrs(:,3);
        newNs = imgNs(firstidx:lastidx,:);
        parfor j=1:length(newDs) 
            gradient = newNs(j,:);
            npts = find(lDs(j,:)<fRad,1,'last');
            if npts<minPtcutoff, continue; end
            lidx = lIDX(j,1:npts); %#ok<PFBNS>
            npts = lPts(lidx,:); nns = lNs(lidx,:); %#ok<PFBNS>
            projPt = lCtrs(j,:); dZ = 1.0e+09;
            projPt(3) = median(npts(:,3));
            
            its=0;
            while (abs(dZ) > minDelta && its < projIter)
                [flag,dZ,grad] = calculateRIMLS(projPt,npts,nns,fRad);
                if flag
                    projPt(3) = projPt(3) - dZ;
                    gradient = grad;
                    its = its+1;
                else
                    if its<1
                        projPt(3) = minDepth-1;
                    end
                    its = projIter;
                end
            end
            newDs(j) = projPt(3);
            newNs(j,:) = gradient;
        end
        lIDX = imgIDX(firstidx:lastidx);
        imgout.Depths(lIDX) = newDs;
        for k=1:3, imgout.Normals(lIDX+(k-1)*prod(imgSz)) = newNs(:,k); end
        
        if ~mod(i,rStep), progT.update(i/prcLin); end
    end
    progT.done;
    
    %Mask out the pixels contiguous with the border of the image, turning on the rest
    emptyL = bwlabel(~imgout.Mask);
    edgeGroups = nonzeros(unique([emptyL(1,:),emptyL(end,:),emptyL(:,1)',emptyL(:,end)']));
    extMask = ismember(emptyL,edgeGroups);
    extMask = imdilate(extMask,strel('disk',2));
    extMask = ~imdilate(~extMask,strel('disk',2));
    patchMask = extMask | ~imdilate(imgout.Mask,strel('disk',2));
    
    extMask = find(extMask(:)); %Handle by linear indexing
    imgout.Depths(extMask) = minDepth;
    for k=1:3, imgout.Normals(extMask+(k-1)*prod(imgSz)) = 0; end
    
    %Isolate holes remaining inside the model area to be patched
    patchMask = patchMask | imgout.Depths < minDepth;
    patchMask(extMask) = false;
    emptyL = bwlabel(patchMask);
    holeMask = histc(emptyL(:), nonzeros(unique(emptyL)));
    %Holes are gaps inside the model of at least 16 contiguous pixels
    holeMask = ismember(emptyL,find(holeMask>16));
    patchMask = patchMask & ~holeMask;
    
    %Apply new masks to zero out holes and save transparency layers
    imgout.Mask = ~holeMask;
    imgout.Mask(extMask) = false;
    imgout.Transparency(~imgout.Mask) = 0;
    holeMask = find(holeMask(:));
    imgout.Depths(holeMask) = minDepth;
    for k=1:3, imgout.Normals(holeMask+(k-1)*prod(imgSz)) = 0; end
    
    %Flood the remaining areas
    imgout.Depths = regionfill(imgout.Depths,patchMask);
    for k=1:3, imgout.Normals(:,:,k) = regionfill(imgout.Normals(:,:,k),patchMask); end
    imgout.Normals = normalizeCells(imgout.Normals);
    
    fprintf('\b\b complete.\n');
    if ~isempty(paralleljob), delete(paralleljob); end
    clear paralleljob;
end

function [resolved,potential,gradient] = calculateRIMLS(pt,nearPts,nearNs,sRads)
    squ_norm = @(a) sum(a.^2,2);
    %Initialize parameters
    maxIts = 3; refitMin = 1.0e-03;
    sigmaR = 0.5; %Set as a fraction of the local spatial radius
    sigmaN = 0.6; %A low value (0.5?1.5) to sharpen the surface
    invsN2 = 1/sigmaN^2; invsR2 = 1/sigmaR^2;
    
    %Initialize calculated values
    resolved = true;
    potential = 1.0e+09;
    gradient = zeros(1,3);
    
    %Initialize weights, derivatives, and gradient
    npts = size(nearPts,1);
    diffs = repmat(pt,npts,1)-nearPts;
    S = 1./sRads.^2;
    ptwt = 1-sum(diffs.^2,2).*S;
    ptwt(ptwt<0) = 0;
    if sum(ptwt)<10*eps
        resolved = false;
        return;
    end
    initPtwts = ptwt.^4; %npts x 1: initial point weights
    init1Dwts = -8 * S * ptwt.^3; %npts x 1: 1st derivative
    initGwts = diffs.*repmat(init1Dwts,1,3); %npts x 3: gradient weights
    
    %RIMLS iterations
    fs = sum(diffs.*nearNs,2);
    lastGrad = gradient; %Gradient / surface normal
    As = ones(npts,1); %refitting weights
    its=0;
    while (its<1) || (squ_norm((gradient-lastGrad))>refitMin && its <= maxIts)
        lastGrad = gradient;
        if its>0
            As = exp(-invsR2*(fs-potential).^2) .* ...
                exp(-squ_norm(nearNs-repmat(lastGrad,npts,1))*invsN2);
        end
        ptwt = repmat(As.*initPtwts,1,3);
        ptwtsum = sum(ptwt,1);
        if abs(ptwtsum(1)) > 10*eps
            gw = repmat(As,1,3).*initGwts;

            potential = sum(ptwt(:,1).*fs)/ptwtsum(1);
            gradient = (sum(repmat(fs,1,3).* gw,1)-potential*sum(gw,1)...
                + sum(ptwt.*nearNs)) ./ ptwtsum;
            its = its+1;
        else
            if its<2, resolved = false; end
            return;
        end
    end
end