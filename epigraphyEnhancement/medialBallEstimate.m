function MATobj = medialBallEstimate(ptcloud)
    delChars = fprintf('Downsampling scan...');
    
    ptTarget = 2000000; %Reduced MAT model to avoid long computation
    ptTarget = floor(ptcloud.Count/ptTarget);
    if ptTarget<8, ptTarget = 8; end
    decicloud = pcdownsample(ptcloud,'nonuniformGridSample',ptTarget);
    pts = double(decicloud.Location);
    Ns = normalizeCells(double(decicloud.Normal));
    clear decicloud ptcloud;
    
    delChars = delChars + fprintf('\b\b Removing narrow, sharp features');
    for its=1:5
        %Average normals from the nearest 20 points
        [nearIDX,nearD] = knnsearch(KDTreeSearcher(pts),pts,'k',21);
        nearIDX = nearIDX(:,2:end); nearD = nearD(:,2:end);
        meanNs = zeros(size(Ns));
        for dim=1:size(nearIDX,2), meanNs = meanNs + Ns(nearIDX(:,dim),:); end
        meanNs = normalizeCells(meanNs);
        normDiff = sum(Ns.*meanNs,2);
        %Delete points which do not match their neighborhood well
        keepIDX = normDiff > 0.85;
        %Delete points more than twice the distance from their nearest
        %neighbor as that neighbor's nearest point
        keepIDX = keepIDX & nearD(:,1) < 2*nearD(nearIDX(:,1),1);
        pts = pts(keepIDX,:); Ns = Ns(keepIDX,:);
    end
    np = size(pts,1);
    [~,ptIDX] = sort(pts(:,1)); %Sort by x values for later calculation
    pts = pts(ptIDX,:); Ns = Ns(ptIDX,:);
    clear delIDX nearIDX normDiff meanNs;
    
    ptTarget = 50000; %Low resolution model for general area
    ptslo = pts;
    if floor(np/ptTarget)>2
        ptslo = pcdownsample(pointCloud(ptslo),'nonuniformGridSample',floor(np/ptTarget));
        ptslo = double(ptslo.Location);
        [~,ptIDX] = sort(ptslo(:,1));
        ptslo = ptslo(ptIDX,:);
    end
    xyzLim = double(max(pts,[],1)-min(pts,[],1));
    
    ptTarget = 30000; %for the region of interest
    Tpreserve = 10; %preservation angle in degrees; see Peters & Ledoux 2016, p. 127
    minDelta = max(xyzLim)/200000; %when the readius changes less than 1 part
        %in 200k relative to the maximum dimension of the model, stop
    
    cellLen = nthroot(prod(xyzLim)*ptTarget/np,3);
    initR = 5*cellLen;
    roiLBs = double(min(pts,[],1)) - (cellLen*ceil(xyzLim/cellLen) - xyzLim)/2;
    roiUBs = roiLBs + cellLen*(ceil(xyzLim/cellLen)-1);
    [Xlo,Ylo,Zlo] = meshgrid(roiLBs(1):cellLen:roiUBs(1), ...
        roiLBs(2):cellLen:roiUBs(2), roiLBs(3):cellLen:roiUBs(3));
    
    rois = [Xlo(:)';Ylo(:)';Zlo(:)'];
    innerRs = zeros(np,1); outerRs = zeros(np,1);
    ROIobj = struct;
    
    fprintf(repmat(char(8),delChars));
    delChars = fprintf('Calculating Medial Axis Transform: ');
    progT = progressTimer(size(rois,2),0);
    rStep = progT.reportStep;
    for i=1:size(rois,2)
        thisROI = transpose([rois(:,i)-0.2*cellLen,rois(:,i)+1.2*cellLen]);
        knnpts = collectROIpts(thisROI,pts,true);
        outerpts = collectROIpts(thisROI,ptslo,false);
        ROIobj.knnPt = [knnpts;outerpts];
        ROIobj.kdOBJ = KDTreeSearcher(ROIobj.knnPt,'BucketSize',100);
        
        thisROI = transpose([rois(:,i),rois(:,i)+cellLen]);
        [ROIobj.roiPts, thisIDX] = collectROIpts(thisROI,pts,true);
        ROIobj.roiNs = Ns(thisIDX,:);
        
        innerRs(thisIDX) = calculateMAD(ROIobj,initR,Tpreserve,minDelta,true);
        outerRs(thisIDX) = calculateMAD(ROIobj,initR,Tpreserve,minDelta,false);
        if ~mod(i,rStep), progT.update(i); end
    end
    progT.done;
    inCtrs = pts - repmat(innerRs,1,3).*Ns;
    outCtrs = pts + repmat(outerRs,1,3).*Ns;
    
    %delete any radius that was not adjusted
    inKeepIDX = innerRs < initR;
    outKeepIDX = outerRs < initR;
    MATobj.inCs = inCtrs(inKeepIDX,:);
    MATobj.inRs = innerRs(inKeepIDX);
    MATobj.outCs = outCtrs(outKeepIDX,:);
    MATobj.outRs = outerRs(outKeepIDX,:);
    fprintf(repmat(char(8),1,delChars));
    fprintf('Medial Axes calculated.\n');
end

function [ptout,ptsIDX] = collectROIpts(ROI,pt,insideROI)
    ptsIDX = pt <= ROI(2,:);
    ptsIDX = ptsIDX & pt > ROI(1,:);
    ptsIDX = all(ptsIDX');
    if ~insideROI, ptsIDX = ~ptsIDX; end
    ptout = pt(ptsIDX,:);
end

function [ballRs] = calculateMAD(obj,initR,Tpres,minD,inside)
    Pts = obj.roiPts; Ns = obj.roiNs;
    if ~inside, Ns = -Ns; end
    ballRs = repmat(initR,size(Pts,1),1);
    ballCs = Pts - Ns.*repmat(ballRs,1,3);
    roiTestidx = 1:size(Pts,1);
    iter = 1; maxIts = 30;
    while iter<maxIts
        roiQidx = knnsearch(obj.kdOBJ,ballCs(roiTestidx,:));
        roiQs = obj.knnPt(roiQidx,:);
        ds = sqrt(sum((Pts-roiQs).^2,2));
        completedMADs = ds < 2*eps;
        ds(completedMADs) = 1;
        cosTheta = sum(Ns.*(Pts-roiQs),2)./ds;
        newRs = ds./(2*cosTheta);

        completedMADs = completedMADs | newRs < 0; %Also stop testing if R appears to be zero
        completedMADs = completedMADs | newRs > initR; %and if it is expanding

        newCs = Pts - Ns.*repmat(newRs,1,3);

        %Preservation angle filtering
        completedMADs = completedMADs | cosAngle(Pts-newCs,roiQs-newCs) < Tpres;
        completedMADs = completedMADs | ballRs(roiTestidx)-newRs < minD;

        %Delete the invalid or completed points
        Ns(completedMADs,:)=[];
        Pts(completedMADs,:)=[];
        roiTestidx(completedMADs)=[];

        ballRs(roiTestidx) = newRs(~completedMADs);
        ballCs(roiTestidx,:) = newCs(~completedMADs,:);

        iter = iter+1;
        if isempty(roiTestidx), iter=maxIts; end
    end
end

function [angAB] = cosAngle(A,B)
    dA = sqrt(sum(A.^2,2)); dA(dA<2*eps) = 1;
    dB = sqrt(sum(A.^2,2)); dB(dB<2*eps) = 1;
    angAB = acosd(sum(A.*B,2) ./ (dA.*dB));
end