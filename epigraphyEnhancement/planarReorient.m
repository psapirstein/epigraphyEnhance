function [transformedcloud, pixelWidth] =  planarReorient(ptcloud,samples)
    median_normal = @(x) median(x)'/sqrt(sum(median(x).^2));
    
    k = 200; TargetPts = 250000; ptdeci = ptcloud; delChar = 0;
    if TargetPts/ptdeci.Count < 1
        delChar = fprintf('Decimating cloud...');
        sample = floor(ptdeci.Count/TargetPts);
        if sample<6, sample = 6; end
        ptdeci = pcdownsample(ptdeci,'nonuniformGridSample',sample);
    end
    Pt = ptdeci.Location;
    Ns = ptdeci.Normal;
    NP = ptdeci.Count;
    clear PCdeci;
    
    fprintf(repmat(char(8),delChar));
    delChar = fprintf('Isolating flattest areas in the point cloud...');
    kdOBJ = KDTreeSearcher(Pt);
    KNNidx = knnsearch(kdOBJ,Pt,'K',k);
    DPlocal = zeros(size(Pt,1),1);
    for i=1:NP
        nearbyVs = Pt(KNNidx(i,:),:);
        [~,~,V] = svd(nearbyVs-repmat(mean(nearbyVs,1),k,1),0);
        DPlocal(i) = Ns(i,:)*V(:,3);
    end
    DPlocal = abs(DPlocal);
    
    %Now isolate the dominant points on the face to be analyzed
    flatIDX = DPlocal>prctile(DPlocal,95); %Only monitor top 5% flattest points
    MN = median_normal(Ns(flatIDX,:));
    DPmedian = Ns*MN;
    Q = DPmedian.*DPlocal; Q(Q<0) = 0;
    
    %Isolate the points whose normals are close to the dominant face
    hiQpt = Pt(Q>prctile(Q,95),:);
    MN = median_normal(Ns(Q>prctile(Q,95),:));
    hiQdist = (hiQpt - repmat(mean(hiQpt),size(hiQpt,1),1))*MN;
    %Select the more forward points
    testIDX = hiQdist>-std(hiQdist);
    testSpan = std(hiQdist)/25;
    hiQpt = hiQpt(testIDX,:);
    
    fprintf(repmat(char(8),delChar));
    fprintf('Fitting plane to primary inscribed surface with RANSAC...');
    trials = 50000;
    [PN,newCtr] = ransacPlane(hiQpt,trials,testSpan);
    if PN*MN < 0, PN = -PN; end
    
    finalTrans = transpose(makehgtform('axisrotate',cross(PN,[0,0,1]),acos(dot(PN,[0,0,1]))));
    finalTrans = transpose(makehgtform('translate',-newCtr))*finalTrans;
    Pt = transformPointsForward(affine3d(finalTrans),Pt);
    [~,angZ] = minBoundingBox(double(Pt(:,1:2)'));
    zRot = transpose(makehgtform('zrotate',pi/2-angZ));
    Pt = transformPointsForward(affine3d(zRot),Pt);
    if max(Pt(:,1))-min(Pt(:,1)) < max(Pt(:,2))-min(Pt(:,2))
        zRot = transpose(makehgtform('zrotate',pi-angZ));
    end
    finalTrans = finalTrans*zRot;
    
    transformedcloud = pctransform(ptcloud,affine3d(finalTrans));
    %Recenter the resulting cloud
    newXYctr = [sum(transformedcloud.XLimits), sum(transformedcloud.YLimits)];
    newXYctr = [-newXYctr/2,0];
    transformedcloud = pctransform(transformedcloud, ...
        affine3d(transpose(makehgtform('translate',newXYctr))));
    
    %Calculate the raster pixel spacing from the projected scan data
    spaces = samplePointSpacing(double(transformedcloud.Location));
    pixelWidth = getPPD(median(spaces)*sqrt(samples));
end

function [PN, CTR] = ransacPlane(datapt, ntrials, tol)
% Simplified from code (c) 2003-2013 Peter Kovesi
% Centre for Exploration Targeting, The University of Western Australia
% peter.kovesi at uwa edu au
% http://www.csse.uwa.edu.au/~pk
    isdegenerate = @(X) norm(cross(X(2,:)-X(1,:), X(3,:)-X(1,:))) < eps;
    
    maxDataTrials = 100;
    npts = size(datapt,1);
    
    bestscore =  0;
    besttrial = 0; % Sentinel value allowing detection of solution failure.
    
    for trial=1:ntrials
        % Select at random s datapoints to form a trial model, M.
        % In selecting these points we have to check that they are not in
        % a degenerate configuration.
        degenerate = 1;
        count = 1;
        while degenerate
            ind = randsample(npts, 3);
            % Test that these points are not a degenerate configuration.
            degenerate = isdegenerate(datapt(ind,:));
            if ~degenerate
                Ptri = datapt(ind,:);
                if isempty(Ptri), degenerate = 1; end
            end
            
            count = count + 1;
            if count > maxDataTrials
                % Safeguard against being stuck in this loop forever
                warning(' Unable to select a nondegenerate data set\n');
                break
            end
        end
        
        % Once we are out here we should have some kind of model...
        % Evaluate distances between points and model returning the indices
        % of elements in x that are inliers.
        planeN = cross(Ptri(2,:)-Ptri(1,:), Ptri(3,:)-Ptri(1,:));
        planeN = (planeN/sqrt(sum(planeN.^2)))';
        inliers = abs((datapt-repmat(Ptri(1,:),size(datapt,1),1))*planeN) < tol;
        if sum(inliers) > bestscore    % Largest set of inliers so far...
            bestscore = sum(inliers);  % Record data for this model
            bestinliers = inliers;
            besttrial = trial;
        end
    end
    
    if besttrial>0
        CTR = mean(datapt(bestinliers,:));
        [~,~,PN] = svd(datapt(bestinliers,:)-repmat(CTR,sum(bestinliers),1));
        PN = PN(:,3)';
        fprintf('\b\b\b: fitted %d/%d points in trial %d\n',bestscore,npts,besttrial);
    else
        PN = [];
        CTR = [];
        warning(' Unable to find a useful solution\n');
    end
end

function pW = getPPD(coarsePPD)
% Fix the pixel scale to a common ratio like 10:2, 4:1, 10:4, 2:1, etc.
    preferredscales = [1,2,2.5,4,5,8,10];
    basis = floor(log10(coarsePPD));
    psidx = find(preferredscales * 10^basis > coarsePPD,1,'first');
    pW = preferredscales(psidx)*10^basis;
end

function spc = samplePointSpacing(ptsIn)
    boundsX = 5:10:95; boundsY = boundsX;
    minX = min(ptsIn(:,1)); minY = min(ptsIn(:,2));
    boundsX = minX + (max(ptsIn(:,1))-minX)*[boundsX'-0.2,boundsX']/100;
    boundsY = minY + (max(ptsIn(:,2))-minY)*[boundsY'-0.2,boundsY']/100;
    [~,xidx] = sort(ptsIn(:,1));
    ptsIn = ptsIn(xidx,:);
    spc = zeros(100,1);
    for i=1:10
        firstidx = find(ptsIn(:,1) > boundsX(i,1),1,'first');
        nextidx = find(ptsIn(:,1) <= boundsX(i,2),1,'last');
        ptsLocal = ptsIn(firstidx:nextidx,:);
        [~,yidx] = sort(ptsLocal(:,2));
        ptsLocal = ptsLocal(yidx,:);
        for j=1:10
            firstidx = find(ptsLocal(:,2) > boundsY(j,1),1,'first');
            nextidx = find(ptsLocal(:,2) <= boundsY(j,2),1,'last');
            if nextidx-firstidx > 100
                pts = ptsLocal(firstidx:nextidx,:);
                [~,neighbors] = knnsearch(KDTreeSearcher(pts),pts,'k',2);
                spc((i-1)*10+j) = median(neighbors(:,2));
            end
        end
    end
    spc = nonzeros(spc);
end

function [bb,Rout] = minBoundingBox(X)
    % Adapted from https://www.mathworks.com/matlabcentral/fileexchange/...
    %       31126-2d-minimal-bounding-box/content/minBoundingBox.m
    % Created by Julien Diener: https://www.mathworks.com/matlabcentral/...
    %       profile/authors/1978578-julien-diener
    % Copyright (c) 2011, julien diener
    % All rights reserved.
    %
    % Compute the minimum bounding box of a set of 2D points
    %   Use:   boundingBox = minBoundingBox(point_matrix)
    % Input:  2xn matrix containing the [x,y] coordinates of n points
    %         *** there must be at least 3 points which are not collinear
    % output: 2x4 matrix containing the coordinates of the bounding box corners
    
    % compute the convex hull (CH is a 2*k matrix subset of X)
    k = convhull(X(1,:),X(2,:));
    CH = X(:,k);
    % compute the angle to test, which are the angle of the CH edges as:
    %   "one side of the bounding box contains an edge of the convex hull"
    E = diff(CH,1,2);           % CH edges
    T = atan2(E(2,:),E(1,:));   % angle of CH edges (used for rotation)
    T = unique(mod(T,pi/2));    % reduced to the unique set of first quadrant angles
    % create rotation matrix which contains the 2x2 rotation matrices for *all* angles in T
    % R is a 2n*2 matrix
    R = cos( reshape(repmat(T,2,2),2*length(T),2) ... % duplicate angles in T
           + repmat([0 -pi ; pi 0]/2,length(T),1));   % shift angle to convert sine in cosine
    % rotate CH by all angles
    RCH = R*CH;
    % compute border size [w1;h1;w2;h2;....;wn;hn] & area of bounding box for all possible edges
    bsize = max(RCH,[],2) - min(RCH,[],2);
    area  = prod(reshape(bsize,2,length(bsize)/2));
    % find minimal area, thus the index of the angle in T 
    [~,i] = min(area);
    % compute the bound (min and max) on the rotated frame
    Rf    = R(2*i+[-1 0],:);   % rotated frame
    bound = Rf * CH;           % project CH on the rotated frame
    bmin  = min(bound,[],2);
    bmax  = max(bound,[],2);
    % compute the corner of the bounding box
    Rf = Rf';
    bb(:,4) = bmax(1)*Rf(:,1) + bmin(2)*Rf(:,2);
    bb(:,1) = bmin(1)*Rf(:,1) + bmin(2)*Rf(:,2);
    bb(:,2) = bmin(1)*Rf(:,1) + bmax(2)*Rf(:,2);
    bb(:,3) = bmax(1)*Rf(:,1) + bmax(2)*Rf(:,2);
    % Added export of the rotation to minimize the bounding box
    Rout = T(i);
end