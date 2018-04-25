function [img] = modelInscribedSurface(img,incisionW)
    %Convert the desired incision width into an odd number of pixels
    incisPx = ceil(incisionW/img.PixD);
    incisPx = incisPx+mod(incisPx+1,2);
    imgD = img.Depths/img.PixD;
    
    fprintf('Identifying incisions and other concavities: ');
    %Calculate Radiance Scaling and Ambient Occlusion layers
    img.meanCurvature = radianceScaling(imgD,img.Normals);

    img.ambientOcclusion = calcAmbientOcclusion(imgD,img.Mask,incisPx);
    fprintf('\b\b complete.\n');
    
    %Segmentation to assist finding the original plane of the inscription
    imSeg = slic2_5D(imgD, img.Normals, incisPx);
    imSeg = slicMerge(imSeg, imgD, img.Normals);
    
    fprintf('Modeling original uninscribed plane: ');
    [origPl,origNs] = restoreFace(imSeg,imgD,img.Normals,...
        img.ambientOcclusion,incisPx);
    Nscore = 1-acosd(max(sum(origNs.*img.Normals,3),0))/90;
    
    %Identify the stroke depth from depth, normals, and curvature layers
    %Depth: the strokes will all be near the preserved face
    relD = img.Depths - origPl*img.PixD;
    smDisk = strel('disk',ceil(incisPx/3));
    bigDisk = strel('disk',2*incisPx);
    tmpMask = relD > -incisionW/10; %Inflection slightly below surface
    %This initial mask includes strokes, but also pocks and breaks
    strMask = ~tmpMask & imdilate(tmpMask,bigDisk);
    tmpMask = imclose(abs(relD) < incisionW & img.Mask,bigDisk);
    strMask = strMask & imerode(tmpMask,bigDisk);
    %Normals: mask removes some pocks in large forward-facing areas
    strMask = strMask & imdilate(imgaussfilt(Nscore,3) < 0.8,smDisk);
    %Occlusion: only hollows are retained
    strMask = strMask & (imgaussfilt(img.ambientOcclusion,2) < 0.8);
    %Estimate the inscription depth from the masked points
    inscrD = log(abs(relD(strMask))); %Log-normal distribution
    inscrD = 3*median(inscrD)-2*mean(inscrD)+1.4826*mad(inscrD);
    inscrD = 0.8*min(exp(inscrD),abs(prctile(relD(strMask),5)));
    inscrD = max(min(inscrD,incisionW*2),incisionW/10);
    %Save the results
    img.originalPlane = origPl*img.PixD;
    img.relativeNormals = Nscore;
    img.incisionDepth = inscrD;
    
    %Estimate the central profile line in the x and y directions
    faceProfile = cell(2,1); faceRestored = cell(2,1);
    [faceRestored{1},faceProfile{1}] = mapCenterLine(img,0);
    [faceRestored{2},faceProfile{2}] = mapCenterLine(img,1);
    img.faceProfile = faceProfile;
    img.faceRestored = faceRestored;
    fprintf('\b\b complete.\n');
end

function [faceCL,scanCL] = mapCenterLine(img,rot)
    iMask = rot90(img.Mask,rot);
    imgSz = size(iMask); iPixD = img.PixD;
    iDepth = rot90(img.Depths,rot)/iPixD;
    fDepth = rot90(img.originalPlane,rot)/iPixD;
    relD = iDepth-fDepth;
    width = ceil(1.5*img.incisionDepth/iPixD);
    
    faceMask = imopen(relD > -width & relD < width/2, ...
        strel('disk',ceil(width/3)));
    imMaskY = repmat(1:imgSz(2),imgSz(1),1); fcMaskY = imMaskY;
    fcMaskY(~faceMask) = 0; imMaskY(~iMask) = 0;
    %Set up the sample points for modeling the center line
    xpos = ceil(width/2:width:imgSz(1)-width/2);
    ypos = zeros(1,length(xpos)); lnum = length(xpos);
    for i=1:lnum
        xidx = max(xpos(i)-width,1):min(xpos(i)+width,imgSz(1));
        ypos(i) = mean( [nonzeros(imMaskY(xidx,:)); ...
            repmat(nonzeros(fcMaskY(xidx,:)),2,1)] );
    end
    %Smooth initial y positions estimated from the scan
    yshift = median(ypos); %Shift to preserve median center at edges
    ysm = medfilt1(ypos-yshift,max(ceil(length(ypos)/20),5)) + yshift;
    for i=1:lnum, ypos(i) = mean(ysm(max(i-5,1):min(i+5,lnum))); end
    
    %Calculate depths at the resulting y positions
    faceDctr = ypos; scanDctr = ypos;
    for i=1:lnum
        xidx = max(xpos(i)-width,1):min(xpos(i)+width,imgSz(1));
        yps = round(ypos(i));
        faceDctr(i) = mean(fDepth(xidx,yps)); %Ideal surface depth at y
        
        %Fit a line to the preserved surface near the y position
        [yidx,xidx] = meshgrid(1:imgSz(2),xidx);
        ptsIDX = sub2ind(imgSz,xidx(:),yidx(:));
        ptsIDXscan = ptsIDX(iMask(ptsIDX)); %All scan points in the area
        ptsIDXface = ptsIDX(faceMask(ptsIDX)); %Just points near the face
        [~,ptY] = ind2sub(imgSz,ptsIDXface);
        nscanpt = length(ptsIDXscan); nfacept = length(ptsIDXface);
        %If at least 10% of the scan points are near the original surface,
        %and at least 10% of the nearby points are to either side of the
        %center line, then fit to the nearby points
        if nfacept/(nscanpt+1) > 0.1 && sum(ptY>yps) > 0.1*nfacept ...
            && sum(ptY<yps) > 0.1*nfacept
            wts = max(1 - abs(0.5*relD(ptsIDXface)/width),0);
            wts = wts.*sqrt(1-(yps-ptY)/max(abs(yps-ptY)));
            ptD = iDepth(ptsIDXface);
        else %Otherwise, fit to the scan points
            [~,ptY] = ind2sub(imgSz,ptsIDXscan);
            wts = abs(relD(ptsIDXscan));
            wts = max(1-1.5*wts/max(wts),0).^4;
            ptD = iDepth(ptsIDXscan);
            wtB = sqrt((yps-ptY).^2 + ptD.^2);
            wts = wts.*sqrt(max(1-wtB/max(wtB),0)).^4;
            ptY = ptY(wts>0); ptD = ptD(wts>0); wts = wts(wts>0);
        end
        pvals = [ptY.*wts, wts] \ (ptD.*wts);
        scanDctr(i) = max(pvals(1)*yps + pvals(2), mean(ptD(ptD<0)));
        if faceDctr(i)<scanDctr(i), faceDctr(i) = scanDctr(i); end
    end
    faceCL = iPixD*interp1(xpos,faceDctr,1:imgSz(1),'linear','extrap');
    scanCL = iPixD*interp1(xpos,scanDctr,1:imgSz(1),'linear','extrap');
end

function [faceDs,faceNs] = restoreFace(imSeg,imgD,imgN,imgAO,incisPx)
    imgSz = size(imgD);
    smDsk = strel('disk',2); medDsk = strel('disk',4);
    %Normals: regions facing forward are weighted positively
    tS = min(max(25-acosd(imSeg.planes(:,6)),0)/20,1);
    plScore = floodRegions(imSeg,tS.^2);
    %Scoring by consistency of local normal to trend in segmentation
    tS = imgN; sN = imgN;
    for i=1:3, tS(:,:,i) = imgaussfilt(squeeze(imgN(:,:,i)),1); end
    for i=1:3, sN(:,:,i) = floodRegions(imSeg,imSeg.planes(:,i+3)); end
    tS = min(max(25 - acosd(abs(sum(normalizeCells(tS).*sN,3))),0)/20,1);
    plScore = plScore.*imerode(tS,smDsk);
    %Scoring by size of superpixel (removing small regions)
    tS = cellfun(@length,imSeg.pixelCoords);
    tS(imSeg.regions) = log(tS(imSeg.regions))-log(prctile(tS(imSeg.regions),20));
    plScore = plScore.*floodRegions(imSeg,max(min(tS/max(tS(imSeg.regions)),1),0));
    %Penalizing regions with negative curvature (concave)
    tS = floodRegions(imSeg,imSeg.segmentCurvature);
    plScore = plScore.*max(min(1 + tS/mad(tS(abs(tS)>eps)),1),0);
    %Penalizing highly occluded areas, few of which are at intact regions
    plScore = plScore.*(min(imgAO/0.9,1).^3);
    %Depth: within range established by the higher-scoring regions
    segDs = min(max(restoreDepth(imSeg),min(imgD(:))),max(imgD(:)));
    estDepths = segDs(plScore > 0.3);
    estDepths = estDepths(estDepths > -prctile(estDepths(estDepths>0),95));
    Dcutoff = 3*1.4826*mad(estDepths);
    plScore = plScore.*min(max(1-abs(imgD)/Dcutoff,0),1);
    
    step = max(min(2*incisPx,min(imgSz)),10);
    faceDs = zeros(imgSz); imTrsp = faceDs;
    progT = progressTimer(ceil(imgSz(1)/step),0);
    rStep = progT.reportStep;
    
    plSumPP = sum(plScore(:))/prod(imgSz);
    for i=1:ceil(imgSz(1)/step)
        xr = min(max(floor([(i-3)*step,(i+2)*step]),1),imgSz(1));
        for j=1:ceil(imgSz(2)/step)
            yr = min(max(floor([(j-3)*step,(j+2)*step]),1),imgSz(2));
            ijsz = [xr(2)-xr(1)+1,yr(2)-yr(1)+1];
            ijnp = prod(ijsz);
            [gy,gx] = meshgrid(yr(1):yr(2),xr(1):xr(2));
            xyIDX = sub2ind(imgSz,gx,gy);
            
            ijW = plScore(xyIDX);
            faceIDX = find(ijW > 0.05);
            ijW = ijW(faceIDX);
            ijD = zeros(ijsz);
            plNormal = [0,0,0];
            if length(faceIDX) > 0.01*ijnp
                ijD = imgD(xyIDX);
                ijPt = []; [ijPt(:,1),ijPt(:,2)] = ind2sub(ijsz,faceIDX);
                ijPt(:,3) = ijD(faceIDX);
                %Fit a plane to the points
                plCtr = sum(ijPt.*repmat(ijW,1,3),1)/sum(ijW);
                ijPt = (ijPt-plCtr).*ijW;                
                [plNormal,~] = eig(ijPt'*ijPt); plNormal = plNormal(:,1);
                if plNormal(3)<0, plNormal = -plNormal; end
                if plNormal(3)>eps %Calculate depths from the fitted plane
                    [ijgy,ijgx] = meshgrid(1:ijsz(2),1:ijsz(1));
                    ijD = [ijgx(:),ijgy(:),zeros(ijnp,1)];
                    ijD = ijD - repmat(plCtr,ijnp,1);
                    ijD = reshape(ijD*-plNormal/plNormal(3),ijsz);
                    ijD = max(min(ijD,Dcutoff),-Dcutoff);
                else
                    ijD = zeros(ijsz); plNormal = [0,0,0];
                end
            end
            if plNormal(3) < eps, newTransp = zeros(size(ijD));
            else
                %Transparency begins with the distance from center
                newTransp = sqrt((gx-(i-0.5)*step).^2+(gy-(j-0.5)*step).^2);
                newTransp = max(1-newTransp/(2.5*step),0);
                %Lower weight according to direction of the fitted plane
                plNormScores = 1-min(acosd(abs(plNormal(3)))/45,1);
                newTransp = newTransp * plNormScores;
                %Increase weight according to number of points in region
                newTransp = newTransp * (1+length(faceIDX)/ijnp);
                %Reweight by cumulative score of points
                plSumRatios = sum(ijW)/(ijnp*plSumPP);
                newTransp = newTransp * sqrt(plSumRatios);
                if sum(newTransp(:))
                    %Raise the model to keep it above all scored points
                    ijDorig = min(imgD(xyIDX),Dcutoff/2);
                    lowPts = sum(ijD(:) < ijDorig(:));
                    ijD = max(ijD,ijDorig);
                    ijD = imgaussfilt(imdilate(ijD,medDsk),3);
                    %Reduce the transparency by the percentage of low points
                    newTransp = newTransp * (1-lowPts/ijnp);
                end
            end
            
            faceDs(xyIDX) = faceDs(xyIDX) + ijD.*newTransp;
            imTrsp(xyIDX) = imTrsp(xyIDX) + newTransp;
        end
        if ~mod(i,rStep), progT.update(i); end
    end
    progT.done;
    fillMask = imerode(imTrsp < 2,smDsk) | imdilate(imTrsp < eps,smDsk);
    fillMask = fillMask & imdilate(imTrsp>0,strel('disk',step*3));
    imTrsp(imTrsp<0.01) = 0.01;
    faceDs = regionfill(faceDs./imTrsp,fillMask);
    faceDs = imgaussfilt(faceDs,1.5*incisPx);
    %Normals derived from the reconstructed plane
    plDpad = padarray(faceDs,[1 1],'replicate','both');
    faceNs = (plDpad(3:end,2:end-1)-plDpad(1:end-2,2:end-1))/2;
    faceNs(:,:,2) = (plDpad(2:end-1,3:end)-plDpad(2:end-1,1:end-2))/2;
    faceNs(:,:,3) = ones(imgSz);
    faceNs = normalizeCells(faceNs);
end

function [imout] = restoreDepth(seg)
    %Create an image with the depth channel from the fitted planes
    imsz = size(seg.labels);
    imout= zeros(imsz);
    plCtr = seg.planes(seg.regions,1:3);
    plNrm = seg.planes(seg.regions,4:6);
    for i=1:length(seg.regions)
        segPix = seg.pixelCoords{seg.regions(i)};
        ptXYZ = []; [ptXYZ(:,1),ptXYZ(:,2)] = ind2sub(imsz,segPix);
        ptXYZ(:,2) = -ptXYZ(:,2); ptXYZ(:,3) = zeros(length(segPix),1);
        ptXYZ = ptXYZ - repmat(plCtr(i,:),length(segPix),1);
        imout(segPix) = ptXYZ * transpose(-plNrm(i,:)) / plNrm(i,3);
    end
end

function [imout] = floodRegions(seg,rdata)
    %Create an image with the provided values flooded in the region
    imout = zeros(size(seg.labels));
    for i=1:length(seg.regions)
        ir = seg.regions(i);
        segPix = seg.pixelCoords{ir};
        if ~isempty(segPix), imout(segPix) = rdata(ir); end
    end
end

function [imAO] = calcAmbientOcclusion(imD,imMask,incisPx)
    %The depth layer is provided at the XY pixel scale
    delchars = fprintf('Ambient occlusion: ');
    %Radius assessed around each point, constrained within 5?50 pixels
    radWindow = min(max(ceil(incisPx*2.5),5),50);
    nDivs = 32; %Number of radial segments checked around each point
    if radWindow <= 16
        if radWindow <= 8, nDivs = 8; else, nDivs = 16; end
    end
    tilesz = repmat(2*radWindow+1,1,2); %Width gathered around each point
    
    angRange = 360*((0:nDivs)-0.5)/nDivs;
    tileXs = repmat(transpose(-radWindow:radWindow),1,tilesz(1));
    tileYs = repmat(-radWindow:radWindow,tilesz(1),1);
    tileAngle = atan2d(-tileYs,tileXs);
    tileAngle(tileAngle<0) = tileAngle(tileAngle<0)+360;
    %Precomputed distances around each point
    tileDst = sqrt(tileXs.^2+tileYs.^2);
    
    %Precalculated indices of points in each ray around the center point
    dvPts = cell(nDivs,1);
    ipts = find(tileAngle<=angRange(2) | tileAngle>angRange(1)+360);
    [~,sortIDX] = sort(tileDst(ipts));
    dvPts{1} = ipts(sortIDX(2:end)); %Remove the center point
    for i=2:nDivs
        ipts = find(tileAngle<=angRange(i+1) & tileAngle>angRange(i));
        [~,sortIDX] = sort(tileDst(ipts));
        dvPts{i} = ipts(sortIDX);
    end
    minL = min(cellfun(@length,dvPts));
    dvPts = cellfun(@(x) x(1:minL,:),dvPts,'UniformOutput',false);
    %Keeping an even number of points for each ray, within the same radius
    dvPts = transpose(reshape(cell2mat(dvPts),minL,[]));
    tileDst = 1./tileDst(dvPts); %XY distances within each ray
    
    imgSz = size(imD);
    imD = imgaussfilt(imD,0.5);
    imAO = zeros(imgSz);
    imDpad = padarray(imD,[radWindow,radWindow],min(imD(:)),'both');
    imMpad = double(padarray(imMask,[radWindow,radWindow],'both'));
    
    progT = progressTimer(imgSz(1),0);
    rStep = progT.reportStep;
    for i=1:imgSz(1)
        dvMline = imMpad(i:i+2*radWindow,:);
        dvDline = imDpad(i:i+2*radWindow,:);
        for j=1:imgSz(2)
            dvMsk = dvMline(:,j:j+2*radWindow);
            dvDs = dvDline(:,j:j+2*radWindow) - imD(i,j);
            dvTh = atand(max(tileDst.*dvDs(dvPts).*dvMsk(dvPts),[],2));
            imAO(i,j) = 180-median(dvTh)-mean(dvTh);
        end
        if ~mod(i,rStep), progT.update(i); end
    end
    progT.done;
    
    imAO = imgaussfilt(min(imAO/180,1),0.25);
    imAO(~imMask) = 0;
    imAO = imAO.^2;
    fprintf(repmat(char(8),delchars,1));
end

function [meanCurv] = radianceScaling(imgD,imgN)
    %The depth layer is provided at the XY pixel scale
    %Calculate the curvature by the method proposed in Vergne et al 2010
    epsRSI = 0.01; foreshortening = 0.2;
    gS = imgN(:,:,3);
    gS(gS<epsRSI) = epsRSI;
    gS = gS.^-foreshortening;
    gPx = -imgN(:,:,1).*gS;
    gPy = -imgN(:,:,2).*gS;
    [~, gPxx] = gradient(gPx);
    [gPyy, ~] = gradient(gPy);
    meanCurv = (gPxx-gPyy)/2;
    
    %Sobel filter of the depth isolates silhouettes
    wsPx = imfilter(imgD,fspecial('sobel'),'conv');
    wsPy = imfilter(imgD,fspecial('sobel')','conv');
    wsP = sqrt(wsPx.^2+wsPy.^2);
    %Use Otsu's method to distinguish clusters of points on silhouettes
    grayThresh = multithresh(nonzeros(wsP),3);
    wsP = wsP/grayThresh(1);
    
    wsN1x = imfilter(imgN(:,:,1),fspecial('sobel'),'conv');
    wsN1y = imfilter(imgN(:,:,1),fspecial('sobel')','conv');
    wsN1 = abs(wsN1x)+abs(wsN1y);
    wsN2x = imfilter(imgN(:,:,2),fspecial('sobel'),'conv');
    wsN2y = imfilter(imgN(:,:,2),fspecial('sobel')','conv');
    wsN2 = abs(wsN2x)+abs(wsN2y);
    wsN = wsN1.^2+wsN2.^2;
    grayThresh = multithresh(nonzeros(wsN),3);
    wsN = wsN/grayThresh(2);
    wts = 1-max(wsP,wsN);
    wts(wts<0) = 0;
    
    smoothCurv = anisotropicDiffusion(meanCurv,0);
    smoothCurv = anisotropicDiffusion(smoothCurv,wts);
    
    meanCurv = (meanCurv+smoothCurv)/2;
end

function diffG = anisotropicDiffusion(G,wts)
    diffG = G;
    kappa = 5;
    passBuffer = zeros([size(G),4]);
    if sum(abs(wts(:)))>0, Cp = wts/4; else, Cp = 0.25; end
    
    filt = zeros(3,3,4);
    filt(:,:,1) = [0 1 0; 0 -1 0; 0 0 0];
    filt(:,:,2) = [0 0 0; 0 -1 0; 0 1 0];
    filt(:,:,3) = [0 0 0; 0 -1 1; 0 0 0];
    filt(:,:,4) = [0 0 0; 1 -1 0; 0 0 0];
    for i=1:3
        for j=1:4, passBuffer(:,:,j) = imfilter(diffG,filt(:,:,j),'conv'); end
        passConductance = 1./(1+(passBuffer/kappa).^2);
        diffG = diffG + Cp .* sum(passBuffer.*passConductance,3);
    end
end