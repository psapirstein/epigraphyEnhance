function [img2_5D,trimmedCloud] = rasterize2_5D(ptCloud,MATobj,pixelW,incW)
    fprintf('Rasterizing at %.1f pixels per mm.',round(4/pixelW/1000)/4);
    Pt = double(ptCloud.Location);
    Cs = double(ptCloud.Color);
    Ns = double(ptCloud.Normal);
    clear ptCloud;
    
    %Remove points facing away from the current view
    forwardIDX = Ns(:,3)>0;
    Pt = Pt(forwardIDX,:); Cs = Cs(forwardIDX,:); Ns = Ns(forwardIDX,:);
    %Sort scan points by increasing depth for blending
    [~,PTidx] = sort(Pt(:,3));
    Pt = Pt(PTidx,:); Ns = Ns(PTidx,:); Cs = Cs(PTidx,:);
    %Recenter the point cloud on x-y plane
    maxPt = max(Pt); minPt = min(Pt);
    shiftPt = minPt + maxPt; shiftPt(3) = 0;
    shiftPt = shiftPt/2;
    Pt = Pt - shiftPt;
    
    nPix = ceil(max(abs(Pt))/pixelW)*2+1;
    ctrCoord = [0,0,0];
    
    xCs = ctrCoord(1) + ((1:nPix(1))-1-floor(nPix(1)/2))*pixelW;
    yCs = ctrCoord(2) - ((1:nPix(2))-1-floor(nPix(2)/2))*pixelW;
    xCs = repmat(xCs',1,nPix(2));
    yCs = repmat(yCs,nPix(1),1);
    
    delChars = fprintf(' Mapping interior: ');
    %Sort inside points to begin with those farthest below the inscription
    [~,cidx] = sort(MATobj.inCs(:,3));
    inCs = MATobj.inCs(cidx,:) - shiftPt; inRs = MATobj.inRs(cidx);
    outRs = MATobj.outRs; outCs = MATobj.outCs - shiftPt;
    clear MATobj;
    
    %Delete all MAT points centered significantly outside the image space
    buf = (maxPt-minPt)/20; buf = (buf+mean(buf))/2;
    MATidx = inCs(:,1) < max(Pt(:,1))+buf(1) & inCs(:,1) > min(Pt(:,1))-buf(1);
    MATidx = MATidx & inCs(:,2) < max(Pt(:,2))+buf(2) & inCs(:,2) > min(Pt(:,2))-buf(2);
    MATidx = MATidx & inCs(:,3) < max(Pt(:,3))+buf(3);
    inCs = inCs(MATidx,:); inRs = inRs(MATidx);
    %Delete outside points with large radii, which may intrude below the
    %original surface due to downsampling
    MATidx = outRs < mean(outRs)+4*std(outRs);
    outCs = outCs(MATidx,:); outRs = outRs(MATidx);
    %Delete points that are overlapped by the outside medial balls
    [outIDX,outDist] = knnsearch(KDTreeSearcher(outCs),inCs,'k',50);
    MATidx = zeros(size(outDist));
    for i=1:size(outDist,2)
        MATidx(:,i) = inRs + outRs(outIDX(:,i)) - outDist(:,i);
    end
    MATidx = max(MATidx,[],2);
    MATidx = MATidx < prctile(MATidx, 90); %Discard highest 10% overlaps
    inCs = inCs(MATidx,:); inRs = inRs(MATidx);
    %Delete balls with with very large radii, but not more than top 1/1000
    maxRad = max(inRs)/4;
    maxRad = max(maxRad,prctile(inRs,99.9));
    MATidx = inRs < maxRad;
    inCs = inCs(MATidx,:); inRs = inRs(MATidx);
    
    %Map circle centers to the image space
    MATxy = round(getPointIDX(inCs,nPix,ctrCoord,pixelW));
    ballRad = round(inRs/pixelW);
    
    mindepthMap = minPt(3)*ones(nPix(1:2));
    progT = progressTimer(length(ballRad),0);
    rStep = progT.reportStep;
    for i=1:length(ballRad)
        pixRad = ballRad(i);
        xyp = MATxy(i,:);
        if pixRad<1 %Evaluate just at the current raster point
            if all(xyp(1:2) > 0) && all(xyp(1:2) < nPix(1:2))
                localRad = sum(([xCs(xyp(1),xyp(2)),yCs(xyp(1), ...
                    xyp(2))]-inCs(i,1:2)).^2,2);
                localH = inRs(i)^2 - localRad;
                if localH>0
                    localD = inCs(i,3) + sqrt(localH);
                    if localD > mindepthMap(xyp(1),xyp(2))
                        mindepthMap(xyp(1),xyp(2)) = localD;
                    end
                elseif inCs(i,3) > mindepthMap(xyp(1),xyp(2))
                    mindepthMap(xyp(1),xyp(2)) = inCs(i,3);
                end
            end
        else %Evaluate for the grid of cells within the ball radius
            xs = xyp(1)-pixRad:xyp(1)+pixRad;
            ys = xyp(2)-pixRad:xyp(2)+pixRad;
            xs(xs<1) = []; xs(xs>nPix(1)) = [];
            ys(ys<1) = []; ys(ys>nPix(2)) = [];
            if ~isempty(xs) && ~isempty(ys)
                xGrid = repmat(xs,length(ys),1);
                yGrid = repmat(ys',1,length(xs));
                pixIDX = sub2ind(nPix,xGrid(:),yGrid(:));

                oldDepth = mindepthMap(pixIDX);
                evalXY = [xCs(pixIDX),yCs(pixIDX)];
                localRad = sum((evalXY-inCs(i,1:2)).^2,2);
                localH = inRs(i)^2 - localRad;
                localH(localH<0) = 0;
                localD = inCs(i,3) + sqrt(localH);
                %Ensure that all cells outside the ball radius are not considered
                localD(localH < 10*eps) = minPt(3);

                higherIDX = find((localD-oldDepth)>0);
                if ~isempty(higherIDX)
                    mindepthMap(pixIDX(higherIDX)) = localD(higherIDX);
                end
            end
        end
        if ~mod(i,rStep), progT.update(i); end
    end
    progT.done;
    fprintf(repmat(char(8),1,delChars));
    delChars = fprintf(' Blending points: ');
    
    ptIDX = getPointIDX(Pt,nPix,ctrCoord,pixelW);
    [pixIDX,pixAreas,shiftIDX] = mapPixels(ptIDX(:,1:2),nPix(1:2));
    
    %Select just the points above or slightly below the MAT image
    mindepthPts = zeros(size(pixIDX));
    for i=1:4, mindepthPts = mindepthPts + ...
            pixAreas(:,i).*mindepthMap(pixIDX+shiftIDX(i)); end
    ptIDX = find(Pt(:,3) > mindepthPts-(pixelW+incW)/2);
    Pt = Pt(ptIDX,:); Ns = Ns(ptIDX,:); Cs = Cs(ptIDX,:);
    pixIDX = pixIDX(ptIDX); pixAreas = pixAreas(ptIDX,:);
    
    img2_5D = struct;
    img2_5D.PixD = pixelW;
    img2_5D.xCoords = xCs;
    img2_5D.yCoords = yCs;
    
    %Default image depth is the lowest depth
    imgDepths = minPt(3)*ones(nPix(1:2));
    imgTransparency = zeros(nPix(1:2));
    imgNormals = zeros([nPix(1) nPix(2) 3]); imgColors = zeros(size(imgNormals));
    
    progT = progressTimer(length(pixIDX),0);
    rStep = progT.reportStep;
    for i=1:length(pixIDX)
        localIDX = pixIDX(i) + shiftIDX;
        newTransp = pixAreas(i,:);
        oldTransp = imgTransparency(localIDX);        
        imgTransparency(localIDX) = oldTransp + newTransp;
        
        localWts = newTransp ./ (newTransp + oldTransp);
        imgDepths(localIDX) = localWts.*Pt(i,3) + ...
            (1-localWts).*imgDepths(localIDX);
        for j=1:3
            localIDX = pixIDX(i) + (j-1)*nPix(1)*nPix(2) + shiftIDX;
            imgColors(localIDX) = localWts.*Cs(i,j) + ...
                (1-localWts).*imgColors(localIDX);
            imgNormals(localIDX) = localWts.*Ns(i,j) + ...
                (1-localWts).*imgNormals(localIDX);
        end

        if ~mod(i,rStep), progT.update(i); end
    end
    progT.done;
    
    %Mask out the pixels contiguous with the border of the image, turning on the rest
    emptyL = bwlabel(imgTransparency<10*eps);
    edgeGroups = nonzeros(unique([emptyL(1,:),emptyL(end,:),emptyL(:,1)',emptyL(:,end)']));
    imgMask = ~ismember(emptyL,edgeGroups);
    imgMask = imdilate(imgMask,strel('disk',2));
    imgMask = ~imdilate(~imgMask,strel('disk',2));
        
    %Flood gaps in scan with MAT values, when available
    floodMask = imgTransparency<10*eps & imgMask;
    MATMask = mindepthMap > minPt(3);
    imgMask(floodMask & ~MATMask) = false;
    floodMask = floodMask & MATMask;
    imgDepths(floodMask) = mindepthMap(floodMask);    
    for j=1:3
        imgNormals(:,:,j) = regionfill(imgNormals(:,:,j),floodMask);
        imgColors(:,:,j) = regionfill(imgColors(:,:,j),floodMask);
    end
    imgTransparency = imgTransparency/median(nonzeros(imgTransparency));
    imgTransparency(imgTransparency>1) = 1;
    imgTransparency = imgTransparency.^0.25; %redistribute values toward 1
    imgTransparency(floodMask) = 0.5;
    
    smN = medianSmoothing(normalizeCells(imgNormals),2);
    smN(repmat(~imgMask,1,1,3)) = 0;
    
    img2_5D.Mask = imgTransparency > 10*eps;
    img2_5D.Transparency = imgTransparency;
    img2_5D.Colors = uint8(imgColors);
    img2_5D.Normals = normalizeCells(smN);
    img2_5D.Depths = imgDepths;
    
    fprintf(repmat(char(8),1,delChars));
    fprintf(' Image complete.\n');
    
    %Save the transformed point cloud, trimmed to only the visible points
    trimmedCloud = pointCloud(Pt,'Color',Cs,'Normal',Ns);
end

function xyzidx = getPointIDX(pts,npix,ctr,pixw)
    xyzidx = (pts-ctr)/pixw; xyzidx(:,2) = -xyzidx(:,2);
    xyzidx = ceil(npix/2) + xyzidx;
end

function [pixidx,pixarea,shift] = mapPixels(ptidx,npix)
    ptDiff = ptidx - floor(ptidx);
    %Contribution of each 3D point to four adjacent pixels
    xSpan = [mod(0.5+ptDiff(:,1),1),mod(0.5-ptDiff(:,1),1)];
    ySpan = [mod(0.5+ptDiff(:,2),1),mod(0.5-ptDiff(:,2),1)];
    pixidx = sub2ind(npix,floor(ptidx(:,1)),floor(ptidx(:,2)));
    pixarea = [xSpan(:,1).*ySpan(:,1),... %LL areas
        xSpan(:,2).*ySpan(:,1), xSpan(:,1).*ySpan(:,2),... %LR, UL
        xSpan(:,2).*ySpan(:,2)]; %UR
    shift = [0, 1, npix(1), npix(1)+1];    
end

function [normSmooth,wts] = medianSmoothing(norm,pad)
    normSmooth = norm;
    for dim=1:3
        normSmooth(:,:,dim) = medfilt2(norm(:,:,dim),[pad pad]);
    end
    normSmooth = normalizeCells(normSmooth);
    
    wts = zeros(size(norm,1),size(norm,2));
    for i=1:size(norm,1)
        wts(i,:) = sum(squeeze(norm(i,:,:).*normSmooth(i,:,:)),2);
    end
    wts(wts>1) = 1;
    %Any normal more than 45 degrees from the smoothed set is weighted zero
    wts = (45-acosd(wts))/45;
    wts(wts<0) = 0;
    %Sharper normals are kept when they are close to the median values
    normSmooth = norm.*repmat(wts,1,1,3) + normSmooth.*repmat(1-wts,1,1,3);
    normSmooth = normalizeCells(normSmooth);
end