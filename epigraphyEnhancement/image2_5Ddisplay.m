classdef image2_5Ddisplay < handle
properties
    handles
    imgBorder
    imgDepth
    imgRelief
    imgData
    userVals
end

methods
    function obj = image2_5Ddisplay(img)
        %Load the image color data and mask
        imgcol = permute(img.Colors,[2 1 3]);
        imgmask = img.Mask';
        if any(strcmp('faceRestored',fieldnames(img)))
            imgborder = addProfiles(img.faceProfile,img.faceRestored,img.PixD);
            imgmask = renderBordered(imgmask,true(size(imgborder(:,:,1))));
            imgcol = renderBordered(imgcol,imgborder);
        else
            imgborder = 1;
        end
        hfig = figure('Units','normalized');
        hfig.MenuBar = 'none';
        hfig.Name = [img.name,' viewer'];
        warning('off','images:initSize:adjustingMag');
        h.image = imshow(imgcol);
        warning('on','images:initSize:adjustingMag');
        colormap('gray');
        h.image.CDataMapping = 'direct';
        h.image.AlphaData = imgmask;
        hfig.Position = [0.1 0.1 0.8 0.8];

        %% Variables the user can modify through the interface
        uv.DRweight = 0.5; %Weighting between depth and relief
        uv.inscriptionDepth = img.incisionDepth; %Inscription depth setting
        uv.pngCount = 1;
        uv.colorMode = true;
        uv.rotation = 0;
        %Default light settings for the shader
        ldir = [0 0];
        uv.lighting.dir = ldir;
        uv.lighting.ambient = 0.05;
        uv.lighting.specular = 0.50;
        uv.lighting.alpha = 0.3;

        %% Create UI elements for manipulating and saving the figure
        scrollP = imscrollpanel(hfig, h.image);
        bgColor = [0.98 0.98 0.98];
        controlpanel = uipanel(hfig,'Units','pixels','Position',[40 18 163 68],...
            'BackgroundColor',bgColor,'BorderType','beveledout');
        uicontrol('Parent',controlpanel, 'Style','text',...
            'String','Zoom %', 'FontSize',8, 'BackgroundColor',bgColor,...
            'Units','pixels', 'Position', [72 50 45 15]);
        zoomBox = uicontrol('Parent',controlpanel, 'Style','edit',...
            'String','100', 'Value',100,...
            'Position', [120 50 40 16], 'Callback', {@obj.changeZoom,scrollP});
        % The first buttons are always visible
        h.layerSet = uicontrol('Parent',controlpanel, 'Style','popup',...
            'String',{'Colors','Normals','Enhanced','Shaded'},...
            'Position',[65 18 102 28], 'Callback', @obj.updateLayer);
        uicontrol('Parent',controlpanel, 'Style','pushbutton',...
            'String','Save PNG #1', 'Position', [0 0 85 25],...
            'Callback', @obj.savePNG);
        uicontrol('Parent',controlpanel, 'Style','pushbutton',...
            'String','Gray mode', 'Position', [85 0 77 25],...
            'Callback', @obj.toggleRGB);
        % Rotation buttons based on an icon
        [~,~,rotIcon] = imread('roticon.png');
        RI = imresize(double(rotIcon)/255, [32 32]);
        RI(RI>1) = 1; RI(RI<0) = 0;
        RI = repmat(1-RI,1,1,3);
        uicontrol('Parent',controlpanel,'Style','pushbutton','Position',...
            [3 30 32 32],'CData',RI,'Callback',{@obj.rotateImg,1,scrollP});
        uicontrol('Parent',controlpanel,'Style','pushbutton','Position',...
            [35 30 32 32], 'CData',rot90(flipud(RI),2),...
            'Callback', {@obj.rotateImg,-1,scrollP});
        
        % Enhanced image control panels are shown when activated
        bgColor = [0.92 0.92 0.92];
        h.enhancedPanel = uipanel('Parent',hfig, 'BackgroundColor',bgColor,...
            'Units','pixels', 'Position',[203 18 92 68], 'BorderType','beveledout');
        %[0-1] for balancing between the Depth & Relief layers
        h.opacityTitle = uicontrol('Parent',h.enhancedPanel, 'Style','text',...
            'FontSize',9, 'BackgroundColor',bgColor, 'Position', [0 50 88 15]);
        h.opacityTitle.String = ['Relief strength: ',num2str(uv.DRweight,'%.1f')];
        uicontrol('Parent',h.enhancedPanel, 'Style','slider',...
            'Min',0,'Max',1,'Value',uv.DRweight, 'SliderStep',[0.1 0.3],...
            'Position',[5 28 80 20], 'Callback', @obj.slideOpacity);
        % Incision depth text entry
        uicontrol('Parent',h.enhancedPanel, 'Style','text', 'String','Incision Depth',...
            'FontSize',8, 'BackgroundColor',bgColor, 'Position',[0 4 40 25]);
        h.inscriptionDepth = uicontrol('Parent',h.enhancedPanel, 'Style','edit',...
            'String',num2str(uv.inscriptionDepth,'%.1e'), 'Value',uv.inscriptionDepth,...
            'Position',[38 5 48 22], 'Callback', @obj.changeIncision);
        h.enhancedPanel.Visible = 'off';
        
        % Directional lighting panel
        h.lightingPanel = uipanel('Parent',hfig, 'BackgroundColor',bgColor,...
            'Units','pixels','Position',[203 18 200 68],'BorderType','beveledout');
        lBpx = 65;
        lB = axes(h.lightingPanel,'Units','pixels','Position',[1 1 lBpx lBpx]);
        [lBpix,lBline,~] = drawLightBall(lBpx,ldir);
        hold on;
        h.lightBall = imagesc(lB,lBpix, 'ButtonDownFcn',@obj.updateLightBall,...
            'BusyAction','cancel');
        h.lightLine = plot(lB,lBline(:,1),lBline(:,2),'Color','b');
        xlim(lB,[1 lBpx]); ylim(lB,[1 lBpx]);
        uicontrol('Parent',h.lightingPanel, 'Style','text','FontSize',8,...
            'HorizontalAlignment','left', 'BackgroundColor',bgColor,...
            'Position',[67 52 70 15], 'String', 'Light direction');
        % Specular lighting controller
        h.specRatio = uicontrol('Parent',h.lightingPanel, 'Style','text',...
            'Position',[143 30 70 15], 'HorizontalAlignment','left',...
            'FontSize',8, 'BackgroundColor',bgColor);
        uicontrol('Parent',h.lightingPanel, 'Style','slider','Min',0,...
            'Max',1,'Value',uv.lighting.specular,'SliderStep',[0.05 0.25],...
            'Position',[70 30 70 15], 'Callback', {@obj.updateSpecular,1});
        h.specRatio.String = ['Specular ',...
            num2str(round(100*uv.lighting.specular),'%u'),'%'];
        h.specPower = uicontrol('Parent',h.lightingPanel, 'Style','text',...
            'Position',[143 3 70 15], 'HorizontalAlignment','left',...
            'FontSize',8, 'BackgroundColor',bgColor);
        uicontrol('Parent',h.lightingPanel, 'Style','slider','Min',0.1,...
            'Max',0.9,'Value',uv.lighting.alpha, 'SliderStep',[0.125 0.25],...
            'Position',[70 3 70 15], 'Callback', {@obj.updateSpecular,2});
        h.specPower.String = ['Power: ',...
            num2str(round(0.5*2^(10*uv.lighting.alpha)),'%u')];
        h.lightingPanel.Visible = 'off';
        
        obj.handles = h;
        obj.userVals = uv;
        
        imgdata = struct;
        imgdata.color = permute(img.Colors,[2 1 3]);
        imgdata.ambientOcclusion = img.ambientOcclusion';
        imgdata.mask = img.Mask';
        imgdata.normals = permute(img.Normals,[2 1 3]);
        imgdata.relativeNormals = permute(img.relativeNormals,[2 1 3]);
        imgdata.relativeDepth = img.originalPlane' - img.Depths';
        imgdata.RSI = prepareRSI(img.meanCurvature');
        if sum(size(imgborder))>3
            imgdata.faceProfile = img.faceProfile;
            imgdata.faceRestored = img.faceRestored;
        end
        imgdata.initialInscriptionDepth = uv.inscriptionDepth;
        imgdata.name = img.name;
        imgdata.pixD = img.PixD;
        
        obj.imgData = imgdata;
        obj.imgBorder = imgborder;
        obj.imgDepth = renderDepth(obj);
        imgRelief = img.ambientOcclusion.*img.relativeNormals.^2;
        expon = max(log(0.7)/log(median(nonzeros(imgRelief))),0.5);
        obj.imgRelief = (imgRelief.^expon)';
        changeZoom(obj,zoomBox,0,scrollP);
        clear -regexp ^img;
    end
    
    function changeZoom(~,src,~,scrollp)
        apiSP = iptgetapi(scrollp);
        newzoom = str2double(src.String);
        if ~isnan(newzoom)
            newzoom = round(newzoom);
            if newzoom<5, newzoom = 5; end
            if newzoom>2000, newzoom = 2000; end
            apiSP.setMagnification(newzoom/100);
            src.Value = newzoom;
        end
        src.String = num2str(src.Value,'%u');
    end

    function toggleRGB(obj,src,~)
        colmode = ~obj.userVals.colorMode;
        obj.userVals.colorMode = colmode;
        if colmode, src.String = 'Gray mode';
        else, src.String = 'Color mode'; end
        updateLayer(obj,obj.handles.layerSet,0);
    end

    function changeIncision(obj,src,~)
        newID = str2double(src.String);
        if ~isnan(newID)
            %The original inscription depth is saved in the 2.5D image
            origID = obj.imgData.initialInscriptionDepth;
            if ~newID, newID = origID; end %Entering 0 resets value
            if newID < origID/10, newID = origID/10; end
            if newID > origID*10, newID = origID*10; end
            mag = 2-floor(log(newID)/log(10));
            newID = round(newID*10^mag)/10^mag; %Round to 3 effective digits

            obj.userVals.inscriptionDepth = newID;
            obj.handles.inscriptionDepth.String = num2str(newID,'%.1e');
            obj.handles.inscriptionDepthBW.String = num2str(newID,'%.1e');
            obj.imgDepth = renderDepth(obj);
            updateLayer(obj,obj.handles.layerSet,0);
        end
    end
    
    function slideOpacity(obj,src,~)
        newval = round(src.Value*10)/10;
        if newval ~= obj.userVals.DRweight
            obj.handles.opacityTitle.String = ['Relief strength: ',...
                num2str(newval,'%.1f')];
            obj.userVals.DRweight = newval;
            updateLayer(obj,obj.handles.layerSet,0);
        end
    end
    
    function updateLightBall(obj,LBim,hitpt)
        newpt = hitpt.IntersectionPoint(1:2);
        [LBim.CData,newline,newldir] = drawLightBall(length(LBim.CData),newpt);
        obj.handles.lightLine.delete;
        obj.handles.lightLine = plot(obj.handles.lightBall.Parent,...
            newline(:,1), newline(:,2), 'Color','b');
        obj.userVals.lighting.dir = newldir;
        updateLayer(obj,obj.handles.layerSet,0);
    end

    function updateSpecular(obj,src,~,mode)
        redraw = false;
        if mode==1 %ratio to ambient
            newval = round(src.Value*20)/20;
            if newval<0, newval = 0; end
            if newval>1, newval = 1; end
            if newval ~= obj.userVals.lighting.specular
                obj.userVals.lighting.specular = newval;
                obj.handles.specRatio.String = ['Specular ',...
                    num2str(round(100*newval),'%u'),'%'];
                redraw = true;
            end
        else %alpha slider
            newval = round(src.Value*10)/10;
            if newval>0.8, newval = 0.9; end
            if newval<0.1, newval = 0.1; end
            if newval ~= obj.userVals.lighting.alpha
                obj.userVals.lighting.alpha = newval;
                obj.handles.specPower.String = ['Power: ',...
                    num2str(round(0.5*2^(10*newval)),'%u')];
                redraw = true;
            end
        end
        if redraw, updateLayer(obj,obj.handles.layerSet,0); end
    end
    
    function savePNG(obj,src,~)
        oUV = obj.userVals;
        PPU = num2str(1/obj.imgData.pixD,'%.0f');
        hLS = obj.handles.layerSet;
        fname = [obj.imgData.name,' ',num2str(oUV.pngCount,'%d'),' (',...
            PPU,' ppu) ', hLS.String{hLS.Value}];
        if hLS.Value == 3 %Append enhancement settings
            ID = num2str(oUV.inscriptionDepth,'%.1e');
            DS = num2str(round(oUV.DRweight*10)/10,'%.1f');
            fname = [fname,' (Incision Depth ',ID, ', Relief Weight ',DS,')'];
        end
        if hLS.Value == 4 %Specular
            SP = num2str(round(oUV.lighting.specular*100),'%u');
            PW = num2str(round(0.5*2^(10*oUV.lighting.alpha)),'%u');
            fname = [fname,' (Specular 0.',SP,', power ',PW,')'];
        end
        fname = [fname,'.png'];

        src.String = ['Save PNG #',num2str(oUV.pngCount+1,'%d')];
        obj.userVals.pngCount = oUV.pngCount+1;

        hIm = obj.handles.image;
        imwrite(hIm.CData,fname,'png','Alpha',double(hIm.AlphaData));
    end
    
    function rotateImg(obj,~,~,dir,scrollp)
        obj.userVals.rotation = mod(obj.userVals.rotation + dir,4);
        imgmask = rot90(obj.imgData.mask,-obj.userVals.rotation);
        if sum(size(obj.imgBorder))>3        
            fP = cellfun(@transpose,flip(obj.imgData.faceProfile),'un',0);
            fR = cellfun(@transpose,flip(obj.imgData.faceRestored),'un',0);
            if dir < 0, fP{1} = flip(fP{1}); fR{1} = flip(fR{1});
            else, fP{2} = flip(fP{2}); fR{2} = flip(fR{2}); end
            obj.imgData.faceProfile = fP;
            obj.imgData.faceRestored = fR;
            obj.imgBorder = addProfiles(fP, fR, obj.imgData.pixD);
            imgmask = renderBordered(imgmask,true(size(obj.imgBorder(:,:,1))));
        end
        apiSP = iptgetapi(scrollp);
        apiSP.replaceImage(obj.imgBorder,'PreserveView',true);
        obj.handles.image.AlphaData = imgmask;
        
        updateLayer(obj,obj.handles.layerSet,0);
    end
    
    function updateLayer(obj,src,~)
        colmode = obj.userVals.colorMode;
        imD = obj.imgData;
        obj.handles.enhancedPanel.Visible = 'off';
        obj.handles.lightingPanel.Visible = 'off';
        switch src.Value
            case 1
                cdata = imD.color;
                if ~colmode, cdata = repmat(uint8(mean(double(cdata),3)),1,1,3);
                end
            case 2
                cdata = imD.normals(:,:,3).*(imD.relativeNormals.^1.5);
                if colmode, cdata = uint8(cdata.*double(imD.color));
                else, cdata = repmat(uint8(255*cdata),1,1,3);
                end
            case 3
                if ~colmode, cdata = updateEnhancedBW(obj);
                else, cdata = updateEnhancedRGB(obj); end
                obj.handles.enhancedPanel.Visible = 'on';
            case 4
                if colmode, cdata = shadeImage(imD.color,...
                        imD.normals, obj.userVals.lighting);
                else, cdata = shadeImage(255, imD.normals,...
                        obj.userVals.lighting);
                end
                obj.handles.lightingPanel.Visible = 'on';
        end
        cdata = rot90(cdata, -obj.userVals.rotation);
        hIm = obj.handles.image;
        hIm.CData = renderBordered(cdata, obj.imgBorder);
    end
end
end

function RSI = prepareRSI(meanC)
    %Curvature with Radiance Scaling; see Vergne et al 2008, appendix
    A = max(abs(meanC(:)));
    B = 8;
    K = atanh((2^B-2)/(2^B-1));
    RSI = tanh(-meanC*K/A);
    RSI(RSI>0) = RSI(RSI>0)/3;
    RSI = 0.75*(RSI+1);
end

function imgout = renderBordered(img,border)
    if sum(size(border))>3
        imgNx = 1+size(border,1)-size(img,1);
        imgNy = 1+size(border,2)-size(img,2);
        imgout = border;
        imgout(imgNx:end, imgNy:end, :) = img;
    else
        imgout = img;
    end
end

function [lBpix,lBline,ldir] = drawLightBall(lBpx,ldir)
    %Draw the ball slider for the light
    lBctr = ceil(lBpx/2); lBpad = 2;
    lBscale = lBctr-lBpad;
    
    XYmax = 0.99; %Always leave a small positive Z
    %Convert ldir from pixels to coordinates in the unit circle
    ldir = (ldir-lBctr)/lBscale;
    if sum(ldir.^2)>1, ldir = ldir/sqrt(sum(ldir.^2)); end
    if sum(ldir.^2)>XYmax, ldir = ldir*XYmax; end
    ldir(2) = -ldir(2);
    %Draw the light direction on a simulated sphere
    lBpix = 1 -0.5*(repmat( (ldir(2) + ((1:lBpx)-lBctr)/lBscale).^2, lBpx, 1)...
        + repmat( (ldir(1) + (lBctr-(1:lBpx))'/lBscale).^2, 1, lBpx) );
    lBpix(lBpix<0.2) = 0.2;
    %Mask out the background
    lBdist = repmat(((lBctr - (1:lBpx))/lBscale).^2,lBpx,1);
    lBdist = 1-(lBdist + lBdist');
    lBpix(lBdist<0) = 0.8;
    lBpix = lBpix';
    %Save the line vector
    lBline = [lBctr, lBctr; ldir(1)*lBscale+lBctr, lBctr-ldir(2)*lBscale];
end

function imgout = shadeImage(col,nm,lighting)
    %Light direction is supplied as normalized X,Y component, missing Z
    ldir = lighting.dir;
    ldir(3) = sqrt(1-sum(ldir.^2));
    amb = lighting.ambient;
    spec = (1-amb)*lighting.specular;
    diff = 1-amb-spec;
    A = 0.5*2^(10*lighting.alpha);
    
    H = ldir+[0,0,1];
    H = H/sqrt(sum(H.^2));
    ldir = ldir*diff;
    
    nm = permute(nm,[2 1 3]);
    specM = H(1)*nm(:,:,1) + H(2)*nm(:,:,2) + H(3)*nm(:,:,3);
    specM(specM<0) = 0;
    imgout = amb + spec*(specM.^A) + ldir(1)*nm(:,:,1) + ...
        ldir(2)*nm(:,:,2) + ldir(3)*nm(:,:,3);
    imgout = uint8(double(col) .* repmat(imgout',1,1,3));
end

function DS = renderDepth(obj)
    %Generate the depth-enhanced image, such that 1 = uninscribed,
    %0 = inscribed depth, -1 = unknown / broken
    inD = obj.userVals.inscriptionDepth;
    DS = min((inD-obj.imgData.relativeDepth)/inD,1);
    DS(DS<0) = max(DS(DS<0)/4,-1); %Expanding range of shading in breaks
end

function RGBenhanced = updateEnhancedBW(obj)    
    dw = obj.userVals.DRweight; %Sharpening for blended images
    DS = max(obj.imgDepth,-0.5);
    DS = (1-dw)*(abs(DS).^(1+dw/2)) + dw*(obj.imgRelief.^(1.5-dw/2));
    RGBenhanced = uint8(repmat(255*DS,1,1,3));
end

function [RGBenhanced] = updateEnhancedRGB(obj)
    dw = obj.userVals.DRweight/2; %Blending intensity favors relief
    DS = obj.imgDepth;
    RL = obj.imgData.RSI.^2 .* obj.imgRelief; %RSI included in color
    expon = max(log(0.4)/log(median(nonzeros(RL))),0.5);
    RL = RL.^expon;
    RL = (0.5-dw)*(0.6+0.3*DS) + (0.5+dw)*(RL.^(1.5-dw));
    RL = min(RL/prctile(nonzeros(RL),98),1);
    RL(~obj.imgData.mask) = 0;
    
    HSV = ones([size(DS) 3]); %HSV for remapping depths, such that:
    %green = uninscribed; red = inscribed depth; blue = unknown / broken
    HSV(:,:,1) = mod(0.98+0.35*sign(obj.imgDepth).*abs(DS).^1.2,1);
    HSV(:,:,2) = (1-1.5*dw)*(0.6+0.2*DS); %Lower areas are desaturated
    HSV(:,:,3) = RL;
    RGBenhanced = uint8(255*hsv2rgb(HSV));
end

function bimg = addProfiles(fProfile,fTrend,pixD)
    zFac = 4; %Expand the depth range by a factor of 4
    xBorder = drawProfile(fProfile{1},fTrend{1},pixD,zFac);
    yBorder = rot90(drawProfile(fProfile{2},fTrend{2},pixD,zFac),-1);
    xD = size(xBorder); yD = size(yBorder);
    borderIMG = ones(xD(1)+yD(1)+10, yD(2)+xD(2)+10, 3);
    borderIMG(end-xD(1):end,end-yD(2):end,:) = 0; %Black in image area
    borderIMG(yD(1)+11:end,1:xD(2),:) = xBorder; %Transfer XY profiles
    borderIMG(1:yD(1),xD(2)+11:end,:) = yBorder;
    borderIMG([1,yD(1)],xD(2)+1:xD(2)+10,:) = 0; %Draw ticks at min/max values
    borderIMG(yD(1)+1:yD(1)+10,[1,xD(2)],:) = 0;
    borderIMG(yD(1)+10,:,:) = 0; %Draw a frame around the image area
    borderIMG(:,xD(2)+10,:) = 0;
    
    borderIMG = permute(borderIMG,[2 1 3]);
    %Write the z-scaling into the image border
    bimg = insertText(borderIMG,[2,2], [num2str(zFac,'%d'),'x depth'],...
        'AnchorPoint','LeftTop', 'FontSize',14, 'TextColor','black',...
        'BoxColor','white', 'BoxOpacity',0);
    bimg = (bimg + borderIMG)/2;
    bimg = uint8(255*bimg);   
end

function bimg = drawProfile(fP,rP,PD,zF)
    function [px,wt] = wtpx(pin)
        px = [floor(pin),ceil(pin)];
        wt = repmat(abs(px-pin)',1,2);
    end
    RS = 4.0; %Rescaling factor
    fP = fP(:); rP = rP(:);
    dRange = [min([fP;rP]), max([fP;rP])];
    if dRange(1)>0 %Function shifts profiles that do not cross the origin
        fP = fP-dRange(1);
        rP = rP-dRange(1);
        dRange = dRange-dRange(1);
    elseif dRange(2)<0
        fP = fP-dRange(2);
        rP = rP-dRange(2);
        dRange = dRange-dRange(2);
    end
    dDpix = zF*(dRange(2)-dRange(1))/PD;
    dDpix = ceil(dDpix)+1;
    d0 = dDpix-zF*dRange(2)/PD;
    bimg = ones(length(fP),dDpix,3);
    [px0,wt0] = wtpx(d0);
    %Medium green marking the zero line
    bimg(:,px0(1),[1 3]) = repmat(0.5+0.5*wt0(1,1),length(fP),1,2);
    bimg(:,px0(2),[1 3]) = repmat(0.5+0.5*wt0(2,1),length(fP),1,2);
    for i=1:length(fP)
        [fPx,fWt] = wtpx(d0 + zF*fP(i)/PD);
        [tPx,tWt] = wtpx(d0 + zF*rP(i)/PD);
        %Divvy up amounts per pixel
        bimg(i,fPx,2:3) = fWt; %Face profile shown in red
        bimg(i,tPx,1:2) = tWt; %Trend shown in blue
    end
    %The loop will emphasize horizontals at the expense of verticals, so a
    %bresenham line is generated to fill in the verticals of the face profile.
    %The trend line, which is already highly smoothed, is not altered.
    Bt = zeros(RS*length(fP), RS*dDpix);
    fP1 = RS*(d0+zF*fP(1)/PD);
    for i=1:length(fP)-1
        fP2 = RS*(d0+zF*fP(i)/PD);
        [fPx,fPy] = bresenham(RS*(i-0.5),fP1,RS*(i+0.5),fP2);
        Bt(fPx,fPy) = 1;
        fP1 = fP2;
    end
    Bt = imresize(Bt,1/RS,'bilinear');
    bimg(:,:,2:3) = bimg(:,:,2:3) - repmat(Bt,1,1,2);
    bimg(bimg>1) = 1; bimg(bimg<0) = 0;
end

function [x,y]=bresenham(x1,y1,x2,y2)
    dx=abs(x2-x1); dy=abs(y2-y1);
    steep=abs(dy)>abs(dx);
    if steep, t=dx;dx=dy;dy=t; end
    if dy==0
        q=zeros(dx+1,1);
    else
        q=[0;diff(mod((floor(dx/2):-dy:-dy*dx+floor(dx/2))',dx))>=0];
    end
    if steep
        x=x1+cumsum(q);
        if y1<=y2
            y=(y1:y2)';
        else
            y=(y1:-1:y2)';
        end
    else
        x=(x1:x2)';
        if y1<=y2
            y=y1+cumsum(q);
        else
            y=y1-cumsum(q);
        end
    end
    x = round(x); y = round(y);
end