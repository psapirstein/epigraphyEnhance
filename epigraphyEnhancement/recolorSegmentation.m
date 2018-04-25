function [imout] = recolorSegmentation(seg,rdata,ridx)
    if nargin < 2
        depth = 3;
        rdata = rand(max(seg.labels(:)),3);
    else
        depth = size(rdata,2);
    end
    if nargin < 3, ridx = seg.regions; end
    
    incolor = false;
    if depth == 3, incolor = true; end
    
    imsz = size(seg.labels);
    imout = zeros([imsz,depth]);
    jshift = prod(imsz);
    for i=1:length(ridx)
        ir = ridx(i);
        segPix = seg.pixelCoords{ir};
        if ~isempty(segPix)
            for j=1:depth
                imout(segPix) = rdata(ir,j);
                segPix = segPix+jshift;
            end
        end
    end
    if incolor, imout = uint8(255*imout); end
end