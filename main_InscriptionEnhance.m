clear variables;
%Operates on a PLY formatted point cloud, which must include normals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Specify file and parameters here
basepath = './'; %Specify the path for this file
filename = ''; %The filename of the scan, without the .ply suffix
incisionWidth = 0.010; %Specify in the units of the scan
samplesPerPixel = 4; %Approximate number of scan points per pixel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath([basepath,'epigraphyEnhancement'],'-end');
fprintf('Loading scan, recentering. ');
PCfull = pcread([filename,'.ply']);
if isempty(PCfull.Normal)
    error('Error: a point cloud with precomputed normals is required.');
end

%Reorient the scan at the origin, with the major plane up on the Z
%axis, and with the maximum width along X.
[PCfull,pixW] = planarReorient(PCfull,samplesPerPixel);
%Medial balls approximate the interior and exterior of the object
MATobj = medialBallEstimate(PCfull);
%Resampling the inscription in a 2.5D image space
[PCimg,PCfull] = rasterize2_5D(PCfull,MATobj,pixW,incisionWidth);
save([filename,'-PointCloud.mat'],'PCfull','-v7.3');
PCimg = filterRIMLS2_5D(PCimg,PCfull);
%Restore the original depth of the surface before incisions or breakage
PCimg = modelInscribedSurface(PCimg,incisionWidth);
PCimg.name = filename;
save([filename,'-IMG.mat'],'PCimg','-v7.3');
image2_5Ddisplay(PCimg);