function [xc,yc] = circleCenterMask(mask)

bw = mask;
[Ny, Nx] = size(bw);

cc = bwconncomp(bw);                   % finds connected components
stats = regionprops(cc, 'Centroid');   % calculates centroid for each component
centroids = cat(1, stats.Centroid);   

img_center = [Nx/2, Ny/2];
dist2 = sum((centroids - img_center).^2, 2); % Component distance from image center
[~, idxSphere] = min(dist2);                 % component closer to the center (circle)

% Creates new mask with sphere only
bw_sphere = false(size(bw));
bw_sphere(cc.PixelIdxList{idxSphere}) = true;
bw = bw_sphere;

props = regionprops(bw, 'Centroid');
center = props.Centroid;  

xc = center(1); 
yc = center(2);  



