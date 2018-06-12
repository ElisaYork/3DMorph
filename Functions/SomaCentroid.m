function [ CentroidCoord ] = SomaCentroid( CellImg )
%SomaCentroid Finds the coordinates of the centroid pixel of the soma
%   Erodes cell so that the centroid is not influenced by the 'weight' of
%   processes, but considers only the soma to find the true nucleus. 
%OUTPUT: matrix of x, y, z coordinates for the pixel representing the
%centroid of the soma
%INPUT: 3D array of binary image of one cell at a time

    se=strel('diamond',3);
    nucmask=imerode(CellImg,se);%Erode to find nuclei only. This may give more than one spot.
    nm=bwconncomp(nucmask,26);%use connected components to identify each spot
    nP = cellfun(@numel,nm.PixelIdxList);
    [~,idx] = max(nP);%Find the largest connected component, which is probably the true nuclei
    s = size(CellImg);
    ex1=zeros(s(1),s(2),s(3));%Create blank image of correct size
    nuc=(nm.PixelIdxList{1,idx});%select largest nuc object
    ex1(nuc)=1;%write in only one object to image. Cells are white on black background.
    c = regionprops(ex1,'centroid'); %Find the centroid location of each object (centre of nuclei only so point isn't pulled by 'weight' of processes)
    CentroidCoord = c.Centroid(1,:); %Write the coordinates to the output variable

end

