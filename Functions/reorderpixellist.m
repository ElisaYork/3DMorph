function [distpoint] = reorderpixellist(pxlist,si,ep,stp)
%reorder: finds all pixels in masklist (ie. for each branch) that are 1s,
%and organizes them by connectivity, rather than from smallest to largest.

%OUTPUT: N x 3 matrix of pixel coordinates. N1 is the endpoint, and the
%last row is the centroid. Pixels are arranged in order of connectivity
%based on the distance between them.

%INPUT: 

[mr,mc,mp]=ind2sub(si,pxlist); %Convert indexed values to x y z coordinates (row, column, plane)
distpoint = [mr mc mp];% Concatenate coordinates into single matrix
    
epind = find( (distpoint(:,1) == ep(1)) & (distpoint(:,2) == ep(2)) & (distpoint(:,3) == ep(3))); %Get row index of endpoint in distpoint list

%exchange ep with the first coordinate in the distpoint list
tmp=distpoint(1,:); %Save the current first coordinate as tmp
distpoint(1,:)=ep; %Replace the first coordinate with ep
distpoint(epind,:)=tmp;%Swap the original ep location with the previously first coordinate (so it's not just overwritten)

%Find the index of the closest pixel, and move that pixel to the next row
%in the distpoint list. 
ii = 1; 
    while ii<length(distpoint) && sum(distpoint(ii,:)==stp)~=3; %Run loop until either hit the last row of distpoint, or we hit the centroid
        dist = pdist2(distpoint(ii,:),distpoint(ii:end,:)); %distance from one point to all next pixels in list. Don't want to connect to previous pixel, so ignore all above pixels
        dist(~dist)=nan;%will have 0 distance to itself, so replace 0s with NaN to not interfere with min calculation
        [M I] = min(dist);%I is the column index of the shortest connectoin, so the location of the nearest pixel. This is how many pixels away from current (not the actual row value, bc not comparing all rows)
        I = I+ii-1; %the row of the nearest pixel will be I rows away from ii(-1 because columns start at 1, not 0).
        tmp=distpoint(ii+1,:); %save the pixel coordinates that we're about to erase in the variable tmp
        distpoint(ii+1,:)=distpoint(I,:); %replace the next pixel in the list with the nearest pixel coordinates
        distpoint(I,:)=tmp;%Put the pixel coordinates from tmp in the previous location of the nearest pixel
        ii = ii+1;
    end
    
distpoint=distpoint(1:ii,:); %Remove any pixels below the centroid. These would be 'spur' pixels that weren't needed to complete the path. 
%If leave these pixels in, then it will try to connect the centroid to them
%and give an erronously long length. 

end

