function [ mask ] = ConnectPointsAlongPath( InputImg, StartPt, EndPt )
%ConnectPointsAlongPath: Traces path between two points while constrained to 
%path of pixels with value one.

%Find the shortest path connecting the centroid and the endoiints.
%This script duplicates the image, and finds the path (following only
%along line of ones) from endpoint to center. In the second image, it
%finds the path from center to endpoint. When these two lines are
%added, the shortest path will be the one == n. Any other paths of
%connection (that are longer) will be erased. The output image is
%mask, that has ones where a valid shortest connection exists.

%OUTPUT: mask is an image the same size as input image. It contains 1s
%where the shortest path was discovered. 

%INPUT: 
%InputImg must be a black and white image where the path will lie
%along pixels of value 1. At the moment, it only works with 3D, but could
%easily be modified. 
%StartPt is the x y z coordinate of one pixel of interest. Note, this pixel
%must have a value of 1 in the original image. Check this, it could be
%causing your errors.
%EndPt is the second point of interest to be connected. It must also be of
%value 1 in the original image. 

    init = [nan inf]; % set or randomize initial parameter values. nan = not-a-number, inf = infinity.
    D = repmat(init(1+InputImg),[1 1 1 2]); %repmat repeats copies of a matrix. This duplicates the BoundedSkel image. It also multiplies it by init, so all 0s become Nan and all 1s become Inf
    D(StartPt(1),StartPt(2),StartPt(3),1) = 0; %set the initial and end points to 0.
    D(EndPt(1),EndPt(2),EndPt(3),2) = 0;
    mask = D==0; %matrix with 1s (true) where D is 0. All branch pixels are Inf.
    n = 0; %n is the number of steps taken from starting point.
    while isinf(D(EndPt(1),EndPt(2),EndPt(3)))&nnz(mask),% continue stepping while the endpoint (but in the first image) is inf. Once the path is connected, this value will be overwritten to 1, meaning the loop is complete. nnz = number of nonzero elements. If the path is not completed, but there are no remaining connections, the mask will only have 0s. Both of these conditions must be 1 (true) to continue.
        n = n+1;
        mask = convn(mask,ones([3 3 3]),'same')&isinf(D); %Convolves with cube of ones to see attached pixels?
        D(mask) = n;
    end
    if isinf(D(EndPt(1),EndPt(2),EndPt(3))) %If the above loop ended, it's complete or there's no path. If it ended, but the first condition is still true (the endpoint is still Inf), then there must not have been a connecting path.
        error('no path found');
    else
        mask = sum(D,4)==n; % points within minimum-length path(s). sums the two 3D matrices. Only those paths ==n are the shortest and the rest are ignored.
        %sumBW = sum(BoundedSkel(mask)); %Optionally, can calculate the
        %length of the connected mask in pixels.
    end

end

