function [ BoundedImg, right, left, top, bottom ] = BoundingBoxOfCell( InputImg )
%MakeBoundingBoxAroundCell 
%From a black and white image, Find the pixels of value 1 (assuming these
%are the object of interest. Find the outer most points in x and y
%(bounding points) and use these locations to shrink the image so it
%encloses only the cell with no extra white space.
%   OUTPUT: the coordinates of min and max in x and y = top, bottom, left,
%   and right coordinates of cell. These define our bounding limitations. 
%   BoundedImg is an output image of only the object of interest without
%   any excess empty borders.
%   INPUT: InputImg must be a balck and white (1s and 0s) image with 3
%   dimensions. The object of interest should be in 1s.

    s = size(InputImg);
    ImageList = find(InputImg==1); %get x y z coordinates of all points in skeleton == 1
    [x,y,z]=ind2sub([s(1) s(2) s(3)],ImageList);
    ImageListMatrix=[x y z];

    right = max(ImageListMatrix(:,1)); %Find the min max x and y to make bounding box (include all z). This should make next loop faster.
    left = min(ImageListMatrix(:,1));
    top = max(ImageListMatrix(:,2));
    bottom = min(ImageListMatrix(:,2));
   
    BoundedImg = InputImg(left:right,bottom:top,:);%Extract bounding box area from skel2 to get smaller image, skel4

end

