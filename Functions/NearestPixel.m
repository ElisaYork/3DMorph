function [ closestPt ] = NearestPixel( Img,StartingPixel,scale )
%UNTITLED Summary of this function goes here
%   OUTPUT: x y z coordinates of the nearest pixel to the StartingPoint
%   which has a value of 1
%   INPUT: 
%   Img: Black and White image where pixels of interest have value 1 in.
%   StartingPixel = x y z coordinates of pixel of interest
%   scale = x y scale to adjust 


s = size(Img);
List = find(Img==1); %get x y z coordinates of all points in skel 2 == 1
[x,y,z]=ind2sub([s(1) s(2) s(3)],List);
List1=[x y z];

minlength=inf;
    for i=1:length(List1);
     dist= sqrt((((List1(i,2)-StartingPixel(1))*scale)^2)+(((List1(i,1)-StartingPixel(2))*scale)^2)+(((List1(i,3)-StartingPixel(3))*1)^2));
        if dist<minlength;
            minlength=dist;
            closestPt = List1(i,:);
        end
    end

closestPt = [closestPt(1),closestPt(2),closestPt(3)];

end

