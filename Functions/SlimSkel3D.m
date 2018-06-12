function [ SlimSkel ] = SlimSkel3D( CellToSkel,THR )
%SlimSkel3D Skeletonizes and slims output to reduce residual thick branches
%or spurs.
% Uses Skeleton3D and parts of Skel2Graph3D and Graph2Skel3D
%OUTPUT: 3D binary image of slim skeleton
%INPUT: CellToSkel is a 3D binary image containing the cell we want to
%skeletonize

    skel = Skeleton3D(CellToSkel);
    [~,node,link] = Skel2Graph3D(skel,THR); % skel is the input 3D binary image, and "THR" is a threshold for the minimum length of branches, to filter out skeletonization artifacts.   
    % ~ is the adjacency matrix with the length of the links as matrix entries, and node/link are the structures describing node and link properties.
    %Graph2Skel3D converts the network graph back into a cleaned-up voxel skeleton image
    wl = sum(cellfun('length',{node.links}));%this is still in pixels 
    s=size(CellToSkel);
    SlimSkel = Graph2Skel3D(node,link,s(1),s(2),s(3));
    [~,node2,link2] = Skel2Graph3D(SlimSkel,0);
    % calculate new total length of network
    wl_new = sum(cellfun('length',{node2.links}));
    % iterate the same steps until network length changed by less than
    % 0.5%. This should remove any weird spurs on skeleton. 
while(wl_new~=wl)

    wl = wl_new;   
    
     SlimSkel = Graph2Skel3D(node2,link2,s(1),s(2),s(3));
     [~,node2,link2] = Skel2Graph3D(SlimSkel,0);

     wl_new = sum(cellfun('length',{node2.links}));

end;


end

