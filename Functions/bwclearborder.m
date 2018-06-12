function out = bwclearborder(bw,conn)
%Function to clear true objects touching the border of a binary image
%   This function will do the same thing as IMCLEARBORDER but only for 
%   binary images.  It is much faster for them.
%
%SCd 11/02/2010
%
%Input Arguments:
%   -bw: 2 or 3 dimensional binary image with parts you want true
%   
%   -conn: connectivity of connected components
%       -conn can be 4 or 8 for a 2-dimensional image (default: 8)
%       -conn can be 6, 18 or 26 for 3-dimensional images (default: 26)
%
%Output Arguments:
%   -out = bw with all components connected by conn to the border turned
%       off, i.e. = 0. 
%
%See Also: imclearborder, bwconncomp
%

%Error checking and default assignment:
assert(any(nargin==[1 2]),'The number of input arguments is expected to be 1 or 2.\n  It was %d',nargin);
assert(islogical(bw),'The first input argument, bw, is expected to be of class logical\n  It was of class %s',class(bw));
assert(any(ndims(bw)==[2 3]),'The first input argument, bw, is expected to have 2 or 3 dimensions.\n  It had %d',ndims(bw))
dims = ndims(bw);
if nargin == 1
    if dims == 2
        conn = 8;
    else %3d
        conn = 26;
    end  
else
    assert(isnumeric(conn),'The second input argument conn is expected to be numeric');    
end

%Components analysis; if it fails throw error from bwconncomp
CC = bwconncomp(bw,conn);

%Some sizes
NumObj = CC.NumObjects;
szbw = size(bw);
CCndx = false(NumObj,1);

%Loop through objects
if dims == 2
    for ii = 1:NumObj
        [r c] = ind2sub(szbw,CC.PixelIdxList{ii}); %sub indices
        CCndx(ii) = ~(any([r; c] == 1)||any(reshape(bsxfun(@eq,[r c],szbw),2*length(c),1))); %test if on low or high border
    end
else %3d
    for ii = 1:NumObj
        [r c p] = ind2sub(szbw,CC.PixelIdxList{ii});
        CCndx(ii) = ~(any([r; c; p] == 1)||any(reshape(bsxfun(@eq,[r c p],szbw),3*length(c),1)));
    end   
end
out = false(szbw); %binary image save size as bw
out(cell2mat(reshape(CC.PixelIdxList(CCndx),[],1))) = true; %set objects that survived to true
