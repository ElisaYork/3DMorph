%% Semi-automatic 3D morphology analysis: 3DMorph
%First run script in interactive mode to generate a parameters file. It may
%then be run in automatic mode (with or without figure outputs) to batch
%process a series of images. *Note: these images must have the same scaling
%factor and use the same user-defined values. Loads .tiff or .lsm z-stack
%images. 
%User input: Threshold adjustment, maximum cell size for
%segmentation, minimum cell size to remove out-of-frame processes, and
%skeletonization method. 
%Outputs: Before removing small processes: Total territorial (convex)
%volume covered, total empty volume, and percentage of volume covered. On
%only full cells: Average centroid distance between cells (cell dinsity or
%dispersion), territotrial volume, cell volume, cell complexity
%(territorial volume / cell volume), number of endpoints, branch points,
%and the minimum, maximum, and average branch length.

%% Method Selection
% delete(gcp('nocreate'));
% parpool %Open parallel processing. 

addpath(genpath('Functions'));

question = {'How would you like to run the script?','Interactive Mode is necessary to set your parameters. Use automatic if you already have a saved Parameters.mat file.'};
choiceMode = questdlg(question,'Mode Selection','Interactive Mode', 'Automatic Mode', 'Automatic Mode Without Images', 'Interactive Mode');
%Response:
switch choiceMode 
    case 'Interactive Mode'
        Interactive = 1;
        NoImages = 0;
        FileDataGUI;
        addpath(pathname);
        [~,file] = fileparts(file);
        % example: file = 'con1_CD68_2'; scale = 0.46125; %1 pixel = ___ um
        FileList = 1;
    case 'Automatic Mode'
        Interactive = 2;
        NoImages = 0;
        callgui = SelectFilesGUI;
        uiwait(callgui);        
    case 'Automatic Mode Without Images'   
        Interactive = 2;
        NoImages = 1;
        callgui = SelectFilesGUI;
        uiwait(callgui);                
end

for total = 1:numel(FileList)
%% Load file and saved values

clearvars -except file ch ChannelOfInterest scale zscale Parameters FileList PathList Interactive NoImages total

if Interactive == 2
    load(Parameters);
    addpath(PathList);
        if NoImages ==1
        ShowImg = 0; ShowObjImg = 0; ShowCells = 0; ShowFullCells = 0; ConvexCellsImage = 0; OrigCellImg = 0; SkelImg = 0; EndImg = 0; BranchImg = 0;
        end
end

if Interactive == 2
   [~,file] = fileparts(FileList{1,total});
end

%Load in image stack from either .lsm or .tif file
if exist(strcat(file,'.lsm'),'file') == 2
%use bfopen to read in .lsm files to variable 'data'. Within data, we want
%the 1st row and column, which is a list of all channels. From these, we
%want the _#_ channel (ChannelOfInterest).
data = bfopen(strcat(file,'.lsm'));

s = size(data{1,1}{1,1}); % x and y size of the image
l = length(data{1,1});
zs = l(:,1)/ch; %Number of z planes

%To read only green data from a 3 channel z-stack, will write as slice 1
%ch1, 2, 3, then slice 2 ch 1, 2, 3, etc. Need to extract every third image
%to a new array. [] concatenates all of the retrieved data, but puts them
%all in many columns, so need to reshape.
img = reshape([data{1,1}{ChannelOfInterest:ch:end,1}],s(1),s(2),zs); 
end

if exist('img','var')==0 %If both a tiff and lsm file exist, load only the lsm file
    %Open .tif file and reshape channels
    if exist(strcat(file,'.tif'),'file') == 2
       tiffInfo = imfinfo(strcat(file,'.tif'));
       no_frame = numel(tiffInfo);
       data = cell(no_frame,1);
       for iStack = 1:no_frame
       data{iStack} = imread(strcat(file,'.tif'),'Index',iStack);
       end
       data = data(ChannelOfInterest:ch:end,1);
       s = size(data{1,1});
       zs = length(data);
       img = reshape([data{:,1}], s(1), s(2), zs);
    end
end

voxscale = scale*scale*zscale;%Calculate scale to convert voxels into unit^3.
%% Threshold
%Load GUI to set thresholding parameters. 
if Interactive == 1
midslice = round(zs/2,0);
orig = img(:,:,midslice);
callgui = ThresholdGUI;
uiwait(callgui);
end

if Interactive == 2
    BinaryThresholdedImg=zeros(s(1),s(2),zs);
    for i=1:zs %For each slice
    level=graythresh(img(:,:,i)); %Finds threshold level using Otsu's method. Different for each slice b/c no need to keep intensity conistent, just want to pick up all mg processes. Tried adaptive threshold, but did not produce better images.
    level=level*adjust; %Increase the threshold by *1.6 to get all fine processes 
    BinaryThresholdedImg(:,:,i) = im2bw(img(:,:,i),level);% Apply the adjusted threshold and convert from gray the black white image.
    end
    NoiseIm=bwareaopen(BinaryThresholdedImg,noise); %Removes objects smaller than set value (in pixels). For 3D inputs, uses automatic connectivity input of 26. Don't want small background dots left over from decreased threshold.
end

ConnectedComponents=bwconncomp(NoiseIm,26); %returns structure with 4 fields. PixelIdxList contains a 1-by-NumObjects cell array where the k-th element in the cell array is a vector containing the linear indices of the pixels in the k-th object. 26 defines connectivity. This looks at cube of connectivity around pixel.
numObj = numel(ConnectedComponents.PixelIdxList); %PixelIdxList is field with list of pixels in each connected component. Find how many connected components there are.

%show full image compressed to 2D. Use imagesc to make it look 3D. 
progbar = waitbar(0,'Processing your data...');
for i = 1:numObj
    waitbar (i/numObj, progbar);
    ex=zeros(s(1),s(2),zs);
    ex(ConnectedComponents.PixelIdxList{1,i})=1;%write in only one object to image. Cells are white on black background.
    flatex = sum(ex,3);
    allObjs(:,:,i) = flatex(:,:); 
end
if isgraphics(progbar)
close(progbar);
end
DetectedObjs = sum(allObjs,3);
cmapCompress = parula(max(DetectedObjs(:)));  
cmapCompress(1,:) = zeros(1,3);

if ShowImg == 1    
    title = [file,'_Threshold (compressed to 2D)'];
    figure('Name',title);imagesc(DetectedObjs);
    colormap(cmapCompress);
    daspect([1 1 1]);
end

    %Extract list of pixel values and which object they belong to for
    %segmentation viewing (in CellSizeCutoffGUI). 
    for i = 1:numObj
    ObjectList(i,1) = length(ConnectedComponents.PixelIdxList{1,i}); 
    ObjectList(i,2) = i;  
    end
    ObjectList = sortrows(ObjectList,-1);%Sort columns by pixel size.
%     ObjectList = sortrows(ObjectList,'descend'); %Sort columns by pixel size. 
    udObjectList = flipud(ObjectList);%ObjectList is large to small, flip upside down so small is plotted first in blue.

%% Cell Segmentation 
% Decreased threshold may cause cells to be inappropriately connected.
% Choose the threshold for a large cell, and segment larger objects into
% separate cells (by identifying number of nuclei and running fitgmdist
% (fit Gaussian mixture distribution). Function is run 3 times to improve
% accuracy and replicability

if Interactive == 1
    callgui = CellSizeCutoffGUI;
    waitfor(callgui);
end

if ShowObjImg == 1
    if Interactive ~= 1
        num = [1:3:(3*numObj+1)];
        fullimg = ones(s(1),s(2));
        progbar = waitbar(0,'Plotting...');
        for i = 1:numObj;
            waitbar (i/numObj, progbar);
            ex=zeros(s(1),s(2),zs);
            j=udObjectList(i,2);
            ex(ConnectedComponents.PixelIdxList{1,j})=1;%write in only one object to image. Cells are white on black background.
            flatex = sum(ex,3);
            OutlineImage = zeros(s(1),s(2));
            OutlineImage(flatex(:,:)>1)=1;
            se = strel('diamond',4);
            Outline = imdilate(OutlineImage,se); 
            fullimg(Outline(:,:)==1)=1;
            fullimg(flatex(:,:)>1)=num(1,i+1);
        end
        if isgraphics(progbar)
           close(progbar);
        end
    end
    cmapnumObj = jet(max(fullimg(:)));  
    cmapnumObj(1,:) = zeros(1,3);
    title = [file,'_Selected Objects (compressed to 2D)'];
    figure('Name',title);imagesc(fullimg);
    colormap(cmapnumObj); 
    colorbar('Ticks',[1,3*numObj+1], 'TickLabels',{'Small','Large'});
    daspect([1 1 1]);
end

col=1;
    progbar = waitbar(0,'Segmenting...');
for i = 1:numObj %Evaluate all connected components in PixelIdxList.
    waitbar (i/numObj, progbar);
    if  numel(ConnectedComponents.PixelIdxList{1,i}) > CellSizeCutoff %If the size of the current connected component is greater than our predefined cutoff value, segment it.
        ex=zeros(s(1),s(2),zs);%Create blank image of correct size.
        ex(ConnectedComponents.PixelIdxList{1,i})=1;%write object onto blank array so only the cell pixels = 1.
        se=strel('diamond',6); %Set how much, and what shape, we want to erode by. If increase, erosion will be greater.
        nucmask=imerode(ex,se);%Erode to find nuclei only. This erosion is large - don't want any remaining thick branch pieces.
        nucsize = round((CellSizeCutoff/50),0);
        nucmask=bwareaopen(nucmask,nucsize);%Get rid of leftover tiny spots. Increase second input argument to remove more spots (and decrease segmentation)
        indnuc=bwconncomp(nucmask);%Find these connected comonents (number of remaining nuclei).
        nuc = numel(indnuc.PixelIdxList);%Determine the number of nuclei (how many objects to segment into).
        if nuc ==0 %If erosion only detects one nuc, but this should be segmented, increase nuc to at least 2
            error('Error: The program finds 0 nuclei to segement your object into. Adjust se=strel(diamond,4) to a lower number to decrease image erosion.');
        end  
        if nuc ==1 %If erosion only detects one nuc, but this should be segmented, increase nuc to at least 2
            nuc = 2;
        end
        [x,y,z]=ind2sub(size(ex),find(ex));%Find nonzero elements in ex (ie connected microglia cells) and return x y z locations.
        points = [x y z]; %concatenate to one array
        
        GMModel = fitgmdist(points,nuc,'replicates',3); %Fit Gaussian mixture distribution to data 
        idx = cluster(GMModel,points);
        
        for j=1:nuc %Extract location of all pixels for each object.
            obj = (idx == j); % |1| for cluster 1 membership
%             obj = find(clust==j);%return linear index of all values == nuc.
            object = points(obj,1:3); %find x y z data of given linear index.
            objidx=sub2ind(size(NoiseIm),object(:,1),object(:,2),object(:,3)); %convert x y z data to index in the same size as ex.
            ex=zeros(s(1),s(2),zs);%Create blank image of correct size
            ex(objidx)=1;%write in only one object to image. Cells are white on black background.
            ex=bwareaopen(ex,noise);
            individual=bwconncomp(ex,26);
            NumberOfIdentifiedObjects = length(individual.PixelIdxList);
                for m = 1:NumberOfIdentifiedObjects
                Microglia{1,col}=individual.PixelIdxList{1,m}; %Write this separated object to new cell array, Microglia in location 'col'.
                col=col+1; %Increase the col counter so data is not overwritten.
                end
        clear('individual');        
        end
        clear('da','NucCoord','CenterOfNuc','indnuc','nuc','x','y','z','points');
    else %If the size of the current connected component is NOT greater than our predefined cutoff value, keep the current component and rewrite to new cell array.
        Microglia{1,col}=ConnectedComponents.PixelIdxList{1,i}; %Write the object to new cell array, Microglia in location 'col'.
        col=col+1;%Increase the col counter so data is not overwritten.
    end
end

numObjSep = numel(Microglia); %Rewrite the number of objects to include segmented cells.

%Extract list of pixel values and which object they belong to for
%segmentation viewing (in FullCellsGUI).
    for i = 1:numObjSep
    SepObjectList(i,1) = length(Microglia{1,i}); 
    SepObjectList(i,2) = i;  
    end
    SepObjectList = sortrows(SepObjectList,-1); %Sort columns by pixel size. 
    udSepObjectList = flipud(SepObjectList);%ObjectList is large to small, flip upside down so small is plotted first in blue.

%Below is used in FullCellsGUI
    for i = 1:numObjSep
        ex=zeros(s(1),s(2),zs);
        ex(Microglia{1,i})=1;%write in only one object to image. Cells are white on black background.
        flatex = sum(ex,3);
        AllSeparatedObjs(:,:,i) = flatex(:,:); 
    end
if isgraphics(progbar)
close(progbar);
end
 
if Interactive == 1
    question = {'Would you like to see your separated cells?','Warning: A 3D image will take a moment to process! We may downsample it for you...'};
    choiceSegmentImg = questdlg(question,'Output Segmented Image?','3D image please!', '2D image please!', 'No thanks','No thanks');
    %Response:
    switch choiceSegmentImg 
        case '3D image please!'
            ShowCells = 1;
        case '2D image please!'
            ShowCells = 2;    
        case 'No thanks'
            ShowCells = 0;
    end
end

num = [1:3:(3*numObjSep+1)];
cmap = jet(max(num));
cmap(1,:) = zeros(1,3);

if ShowCells == 1
    title = [file,'_Selected Cells'];
    figure('Name',title);
    colormap(cmap);
    progbar = waitbar(0,'Plotting...');
    for i = 1:numObjSep
        waitbar (i/numObjSep, progbar);
        ex=zeros(s(1),s(2),zs);%Create blank image of correct size
        j=udSepObjectList(i,2);
        ex(Microglia{1,j})=1;%write in only one object to image. Cells are white on black background.
        % If it's larger than 512x512, downsample image to increase processing time (this is ONLY to display, doesn't change actual figure). 
        if s(1)>=1024||s(2)>= 1024
            ex = imresize(ex,0.5);
            ex = (ex(:,:,1:2:end));
        end
        if 512>= s(1)&& s(1)<1024||512>=s(2)&& s(2)<1024
            ex = imresize(ex,0.5);
            ex = (ex(:,:,1:2:end));
        end
        ds = size(ex);
        fv=isosurface(ex,0);%display each object as a surface in 3D. Will automatically add the next object to existing image.
        patch(fv,'FaceColor',cmap(i*3,:),'FaceAlpha',1,'EdgeColor','none');%without edgecolour, will auto fill black, and all objects appear black
        axis([0 ds(1) 0 ds(2) 0 ds(3)]);%specify the size of the image
        camlight %To add lighting/shading
        lighting gouraud; %Set style of lighting. This allows contours, instead of flat lighting
        view(0,270); % Look at image from top viewpoint instead of side  
        daspect([1 1 1]);
        colorbar('Ticks',[0,1], 'TickLabels',{'Small','Large'});
    end
    if isgraphics(progbar)
    close(progbar);
    end
end

if ShowCells == 2
    fullimg = ones(s(1),s(2));
    progbar = waitbar(0,'Plotting...');
        for i = 1:numObjSep
            waitbar (i/numObjSep, progbar);
            ex=zeros(s(1),s(2),zs);
            j=udSepObjectList(i,2);
            ex(Microglia{1,j})=1;%write in only one object to image. Cells are white on black background.
            flatex = sum(ex,3);
            OutlineImage = zeros(s(1),s(2));
            OutlineImage(flatex(:,:)>1)=1;
            se = strel('diamond',4);
            Outline = imdilate(OutlineImage,se); 
            fullimg(Outline(:,:)==1)=1;
            fullimg(flatex(:,:)>1)=num(1,i+1);
        end
            if isgraphics(progbar)
            close(progbar);
            end
    title = [file,'_Selected Cells (compressed to 2D)'];
    figure('Name',title);imagesc(fullimg);
    colormap(cmap);
    colorbar('Ticks',[1,max(num)], 'TickLabels',{'Small','Large'});
    daspect([1 1 1]);    
end

%% Territorial volume 
%Uses convhulln to create a 3D polygon around the object's external points.
%For the total occupied vs unoccupied volume, don't want to exclude any
%cells/processes. Use Microglia list here, not FullMg.

ConvexVol = zeros(numObjSep,1);
sz = size(img);
progbar = waitbar(0,'Finding territorial volume...');
for i = 1:numObjSep
    waitbar (i/numObjSep, progbar);
    [x,y,z] = ind2sub(sz,[Microglia{1,i}]); %input: size of array ind values come from, list of values to convert.
    obj = [y,x,z]; %concatenate x y z coordinates.
    [k,v] = convhulln(obj);
    ConvexVol(i,:) = v*voxscale;
end

if isgraphics(progbar)
close(progbar);
end

TotMgVol = sum(ConvexVol); %Calculate total volume of image covered by microglia. 
CubeVol = (s(1)*s(2)*zs)*voxscale; %volume of image cube in um^3. 
EmptyVol = CubeVol-TotMgVol;%And the remaining 'empty space'.
PercentMgVol = ((TotMgVol)/(CubeVol))*100;

%% Full cells
%Option to remove all cells touching the x y border, and all objects below
%size limit. From here on in code, will only be looking at these full
%cells.

if Interactive == 1
%     Don't need this anymore bc I now generate the image inside the GUI
%     if ShowCells == 2
%         %do nothing, the figure is already made
%     else
%         num = 1:3:(3*numObjSep+1);
%         num(1,1) = zeros(1,1);
%             fullimg = ones(s(1),s(2));
%     progbar = waitbar(0,'Plotting...');
%         for i = 1:numObjSep
%             waitbar (i/numObjSep, progbar);
%             ex=zeros(s(1),s(2),zs);
%             j=udSepObjectList(i,2);
%             ex(Microglia{1,j})=1;%write in only one object to image. Cells are white on black background.
%             flatex = sum(ex,3);
%             OutlineImage = zeros(s(1),s(2));
%             OutlineImage(flatex(:,:)>1)=1;
%             se = strel('diamond',4);
%             Outline = imdilate(OutlineImage,se); 
%             fullimg(Outline(:,:)==1)=1;
%             fullimg(flatex(:,:)>1)=num(1,i+1);
%         end
%             if isgraphics(progbar)
%             close(progbar);
%             end
%     end    
callgui2 = FullCellsGUI;
waitfor(callgui2);
end

if KeepAllCells == 1
    FullMg = Microglia;
else   
    col=1;
    progbar = waitbar(0,'Finding All Full Cells...');         
    for i = 1:numObjSep
        waitbar (i/numObjSep, progbar);
        if numel(Microglia{1,i})>=SmCellCutoff %Extract elements of Microglia that are >SmCellCutoff pixels
            if RemoveXY ==1
            ex=zeros(s(1),s(2),zs);%Create blank image of correct size
            ex(Microglia{1,i})=1;%plot it onto original blank image
            antiborder = logical(padarray(ex,[0 0 1],0));%add row of zeros to top and bottom in z axis so no objects are touching this border. Also make logical.
            cleared = bwclearborder(antiborder,26); %Like imclearborder, but MUCH faster! Removes any 1 objects touching border in 26 connectivity (so whole object)
            nonedge = max(cleared(:)); %If the object was removed, this should be 0, so don't add it to the FullMg array
                if nonedge == 1 % If the cell is not touching an x or y edge, keep it in new cell array. 
                   FullMg{1,col} = (Microglia{1,i}); %Suppress not preallocated error.
                   col = col+1;
                end
            else
            FullMg{1,col} = (Microglia{1,i}); %Suppress not preallocated error.
            col = col+1;   
            end  
        end
    end
    if isgraphics(progbar)
    close(progbar);
    end
end
numObjMg = numel(FullMg);% number of microglia after excluding edges and small processes.

%Extract list of pixel values and which object they belong to for
%segmentation viewing (in FullCellsGUI).
for i = 1:numObjMg
MgObjectList(i,1) = length(FullMg{1,i}); 
MgObjectList(i,2) = i;  
end
MgObjectList = sortrows(MgObjectList,-1); %Sort columns by pixel size. 
udMgObjectList = flipud(MgObjectList);%ObjectList is large to small, flip upside down so small is plotted first in blue.


num = [1:3:(3*numObjMg+1)];
cmap = jet(max(num));
cmap(1,:) = zeros(1,3);
    
if Interactive == 1
    %See all full cells?
    question2 = {'Would you like to see remaining full cells?','Warning: A 3D image will take a moment to process! We may downsample it for you...'};
    choiceFullCellsImg = questdlg(question2,'Output Full Cell Image?','3D image please!', '2D image please!', 'No thanks','No thanks');
    %Response:
    switch choiceFullCellsImg 
        case '3D image please!'
            ShowFullCells = 1;
        case '2D image please!'
            ShowFullCells = 2;    
        case 'No thanks'
            ShowFullCells = 0;
    end
end
    
if ShowFullCells == 1
    progbar = waitbar(0,'Plotting...');
    title = [file,'_Full Cells'];
    figure('Name',title);
    colormap(cmap);
    for i = 1:numObjMg
        waitbar (i/numObjMg, progbar);
        ex=zeros(s(1),s(2),zs);%Create blank image of correct size
        j=udMgObjectList(i,2);
        ex(FullMg{1,j})=1;%write in only one object to image. Cells are white on black background.
        % If it's larger than 512x512, downsample image to increase processing time (this is ONLY to display, doesn't change actual figure). 
        if s(1)>=1024||s(2)>= 1024
            ex = imresize(ex,0.5);
            ex = (ex(:,:,1:2:end));
        end
        if 512>= s(1)&& s(1)<1024||512>=s(2)&& s(2)<1024
            ex = imresize(ex,0.5);
            ex = (ex(:,:,1:2:end));
        end
        ds = size(ex);
        fv=isosurface(ex,0);%display each object as a surface in 3D. Will automatically add the next object to existing image.
        patch(fv,'FaceColor',cmap(i*3,:),'FaceAlpha',1,'EdgeColor','none');%without edgecolour, will auto fill black, and all objects appear black
        axis([0 ds(1) 0 ds(2) 0 ds(3)]);%specify the size of the image
        camlight %To add lighting/shading
        lighting gouraud; %Set style of lighting. This allows contours, instead of flat lighting
        view(0,270); % Look at image from top viewpoint instead of side  
        daspect([1 1 1]);
        colorbar('Ticks',[0,1], 'TickLabels',{'Small','Large'});
    end
    if isgraphics(progbar)
    close(progbar);
    end
end

if ShowFullCells == 2
    fullimg = ones(s(1),s(2));
    progbar = waitbar(0,'Plotting...');
        for i = 1:numObjMg
            waitbar (i/numObjMg, progbar);
            ex=zeros(s(1),s(2),zs);
            j=udMgObjectList(i,2);
            ex(FullMg{1,j})=1;%write in only one object to image. Cells are white on black background.
            flatex = sum(ex,3);
            OutlineImage = zeros(s(1),s(2));
            OutlineImage(flatex(:,:)>1)=1;
            se = strel('diamond',4);
            Outline = imdilate(OutlineImage,se); 
            fullimg(Outline(:,:)==1)=1;
            fullimg(flatex(:,:)>1)=num(1,i+1);
        end
            if isgraphics(progbar)
            close(progbar);
            end
    title = [file,'_Full Cells (compressed to 2D)'];
    figure('Name',title);imagesc(fullimg);
    colormap(cmap);
    colorbar('Ticks',[1,max(num)], 'TickLabels',{'Small','Large'});
    daspect([1 1 1]);    
end

%% Volume of Full Cells
% Determines the cell volume by the number of voxels multiplied by the
% voxscale to convert into real world units. Also finds the convex
% territorial volume of only full cells. Cell complexity or extent) is
% calculated as territorial volume / cell volume and represents how
% bushy/amoeboid or branched cells are within their territory.

NumberOfPixelsPerCell = cellfun(@numel,FullMg);
[biggest,idx] = max(NumberOfPixelsPerCell);
CellVolume = (NumberOfPixelsPerCell*voxscale)'; %list of volume of each cell
MaxCellVol = biggest*voxscale; %Volume determined by microscope scale, and voxel number reported in cc.PixelIdxList. Should be in um^3

if Interactive == 1
    ConvexImgQuestion = {'Would you like to see a convex volume image of each full cell?'};
    choiceConvexVolImg = questdlg(ConvexImgQuestion,'Output Convex Volume Images?','Yes please!', 'No thanks','No thanks');
    %Response
    switch choiceConvexVolImg 
        case 'Yes please!'
            ConvexCellsImage = 1;
        case 'No thanks'
            ConvexCellsImage = 0;    
    end
end

% Find the convexvol of only full cells
FullCellTerritoryVol = zeros(numObjMg,1);
for i = 1:numObjMg
    [x,y,z] = ind2sub(sz,[FullMg{1,i}]); %input: size of array ind values come from, list of values to convert.
    obj = [y,x,z]; %concatenate x y z coordinates.
    [k,v] = convhulln(obj);
    FullCellTerritoryVol(i,:) = v*voxscale;
    if ConvexCellsImage == 1
    figure;
    trisurf(k,obj(:,1),obj(:,2),obj(:,3));
    axis([0 s(1) 0 s(2) 0 zs]);
    daspect([1 1 1]);
    end
end

% Find the complexity (or extent) of full cells.
FullCellComplexity = zeros(numObjMg,1);
for i = 1:length(FullCellTerritoryVol)
FullCellComplexity(i,:) = FullCellTerritoryVol(i)/CellVolume(i);
end

%% Distance Between Centroids
%Finds location of centroid of each cell, and measures distance between
%them as a measure of cell density or dispersion.

cent = (zeros(numObjMg,3));
progbar = waitbar(0,'Finding centroids...');
for i=1:numObjMg
    waitbar (i/numObjMg, progbar);
    ex=zeros(s(1),s(2),zs);%Create blank image of correct size
    ex(FullMg{1,i})=1;%write in only one object to image. Cells are white on black background.
    CentroidCoord = SomaCentroid(ex);
    cent(i,:) = CentroidCoord;
end

if isgraphics(progbar)
close(progbar);
end

centum = (zeros(numObjMg,3));
centum(:,1)=cent(:,1)*scale; %Convert pixel location to microns so that distances are in correct scale
centum(:,2)=cent(:,2)*scale;
centum(:,3)=cent(:,3)*zscale;
centdist = pdist2(centum,centum); %Calculate distance from each centroid to all other centroids
centdist=nonzeros(centdist); %Remove all 0s (distance from one centroid to itself)
AvgDist = mean(centdist);

%% 3D Skeleton
%Can use two different methods to keep all processes (more fine
%structures), or only major branches (and ignore small extensions). From
%skeleton, endpoints and branch points are identified. In skeleton image,
%red processes are primary, yellow secondary, green tertiary, and blue are
%connected to endpoints. Branch lengths are measured as the shortest
%distance from each endpoint to the centroid. If the program is unable to
%propoerly identify a centroid or endpoints, it will output a 0 and move to
%the next cell.

kernel(:,:,1) = [1 1 1; 1 1 1; 1 1 1];
kernel(:,:,2) = [1 1 1; 1 0 1; 1 1 1];
kernel(:,:,3) = [1 1 1; 1 1 1; 1 1 1]; 

numendpts = zeros(numel(FullMg),1);
numbranchpts = zeros(numel(FullMg),1);
MaxBranchLength = zeros(numel(FullMg),1);
MinBranchLength = zeros(numel(FullMg),1);
AvgBranchLength = zeros(numel(FullMg),1);

if Interactive == 1
    %Use Skeleton method to include small processes/fillipodia?
    question3 = ['Would you like to include all small processes or fillipodia in your skeleton analysis? ','Note: only keep small processes if you want all fine fillipodia - this will increase processing time.'];
    choiceSkelMethod = questdlg(question3,'Skeletonization Method','Keep small processes', 'Only major branches', 'Only major branches');
    %Response:
    switch choiceSkelMethod 
        case 'Keep small processes'
            SkelMethod = 1;
        case 'Only major branches'
            SkelMethod = 2;    
    end
end

if Interactive == 1
    callgui = SkeletonImageGUI;
    uiwait(callgui);
end

if SkelImg||EndImg||BranchImg||OrigCellImg||BranchLengthFile ==1
    folder = mkdir ([file, '_figures']); 
    fpath =([file, '_figures']);
end 

if BranchLengthFile == 1;
BranchLengthList=cell(1,numel(FullMg));
end

parfor i=1:numel(FullMg)
    try
    ex=zeros(s(1),s(2),zs);%#ok<PFBNS> %Create blank image of correct size
    ex(FullMg{1,i})=1;%write in only one object at a time to image. 
    ds = size(ex); 
    
    if OrigCellImg == 1
        title = [file,'_Cell',num2str(i)];
        figure('Name',title);
        fv=isosurface(ex,0);%display each object as a surface in 3D. Will automatically add the next object to existing image.
        patch(fv,'FaceColor',cmap(i,:),'FaceAlpha',1,'EdgeColor','none');%without edgecolour, will auto fill black, and all objects appear black
        axis([0 ds(1) 0 ds(2) 0 ds(3)]);%specify the size of the image
        camlight %To add lighting/shading
        lighting gouraud; %Set style of lighting. This allows contours, instead of flat lighting
        view(0,270); % Look at image from top viewpoint instead of side  
        daspect([1 1 1]);
        filename = ([file '_Original_cell' num2str(i)]);
        saveas(gcf, fullfile(fpath, filename), 'jpg');
    end
    
    if SkelMethod == 1
        WholeSkel = SlimSkel3D(ex,100);
        DownSampled = 0;
        adjust_scale = scale;
    end
    
    if SkelMethod == 2
        if s(1)>512 %convert large cells to 512x512 to speed up skeletonization. The branch lengths are later adjusted to account for this down-sampling.
            ex = imresize(ex,0.5);
            adjust_scale = 2*scale;
            DownSampled = 1;
        else 
            DownSampled = 0;
            adjust_scale = scale;
        end

        SmoothEx = imgaussfilt3(ex); %Smooth the cell so skeleton doesn't pick up many fine hairs

        FastMarchSkel = skeleton(SmoothEx);%Find the skeleton! This uses msfm3d and rk4 files, which have been compiled and the .mexw64 versions included. If errors, re-run compilation of these files (in FastMarching_version3b folder), and add the folder and subfolders to path. 

        %Convert cell output of branches into one image for further processing.
          WholeSkel=zeros(size(ex));
          WholeList = round(vertcat(FastMarchSkel{:}));
          SkelIdx = sub2ind(size(ex),WholeList(:,1),WholeList(:,2),WholeList(:,3));
          WholeSkel(SkelIdx)=1;
    end
    
      [BoundedSkel, right, left, top, bottom]  = BoundingBoxOfCell(WholeSkel); %Create a bounding box around the skeleton and only analyze this area to significantly increase processing speed. 
      si = size(BoundedSkel);

% Find endpoints, and trace branches from endpoints to centroid    
    i2 = floor(cent(i,:)); %From the calculated centroid, find the nearest positive pixel on the skeleton, so we know we're starting from a pixel with value 1.
    if DownSampled == 1
       i2(1) = round(i2(1)/2);
       i2(2) = round(i2(2)/2);
    end
    closestPt = NearestPixel(WholeSkel,i2,scale);
    i2 = closestPt; %Coordinates of centroid (endpoint of line).
    i2(:,1)=(i2(:,1))-left+1;
    i2(:,2) = (i2(:,2))-bottom+1;

    endpts = (convn(BoundedSkel,kernel,'same')==1)& BoundedSkel; %convolution, overlaying the kernel cube to see the sum of connected pixels.      
    EndptList = find(endpts==1);
    [r,c,p]=ind2sub(si,EndptList);%Output of ind2sub is row column plane
    EndptList = [r c p];
    numendpts(i,:) = length(EndptList);

    masklist =zeros(si(1),si(2),si(3),length(EndptList));
    ArclenOfEachBranch = zeros(length(EndptList),1);
    for j=1:length(EndptList)%Loop through coordinates of endpoint.
        i1 = EndptList(j,:); 
        mask = ConnectPointsAlongPath(BoundedSkel,i1,i2);
        masklist(:,:,:,j)=mask;
        % Find the mask length in microns
        pxlist = find(masklist(:,:,:,j)==1);%Find pixels that are 1s (branch)
        distpoint = reorderpixellist(pxlist,si,i1,i2); %Reorder pixel lists so they're ordered by connectivity
        %Convert the pixel coordinates by the scale to calculate arc length in microns.
        distpoint(:,1) = distpoint(:,1)*adjust_scale; %If 1024 and downsampled, these scales have been adjusted
        distpoint(:,2) = distpoint(:,2)*adjust_scale; %If 1024 and downsampled, these scales have been adjusted
        distpoint(:,3) = distpoint(:,3)*zscale;
        [arclen,seglen] = arclength(distpoint(:,1),distpoint(:,2),distpoint(:,3));%Use arc length function to calculate length of branch from coordinates
        ArclenOfEachBranch(j,1)=arclen; %Write the length in microns to a matrix where each row is the length of each branch, and each column is a different cell.
    end
  
    %Find average min, max, and avg branch lengths
    MaxBranchLength(i,1) = max(ArclenOfEachBranch);
    MinBranchLength(i,1) = min(ArclenOfEachBranch);
    AvgBranchLength(i,1) = mean(ArclenOfEachBranch);  
    
    %Save branch lengths list
    if BranchLengthFile == 1
       BranchLengthList{1,i} = ArclenOfEachBranch;
    end   
    
    fullmask = sum(masklist,4);%Add all masks to eachother, so have one image of all branches.
    fullmask(fullmask(:,:,:)>3)=4;%So next for loop can work, replace all values higher than 3 with 4. Would need to change if want more than quaternary connectivity.

    % Define branch level and display all on one colour-coded image.
    pri = (fullmask(:,:,:))==4;
    sec = (fullmask(:,:,:))==3;
    tert = (fullmask(:,:,:))==2;
    quat = (fullmask(:,:,:))==1;
    
    if SkelImg == 1
    title = [file,'_Cell',num2str(i)];
    figure('Name',title); %Plot all branches as primary (red), secondary (yellow), tertiary (green), or quaternary (blue). 
    hold on
    fv1=isosurface(pri,0);%display each object as a surface in 3D. Will automatically add the next object to existing image.
    patch(fv1,'FaceColor',[1 0 0],'FaceAlpha',0.5,'EdgeColor','none');%without edgecolour, will auto fill black, and all objects appear black
    camlight %To add lighting/shading
    lighting gouraud; %Set style of lighting. This allows contours, instead of flat lighting
    fv1=isosurface(sec,0);%display each object as a surface in 3D. Will automatically add the next object to existing image.
    patch(fv1,'FaceColor',[1 1 0],'FaceAlpha',0.5,'EdgeColor','none');%without edgecolour, will auto fill black, and all objects appear black
    camlight %To add lighting/shading
    lighting gouraud; %Set style of lighting. This allows contours, instead of flat lighting
    fv1=isosurface(tert,0);%display each object as a surface in 3D. Will automatically add the next object to existing image.
    patch(fv1,'FaceColor',[0 1 0],'FaceAlpha',0.5,'EdgeColor','none');%without edgecolour, will auto fill black, and all objects appear black
    camlight %To add lighting/shading
    lighting gouraud; %Set style of lighting. This allows contours, instead of flat lighting
    fv1=isosurface(quat,0);%display each object as a surface in 3D. Will automatically add the next object to existing image.
    patch(fv1,'FaceColor',[0 0 1],'FaceAlpha',0.5,'EdgeColor','none');%without edgecolour, will auto fill black, and all objects appear black
    camlight %To add lighting/shading
    lighting gouraud; %Set style of lighting. This allows contours, instead of flat lighting
    view(0,270); % Look at image from top viewpoint instead of side
    daspect([1 1 1]);
    hold off
    filename = ([file '_Skeleton_cell' num2str(i)]);
    saveas(gcf, fullfile(fpath, filename), 'jpg');
    end
    
    % Find branchpoints
    brpts =zeros(si(1),si(2),si(3),4);
    for kk=1:3 %For branchpoints not connected to end branches (ie. not distal branches). In fullmask, 1 is branch connected to end point, so anything greater than that is included. 
    temp = (fullmask(:,:,:))>kk;
    tempendpts = (convn(temp,kernel,'same')==1)& temp; %Get all of the 'distal' endpoints of kk level branches
    brpts(:,:,:,kk+1)=tempendpts;
    end

    % Find any branchpoints of 1s onto 4s (ie. final branch coming off of main trunk). 
    quatendpts = (convn(quat,kernel,'same')==1)& quat; %convolution, overlaying the kernel cube onto final branches only.
    quatbrpts = quatendpts - endpts; %Have points at both ends of final branches. Want to exclude any distal points (true endpoints)
    %Only want to keep these quant branchpoints if they're connected to a 4(primary branch). Otherwise, the branch point will have been picked up in the previous for loop. 
    fullrep= fullmask;
    fullrep(fullrep(:,:,:)<4)=0;%Keep only the 4s, as 4s (don't convert to 1)
    qbpts = fullrep+quatbrpts;%Add the two vectors, so should have 4s and 1s.
    qbpts1 = convn(qbpts,ones([3 3 3]),'same'); %convolve with cube of ones to get 'connectivity'. All 1s 
    brpts(:,:,:,1) = (quatbrpts.*qbpts1)>= 5;

    allbranch = sum(brpts,4); %combine all levels of branches
    BranchptList = find(allbranch==1);%Find how many pixels are 1s (branchpoints)
    [r,c,p]=ind2sub(si,BranchptList);%Output of ind2sub is row column plane
    BranchptList = [r c p];
    numbranchpts(i,:) = length(BranchptList);
    
    if EndImg == 1
        title = [file,'_Cell',num2str(i)];
        figure('Name',title); %Plot all branches with endpoints
        fv1=isosurface(fullmask,0);%display each object as a surface in 3D. Will automatically add the next object to existing image.
        patch(fv1,'FaceColor',[0 0 1],'FaceAlpha',0.1,'EdgeColor','none');%without edgecolour, will auto fill black, and all objects appear black
        camlight %To add lighting/shading
        lighting gouraud; %Set style of lighting. This allows contours, instead of flat lighting
        view(0,270); % Look at image from top viewpoint instead of side 
        fv2=isosurface(endpts,0);%display each object as a surface in 3D. Will automatically add the next object to existing image.
        patch(fv2,'FaceColor',[1 0 0],'EdgeColor','none');%without edgecolour, will auto fill black, and all objects appear black
        camlight %To add lighting/shading
        lighting gouraud; %Set style of lighting. This allows contours, instead of flat lighting
        view(0,270);
        daspect([1 1 1]);
        filename = ([file '_Endpoints_cell' num2str(i)]);
        saveas(gcf, fullfile(fpath, filename), 'jpg');
    end
    
    if BranchImg == 1
        title = [file,'_Cell',num2str(i)];
        figure('Name',title); %Plot all branches with branchpoints
        fv1=isosurface(fullmask,0);%display each object as a surface in 3D. Will automatically add the next object to existing image.
        patch(fv1,'FaceColor',[0 0 1],'FaceAlpha',0.1,'EdgeColor','none');%without edgecolour, will auto fill black, and all objects appear black
        camlight %To add lighting/shading
        lighting gouraud; %Set style of lighting. This allows contours, instead of flat lighting
        view(0,270); % Look at image from top viewpoint instead of side 
        fv2=isosurface(allbranch,0);%display each object as a surface in 3D. Will automatically add the next object to existing image.
        patch(fv2,'FaceColor',[1 0 0],'EdgeColor','none');%without edgecolour, will auto fill black, and all objects appear black
        camlight %To add lighting/shading
        lighting gouraud; %Set style of lighting. This allows contours, instead of flat lighting
        view(0,270);
        daspect([1 1 1]);
        filename = ([file '_Branchpoints_cell' num2str(i)]);
        saveas(gcf, fullfile(fpath, filename), 'jpg');
    end
        
   catch
       %do nothing if an error is detected, just write in zeros and
       %continue to next loop iteration. 
   end
   
   disp(['cell ' num2str(i) ' of ' num2str(numel(FullMg))]); %#ok<PFBNS> %To see which cell we are currently prcoessing.
end

%Save Branch Lengths File
 if BranchLengthFile == 1
     names = ["cell1"];
     %Write in headings
    for CellNum = 1:numel(FullMg)
        input = strcat('Cell ',num2str(CellNum));
        names(CellNum,1) = input;
    end   
    BranchFilename = 'BranchLengths';
    xlswrite(fullfile(fpath, BranchFilename),names(:,:),1,'A1');
    %Write in data
    for ColNum = 1:numel(FullMg)
        if numel(BranchLengthList{1,ColNum})>0
            xlswrite(fullfile(fpath, BranchFilename),BranchLengthList{1,ColNum}',1,['B' num2str(ColNum)]);
        end
    end
 end
 
%% Output results
%Creates new excel sheet with file name and saves to current folder.

xlswrite((strcat('Results',file)),{file},1,'B1');
xlswrite((strcat('Results',file)),{'Avg Centroid Distance um'},1,'A2');
xlswrite((strcat('Results',file)),AvgDist,1,'B2');
xlswrite((strcat('Results',file)),{'TotMgTerritoryVol um3'},1,'A3');
xlswrite((strcat('Results',file)),TotMgVol,1,'B3');
xlswrite((strcat('Results',file)),{'TotUnoccupiedVol um3'},1,'A4');
xlswrite((strcat('Results',file)),EmptyVol,1,'B4');
xlswrite((strcat('Results',file)),{'PercentOccupiedVol um3'},1,'A5');
xlswrite((strcat('Results',file)),PercentMgVol,1,'B5');
xlswrite((strcat('Results',file)),{'CellTerritoryVol um3'},1,'D1');
xlswrite((strcat('Results',file)),FullCellTerritoryVol(:,1),1,'E');
xlswrite((strcat('Results',file)),{'CellVolumes'},1,'F1');
xlswrite((strcat('Results',file)),CellVolume(:,1),1,'G');
xlswrite((strcat('Results',file)),{'RamificationIndex'},1,'H1');
xlswrite((strcat('Results',file)),FullCellComplexity(:,1),1,'I');
xlswrite((strcat('Results',file)),{'NumOfEndpoints'},1,'J1');
xlswrite((strcat('Results',file)),numendpts(:,1),1,'K');
xlswrite((strcat('Results',file)),{'NumOfBranchpoints'},1,'L1');
xlswrite((strcat('Results',file)),numbranchpts(:,1),1,'M');
xlswrite((strcat('Results',file)),{'AvgBranchLength'},1,'N1');
xlswrite((strcat('Results',file)),AvgBranchLength(:,1),1,'O');
xlswrite((strcat('Results',file)),{'MaxBranchLength'},1,'P1');
xlswrite((strcat('Results',file)),MaxBranchLength(:,1),1,'Q');
xlswrite((strcat('Results',file)),{'MinBranchLength'},1,'R1');
xlswrite((strcat('Results',file)),MinBranchLength(:,1),1,'S');

if Interactive == 2
disp(['Finished file ' num2str(total) ' of ' num2str(numel(FileList))]);
end

handles=findall(0,'type','figure');

for fig = 1:numel(handles) 
filename = get(handles(fig),'Name');
saveas(handles(fig), fullfile(fpath, filename), 'jpg');
end
close all;

end

%% Parameters file
%Save .mat parameters file for batch processing. Name is "Parameters_file
%name_date(year month day hour)"

if Interactive == 1
time = clock;
name = ['Parameters_',file,'_',num2str(time(1)),num2str(time(2)),num2str(time(3)),num2str(time(4))];
save(name,'ch','ChannelOfInterest','scale','zscale','adjust','noise','s','ShowImg','ShowObjImg','ShowCells','ShowFullCells','CellSizeCutoff','SmCellCutoff','KeepAllCells','RemoveXY','ConvexCellsImage','SkelMethod','SkelImg','OrigCellImg','EndImg','BranchImg','BranchLengthFile');
end

delete(gcp); %close parallel pool so error isn't generated when program is run again.