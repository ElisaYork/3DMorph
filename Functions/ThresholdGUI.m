function varargout = ThresholdGUI(varargin)
% THRESHOLDGUI MATLAB code for ThresholdGUI.fig
%      THRESHOLDGUI, by itself, creates a new THRESHOLDGUI or raises the existing
%      singleton*.
%
%      H = THRESHOLDGUI returns the handle to a new THRESHOLDGUI or the handle to
%      the existing singleton*.
%
%      THRESHOLDGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in THRESHOLDGUI.M with the given input arguments.
%
%      THRESHOLDGUI('Property','ThreshValue',...) creates a new THRESHOLDGUI or raises the
%      existing singleton*.  Starting from the left, property threshvalue pairs are
%      applied to the GUI before ThresholdGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid threshvalue makes property application
%      stop.  All inputs are passed to ThresholdGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ThresholdGUI

% Last Modified by GUIDE v2.5 29-Aug-2017 12:35:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ThresholdGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ThresholdGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ThresholdGUI is made visible.
function ThresholdGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ThresholdGUI (see VARARGIN)

set(handles.ThresholdGUI,'units','normalized','position',[0 1 0.95 0.7]);

im = evalin('base','orig');
axes(handles.OrigImg);
imshow(im);   

% Choose default command line output for ThresholdGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes ThresholdGUI wait for user response (see UIRESUME)
% uiwait(handles.ThresholdGUI);


% --- Outputs from this function are returned to the command line.
function varargout = ThresholdGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'ThreshValue') returns toggle state of checkbox1
box = get(hObject,'Value'); 
% Store application data
setappdata(handles.ThresholdGUI,'ShowImgs',box); 


% --- Executes on button press in DecidedButton.
function DecidedButton_Callback(hObject, eventdata, handles)
% hObject    handle to DecidedButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

assignin('base','ShowImg', getappdata(handles.ThresholdGUI,'ShowImgs'));
assignin('base','NoiseIm', getappdata(handles.ThresholdGUI,'NoiseIm'));
close(handles.ThresholdGUI);


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function NoiseSlider_Callback(hObject, eventdata, handles)
% hObject    handle to NoiseSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'ThreshValue') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = get(hObject,'Value');
val=round(val);

set(handles.PixelValue,'String',val);
setappdata(handles.ThresholdGUI,'noise',val); 



% --- Executes during object creation, after setting all properties.
function NoiseSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NoiseSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in TryThis.
function TryThis_Callback(hObject, eventdata, handles)
% hObject    handle to TryThis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
assignin('base','adjust', getappdata(handles.ThresholdGUI,'adjust'));
adjust = evalin('base','adjust');
zs = evalin('base','zs');
s = evalin('base','s');
img(:,:,:) = evalin('base','img');
BinaryThresholdedImg=zeros(s(1),s(2),zs);
h = waitbar(0,'Calculating...');
for i=1:zs %For each slice
    waitbar (i/zs,h);
    level=graythresh(img(:,:,i)); %Finds threshold level using Otsu's method. Different for each slice b/c no need to keep intensity conistent, just want to pick up all mg processes. Tried adaptive threshold, but did not produce better images.
    level=level*adjust; %Increase the threshold by *1.6 to get all fine processes 
    BinaryThresholdedImg(:,:,i) = im2bw(img(:,:,i),level);% Apply the adjusted threshold and convert from gray the black white image.
end
if isgraphics(h);
close(h);
end
setappdata(handles.ThresholdGUI,'BinaryThresholdedImg',BinaryThresholdedImg);
assignin('base','BinaryThresholdedImg', getappdata(handles.ThresholdGUI,'BinaryThresholdedImg'));


% --- Executes on button press in UpdateButton.
function UpdateButton_Callback(hObject, eventdata, handles)
% hObject    handle to UpdateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
assignin('base','noise', getappdata(handles.ThresholdGUI,'noise'));
BinaryIm(:,:,:) = evalin('base','BinaryThresholdedImg');
midslice = evalin('base','midslice');
noise = evalin('base','noise');

h=msgbox('One moment please...');
NoiseIm=bwareaopen(BinaryIm,noise); %Removes objects smaller than set value (in pixels). For 3D inputs, uses automatic connectivity input of 26. Don't want small background dots left over from decreased threshold.
close(h);

axes(handles.FilteredImg);
imshow(NoiseIm(:,:,midslice));

setappdata(handles.ThresholdGUI,'NoiseIm',NoiseIm);


% --- Executes on button press in HowManyObjects.
function HowManyObjects_Callback(~, eventdata, handles)
% hObject    handle to HowManyObjects (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
assignin('base','NoiseIm', getappdata(handles.ThresholdGUI,'NoiseIm'));
NoiseIm = evalin('base','NoiseIm');

h=msgbox('Counting...');
ConnectedComponents=bwconncomp(NoiseIm,26); %returns structure with 4 fields. PixelIdxList contains a 1-by-NumObjects cell array where the k-th element in the cell array is a vector containing the linear indices of the pixels in the k-th object. 26 defines connectivity. This looks at cube of connectivity around pixel.
close(h);

numObj = numel(ConnectedComponents.PixelIdxList); %PixelIdxList is field with list of pixels in each connected component. Find how many connected components there are.

set(handles.NumberOfObjects,'String',numObj);



% --- Executes on slider movement.
function ThreshSlider_Callback(hObject, eventdata, handles)
% hObject    handle to ThreshSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'ThreshValue') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

val = get(hObject,'Value');
val=round(val,1);

im = evalin('base','orig');
level=graythresh(im); %Finds threshold level using Otsu's method. Different for each slice b/c no need to keep intensity conistent, just want to pick up all mg processes. Tried adaptive threshold, but did not produce better images.
level=level*val;  
BinaryThresholdedImg = im2bw((im),level);% Apply the adjusted threshold and convert from gray the black white image.

axes(handles.ThreshImg);
imshow(BinaryThresholdedImg);
set(handles.ThreshValue,'String',val);
setappdata(handles.ThresholdGUI,'adjust',val);
setappdata(handles.ThresholdGUI,'BinaryThresholdedImg',BinaryThresholdedImg);


% --- Executes during object creation, after setting all properties.
function ThreshSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ThreshSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
