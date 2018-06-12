function varargout = CellSizeCutoffGUI(varargin)
% CELLSIZECUTOFFGUI MATLAB code for CellSizeCutoffGUI.fig
%      CELLSIZECUTOFFGUI, by itself, creates a new CELLSIZECUTOFFGUI or raises the existing
%      singleton*.
%
%      H = CELLSIZECUTOFFGUI returns the handle to a new CELLSIZECUTOFFGUI or the handle to
%      the existing singleton*.
%
%      CELLSIZECUTOFFGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CELLSIZECUTOFFGUI.M with the given input arguments.
%
%      CELLSIZECUTOFFGUI('Property','Value',...) creates a new CELLSIZECUTOFFGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CellSizeCutoffGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CellSizeCutoffGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CellSizeCutoffGUI

% Last Modified by GUIDE v2.5 23-Nov-2017 11:06:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CellSizeCutoffGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @CellSizeCutoffGUI_OutputFcn, ...
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


% --- Executes just before CellSizeCutoffGUI is made visible.
function CellSizeCutoffGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CellSizeCutoffGUI (see VARARGIN)


numObj = evalin('base','numObj');
ConnectedComponents = evalin('base','ConnectedComponents');
s = evalin('base','s');
zs = evalin('base','zs');
voxscale = evalin('base','voxscale');
ObjectList = evalin('base','ObjectList');
udObjectList = flipud(ObjectList);%ObjectList is large to small, flip upside down so small is plotted first in blue.

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
        
cmap = jet(max(fullimg(:)));
cmap(1,:) = zeros(1,3);
axes(handles.AllObjects);
imagesc(fullimg);
colormap(cmap);
colorbar('Ticks',[1,3*numObj+1], 'TickLabels',{'Small','Large'});
set(handles.AllObjects,'visible', 'off');

%Convert pixel sizes to unit^3
ObjectListInUnit = round(ObjectList(:,1)*voxscale);
set(handles.listbox1,'String',ObjectListInUnit(:,1));
setappdata(handles.CellSizeCutoffGUI,'cmap',cmap);
assignin('base','cmap', getappdata(handles.CellSizeCutoffGUI,'cmap'));

setappdata(handles.CellSizeCutoffGUI,'fullimg',fullimg); 
% Choose default command line output for CellSizeCutoffGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CellSizeCutoffGUI wait for user response (see UIRESUME)
% uiwait(handles.CellSizeCutoffGUI);


% --- Outputs from this function are returned to the command line.
function varargout = CellSizeCutoffGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% --- Executes on button press in SelectObjectSize.
function SelectObjectSize_Callback(hObject, eventdata, handles)
% hObject    handle to SelectObjectSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

assignin('base','CellSizeCutoff', getappdata(handles.CellSizeCutoffGUI,'CellSizeCutoff'));
assignin('base','fullimg',getappdata(handles.CellSizeCutoffGUI,'fullimg'));
assignin('base','ShowObjImg', getappdata(handles.CellSizeCutoffGUI,'ShowObjImgs'));

close(handles.CellSizeCutoffGUI);


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

allObjs = evalin('base','allObjs');
ObjectList = evalin('base','ObjectList');
cmap = evalin('base','cmap');

SelectedValue = get(handles.listbox1,'Value');
Object = ObjectList(SelectedValue,2);
axes(handles.SelectedObject);
imshow(allObjs(:,:,Object));
colormap(cmap);
setappdata(handles.CellSizeCutoffGUI,'CellSizeCutoff',ObjectList(SelectedValue,1));



% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in TellMeMoreList.
function TellMeMoreList_Callback(hObject, eventdata, handles)
% hObject    handle to TellMeMoreList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns TellMeMoreList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TellMeMoreList


% --- Executes during object creation, after setting all properties.
function TellMeMoreList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TellMeMoreList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in OutputImage.
function OutputImage_Callback(hObject, eventdata, handles)
% hObject    handle to OutputImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OutputImage

box = get(hObject,'Value'); 
% Store application data
setappdata(handles.CellSizeCutoffGUI,'ShowObjImgs',box); 
