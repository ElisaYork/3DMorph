function varargout = FullCellsGUI(varargin)
% FULLCELLSGUI MATLAB code for FullCellsGUI.fig
%      FULLCELLSGUI, by itself, creates a new FULLCELLSGUI or raises the existing
%      singleton*.
%
%      H = FULLCELLSGUI returns the handle to a new FULLCELLSGUI or the handle to
%      the existing singleton*.
%
%      FULLCELLSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FULLCELLSGUI.M with the given input arguments.
%
%      FULLCELLSGUI('Property','Value',...) creates a new FULLCELLSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FullCellsGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FullCellsGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FullCellsGUI

% Last Modified by GUIDE v2.5 19-Jul-2017 16:31:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FullCellsGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @FullCellsGUI_OutputFcn, ...
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


% --- Executes just before FullCellsGUI is made visible.
function FullCellsGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FullCellsGUI (see VARARGIN)
ShowCells = evalin('base','ShowCells');
numObjSep = evalin('base','numObjSep');
Microglia = evalin('base','Microglia');
s = evalin('base','s');
zs = evalin('base','zs');
voxscale = evalin('base','voxscale');
SepObjectList = evalin('base','SepObjectList');
udSepObjectList = evalin('base','udSepObjectList');%ObjectList is large to small, flip upside down so small is plotted first in blue.

num = [1:3:(3*numObjSep+1)];
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
        
cmap = jet(max(fullimg(:)));
cmap(1,:) = zeros(1,3);
axes(handles.AllObjects);
imagesc(fullimg);
colormap(cmap);
colorbar('Ticks',[1,3*numObjSep+1], 'TickLabels',{'Small','Large'});
set(handles.AllObjects,'visible', 'off');

%Convert pixel sizes to unit^3
SepObjectListInUnit = round(SepObjectList(:,1)*voxscale);
set(handles.listbox1,'String',SepObjectListInUnit(:,1));
setappdata(handles.FullCellsGUI,'cmap',cmap);
assignin('base','cmap', getappdata(handles.FullCellsGUI,'cmap'));
setappdata(handles.FullCellsGUI,'SepObjectList',SepObjectList);
assignin('base','SepObjectList', getappdata(handles.FullCellsGUI,'SepObjectList'));

% Choose default command line output for FullCellsGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FullCellsGUI wait for user response (see UIRESUME)
% uiwait(handles.FullCellsGUI);


% --- Outputs from this function are returned to the command line.
function varargout = FullCellsGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Done.
function Done_Callback(hObject, eventdata, handles)
% hObject    handle to Done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

assignin('base','SmCellCutoff', getappdata(handles.FullCellsGUI,'SmCellCutoff'));
assignin('base','KeepAllCells', getappdata(handles.FullCellsGUI,'KeepAllCells'));
assignin('base','RemoveXY', getappdata(handles.FullCellsGUI,'RemoveXY'));
close(handles.FullCellsGUI);


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
% numObjSep = evalin('base','numObjSep');
cmap = evalin('base','cmap');
SepObjectList = evalin('base','SepObjectList');
AllSeparatedObjs = evalin('base','AllSeparatedObjs');
SelectedValue = get(handles.listbox1,'Value');

Object = SepObjectList(SelectedValue,2);
axes(handles.SelectedObject);
imshow(AllSeparatedObjs(:,:,Object));
colormap(cmap);
setappdata(handles.FullCellsGUI,'SmCellCutoff',SepObjectList(SelectedValue,1));


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



% --- Executes on button press in RemoveXYcells.
function RemoveXYcells_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveXYcells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RemoveXY = get(hObject,'Value');
setappdata(handles.FullCellsGUI,'RemoveXY',RemoveXY);

% Hint: get(hObject,'Value') returns toggle state of RemoveXYcells


% --- Executes on button press in KeepAllObjects.
function KeepAllObjects_Callback(hObject, eventdata, handles)
% hObject    handle to KeepAllObjects (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

KeepAllCells = get(hObject,'Value');
setappdata(handles.FullCellsGUI,'KeepAllCells',KeepAllCells);

% Hint: get(hObject,'Value') returns toggle state of KeepAllObjects
