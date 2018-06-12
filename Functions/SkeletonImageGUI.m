function varargout = SkeletonImageGUI(varargin)
% SKELETONIMAGEGUI MATLAB code for SkeletonImageGUI.fig
%      SKELETONIMAGEGUI, by itself, creates a new SKELETONIMAGEGUI or raises the existing
%      singleton*.
%
%      H = SKELETONIMAGEGUI returns the handle to a new SKELETONIMAGEGUI or the handle to
%      the existing singleton*.
%
%      SKELETONIMAGEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SKELETONIMAGEGUI.M with the given input arguments.
%
%      SKELETONIMAGEGUI('Property','Value',...) creates a new SKELETONIMAGEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SkeletonImageGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SkeletonImageGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SkeletonImageGUI

% Last Modified by GUIDE v2.5 17-Mar-2018 11:30:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SkeletonImageGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SkeletonImageGUI_OutputFcn, ...
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


% --- Executes just before SkeletonImageGUI is made visible.
function SkeletonImageGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SkeletonImageGUI (see VARARGIN)

% Choose default command line output for SkeletonImageGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SkeletonImageGUI wait for user response (see UIRESUME)
% uiwait(handles.SkeletonImageGUI);


% --- Outputs from this function are returned to the command line.
function varargout = SkeletonImageGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in SkeletonImage.
function SkeletonImage_Callback(hObject, eventdata, handles)
% hObject    handle to SkeletonImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SkeletonImage
SkelBox = get(hObject,'Value'); 
% Store application data
setappdata(handles.SkeletonImageGUI,'SkelImg',SkelBox); 



% --- Executes on button press in EndpointImage.
function EndpointImage_Callback(hObject, eventdata, handles)
% hObject    handle to EndpointImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EndpointImage
EndBox = get(hObject,'Value'); 
% Store application data
setappdata(handles.SkeletonImageGUI,'EndImg',EndBox); 


% --- Executes on button press in BranchPointImage.
function BranchPointImage_Callback(hObject, eventdata, handles)
% hObject    handle to BranchPointImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BranchPointImage
BranchBox = get(hObject,'Value'); 
% Store application data
setappdata(handles.SkeletonImageGUI,'BranchImg',BranchBox); 


% --- Executes on button press in ContinueButton.
function ContinueButton_Callback(hObject, eventdata, handles)
% hObject    handle to ContinueButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
assignin('base','SkelImg', getappdata(handles.SkeletonImageGUI,'SkelImg'));
assignin('base','EndImg', getappdata(handles.SkeletonImageGUI,'EndImg'));
assignin('base','BranchImg', getappdata(handles.SkeletonImageGUI,'BranchImg'));
assignin('base','OrigCellImg', getappdata(handles.SkeletonImageGUI,'OrigCellImg'));
assignin('base','BranchLengthFile', getappdata(handles.SkeletonImageGUI,'BranchLengthFile'));
close(handles.SkeletonImageGUI);


% --- Executes on button press in OrigCellImg.
function OrigCellImg_Callback(hObject, eventdata, handles)
% hObject    handle to OrigCellImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OrigCellImg
OrigCellBox = get(hObject,'Value'); 
% Store application data
setappdata(handles.SkeletonImageGUI,'OrigCellImg',OrigCellBox); 


% --- Executes on button press in BranchLengths.
function BranchLengths_Callback(hObject, eventdata, handles)
% hObject    handle to BranchLengths (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
BranchLengthBox = get(hObject,'Value'); 
% Hint: get(hObject,'Value') returns toggle state of BranchLengths
setappdata(handles.SkeletonImageGUI,'BranchLengthFile',BranchLengthBox); 
