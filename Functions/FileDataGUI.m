function varargout = FileDataGUI(varargin)
% FILEDATAGUI MATLAB code for FileDataGUI.fig
%      FILEDATAGUI, by itself, creates a new FILEDATAGUI or raises the existing
%      singleton*.
%
%      H = FILEDATAGUI returns the handle to a new FILEDATAGUI or the handle to
%      the existing singleton*.
%
%      FILEDATAGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FILEDATAGUI.M with the given input arguments.
%
%      FILEDATAGUI('Property','Value',...) creates a new FILEDATAGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FileDataGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FileDataGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FileDataGUI

% Last Modified by GUIDE v2.5 21-Sep-2017 13:26:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FileDataGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @FileDataGUI_OutputFcn, ...
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


% --- Executes just before FileDataGUI is made visible.
function FileDataGUI_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FileDataGUI (see VARARGIN)

% Choose default command line output for FileDataGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FileDataGUI wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FileDataGUI_OutputFcn(~, ~, ~) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

% varargout{1} = handles.output;


% --- Executes on button press in DoneButton.
function DoneButton_Callback(~, ~, handles)
% hObject    handle to DoneButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

assignin('base','ch', getappdata(handles.figure1,'ch'));
assignin('base','scale', getappdata(handles.figure1,'scale'));
assignin('base','zscale', getappdata(handles.figure1,'zscale'));
assignin('base','file', getappdata(handles.figure1,'file'));
assignin('base','pathname', getappdata(handles.figure1,'pathname'));
assignin('base','ChannelOfInterest', getappdata(handles.figure1,'ChannelOfInterest'));

close(handles.figure1);


function ChannelsInput_Callback(hObject, ~, handles)
% hObject    handle to ChannelsInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ChannelsInput as text
% str2double(get(hObject,'String')) returns contents of ChannelsInput as a double
ch = str2double(get(hObject,'String')); 
% Store application data
setappdata(handles.figure1,'ch',ch); 

% --- Executes during object creation, after setting all properties.
function ChannelsInput_CreateFcn(hObject, ~, ~)
% hObject    handle to ChannelsInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function XYScaleInput_Callback(hObject, ~, handles)
% hObject    handle to XYScaleInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of XYScaleInput as text
%        str2double(get(hObject,'String')) returns contents of XYScaleInput as a double
scale = str2double(get(hObject,'String'));  
% Store application data
setappdata(handles.figure1,'scale',scale); 

% --- Executes during object creation, after setting all properties.
function XYScaleInput_CreateFcn(hObject, ~, ~)
% hObject    handle to XYScaleInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ZScaleInput_Callback(hObject, ~, handles)
% hObject    handle to ZScaleInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ZScaleInput as text
%        str2double(get(hObject,'String')) returns contents of ZScaleInput as a double
zscale = str2double(get(hObject,'String'));  
% Store application data
setappdata(handles.figure1,'zscale',zscale); 

% --- Executes during object creation, after setting all properties.
function ZScaleInput_CreateFcn(hObject, ~, ~)
% hObject    handle to ZScaleInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% % 
% % function FilenameInput_Callback(hObject, ~, handles)
% % % hObject    handle to FilenameInput (see GCBO)
% % % eventdata  reserved - to be defined in a future version of MATLAB
% % % handles    structure with handles and user data (see GUIDATA)
% % 
% % % Hints: get(hObject,'String') returns contents of FilenameInput as text
% % %        str2double(get(hObject,'String')) returns contents of FilenameInput as a double
% % file = get(hObject,'String');  
% % % Store application data
% % setappdata(handles.figure1,'file',file); 

% --- Executes during object creation, after setting all properties.
function FilenameInput_CreateFcn(hObject, ~, ~)
% hObject    handle to FilenameInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ChannelOfInterest_Callback(hObject, eventdata, handles)
% hObject    handle to ChannelOfInterest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ChannelOfInterest as text
%        str2double(get(hObject,'String')) returns contents of ChannelOfInterest as a double

ChannelOfInterest = str2double(get(hObject,'string')); 

% Store application data
setappdata(handles.figure1,'ChannelOfInterest',ChannelOfInterest);


% --- Executes during object creation, after setting all properties.
function ChannelOfInterest_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ChannelOfInterest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% % % --- Executes on button press in DoneButton.
% % function DoneButton2_Callback(hObject, eventdata, handles)
% % % hObject    handle to DoneButton (see GCBO)
% % % eventdata  reserved - to be defined in a future version of MATLAB
% % % handles    structure with handles and user data (see GUIDATA)
% % assignin('base','ch', getappdata(handles.figure1,'ch'));
% % assignin('base','scale', getappdata(handles.figure1,'scale'));
% % assignin('base','zscale', getappdata(handles.figure1,'zscale'));
% % assignin('base','file', getappdata(handles.figure1,'file'));
% % assignin('base','ChannelOfInterest', getappdata(handles.figure1,'ChannelOfInterest'));
% % close(handles.figure1);


% --- Executes on key press with focus on DoneButton and none of its controls.
function DoneButton_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to DoneButton (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

if strcmp(eventdata.Key,'return')
    assignin('base','ch', getappdata(handles.figure1,'ch'));
    assignin('base','scale', getappdata(handles.figure1,'scale'));
    assignin('base','zscale', getappdata(handles.figure1,'zscale'));
    assignin('base','file', getappdata(handles.figure1,'file'));
    assignin('base','pathname', getappdata(handles.figure1,'pathname'));
    assignin('base','ChannelOfInterest', getappdata(handles.figure1,'ChannelOfInterest'));
    close(handles.figure1);
end


% --- Executes on button press in SelectFileButton.
function SelectFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to SelectFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile({'*.tif';'*.lsm'});

setappdata(handles.figure1,'file',FileName);
setappdata(handles.figure1,'pathname',PathName);
set(handles.ShowFileName,'String',FileName);
