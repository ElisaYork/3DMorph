function varargout = SelectFilesGUI(varargin)
% SELECTFILESGUI MATLAB code for SelectFilesGUI.fig
%      SELECTFILESGUI, by itself, creates a new SELECTFILESGUI or raises the existing
%      singleton*.
%
%      H = SELECTFILESGUI returns the handle to a new SELECTFILESGUI or the handle to
%      the existing singleton*.
%
%      SELECTFILESGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECTFILESGUI.M with the given input arguments.
%
%      SELECTFILESGUI('Property','Value',...) creates a new SELECTFILESGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SelectFilesGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SelectFilesGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SelectFilesGUI

% Last Modified by GUIDE v2.5 29-Aug-2017 17:09:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SelectFilesGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SelectFilesGUI_OutputFcn, ...
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


% --- Executes just before SelectFilesGUI is made visible.
function SelectFilesGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SelectFilesGUI (see VARARGIN)

% Choose default command line output for SelectFilesGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SelectFilesGUI wait for user response (see UIRESUME)
% uiwait(handles.SelectFilesGUI);


% --- Outputs from this function are returned to the command line.
function varargout = SelectFilesGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in SettingsButton.
function SettingsButton_Callback(hObject, eventdata, handles)
% hObject    handle to SettingsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Parameters = uigetfile('*.mat');
setappdata(handles.SelectFilesGUI,'Parameters',Parameters);


% --- Executes on button press in FilesButton.
function FilesButton_Callback(hObject, eventdata, handles)
% hObject    handle to FilesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileList,PathList] = uigetfile('*.*','MultiSelect','on');
setappdata(handles.SelectFilesGUI,'FileList',FileList);
setappdata(handles.SelectFilesGUI,'PathList',PathList);


% --- Executes on button press in GoButton.
function GoButton_Callback(hObject, eventdata, handles)
% hObject    handle to GoButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

assignin('base','Parameters', getappdata(handles.SelectFilesGUI,'Parameters'));
assignin('base','FileList', getappdata(handles.SelectFilesGUI,'FileList'));
assignin('base','PathList', getappdata(handles.SelectFilesGUI,'PathList'));
close(handles.SelectFilesGUI);
