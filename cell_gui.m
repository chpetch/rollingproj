function varargout = cell_gui(varargin)
% CELL_GUI MATLAB code for cell_gui.fig
%      CELL_GUI, by itself, creates a new CELL_GUI or raises the existing
%      singleton*.
%
%      H = CELL_GUI returns the handle to a new CELL_GUI or the handle to
%      the existing singleton*.
%
%      CELL_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CELL_GUI.M with the given input arguments.
%
%      CELL_GUI('Property','Value',...) creates a new CELL_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cell_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cell_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cell_gui

% Last Modified by GUIDE v2.5 30-Mar-2015 22:12:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cell_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @cell_gui_OutputFcn, ...
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


% --- Executes just before cell_gui is made visible.
function cell_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cell_gui (see VARARGIN)

% Choose default command line output for cell_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cell_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = cell_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in cb_bd.
function cb_bd_Callback(hObject, eventdata, handles)
% hObject    handle to cb_bd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_bd


% --- Executes on button press in pb_openimg.
function pb_openimg_Callback(hObject, eventdata, handles)
% hObject    handle to pb_openimg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filetype = '*.tif';
[FileName,PathName] = uigetfile(filetype,'Select the video file');

handles.pathname = [PathName,FileName];
handles.filename = FileName;

guidata(hObject,handles)

% --- Executes on button press in pb_data.
function pb_data_Callback(hObject, eventdata, handles)
% hObject    handle to pb_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filetype = '*.mat';
[dFileName,dPathName] = uigetfile(filetype,'Select the data file');

handles.dpathname = [dPathName,dFileName];
handles.dfilename = dFileName;

guidata(hObject,handles)

% --- Executes on button press in pb_next.
function pb_next_Callback(hObject, eventdata, handles)
% hObject    handle to pb_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pb_previous.
function pb_previous_Callback(hObject, eventdata, handles)
% hObject    handle to pb_previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pb_clfig.
function pb_clfig_Callback(hObject, eventdata, handles)
% hObject    handle to pb_clfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
