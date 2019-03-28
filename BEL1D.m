function varargout = BEL1D(varargin)
% BEL1D MATLAB code for BEL1D.fig
%      BEL1D, by itself, creates a new BEL1D or raises the existing
%      singleton*.
%
%      H = BEL1D returns the handle to a new BEL1D or the handle to
%      the existing singleton*.
%
%      BEL1D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BEL1D.M with the given input arguments.
%
%      BEL1D('Property','Value',...) creates a new BEL1D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BEL1D_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BEL1D_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BEL1D

% Last Modified by GUIDE v2.5 08-Mar-2019 11:22:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BEL1D_OpeningFcn, ...
                   'gui_OutputFcn',  @BEL1D_OutputFcn, ...
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


% --- Executes just before BEL1D is made visible.
function BEL1D_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BEL1D (see VARARGIN)
addpath([pwd '/GlobalFunctions']);
% Adding the logos of the Universtity of Liège and Ghent and FRS-FNRS:
axes(handles.axes1)
Logo = imread('LogoULiege.jpg');
Logo = double(Logo)/255;
idx1 = Logo(:,:,1) == 1;
idx2 = Logo(:,:,2) == 1;
idx3 = Logo(:,:,3) == 1;
idxWhite = idx1 + idx2 + idx3 == 3;
for idx = 1 : 3,
    rgb = Logo(:,:,idx);
    rgb(idxWhite) = 240/255;
    Logo(:,:,idx) = rgb;
end
image(Logo);
axis off
axis image

axes(handles.axes2)
Logo = imread('LogoFNRS.png');
Logo = double(Logo)/255;
idx1 = Logo(:,:,1) == 0;
idx2 = Logo(:,:,2) == 0;
idx3 = Logo(:,:,3) == 0;
idxWhite = idx1 + idx2 + idx3 == 3;
for idx = 1 : 3,
    rgb = Logo(:,:,idx);
    rgb(idxWhite) = 240/255;
    Logo(:,:,idx) = rgb;
end
image(Logo);
axis off
axis image

axes(handles.axes3)
Logo = imread('Logo_UGent.png');
Logo = double(Logo)/255;
idx1 = Logo(:,:,1) == 1;
idx2 = Logo(:,:,2) == 1;
idx3 = Logo(:,:,3) == 1;
idxWhite = idx1 + idx2 + idx3 == 3;
for idx = 1 : 3,
    rgb = Logo(:,:,idx);
    rgb(idxWhite) = 240/255;
    Logo(:,:,idx) = rgb;
end
image(Logo);
axis off
axis image
% Choose default command line output for BEL1D
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BEL1D wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BEL1D_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbuttonOthers.
function pushbuttonOthers_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonOthers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
BEL1DGENERAL;
close(BEL1D);

% --- Executes on button press in pushbuttonSW.
function pushbuttonSW_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
BEL1DSW;
close(BEL1D);

% --- Executes on button press in pushbuttonSNMR.
function pushbuttonSNMR_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSNMR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
BEL1DSNMR;
close(BEL1D);

% --- Executes on button press in pushbuttonLICENSE.
function pushbuttonLICENSE_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLICENSE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fid = fopen('LICENSE.md','r');
tline = fgetl(fid);
while ischar(tline)
    disp(tline)
    tline = fgetl(fid);
end
fclose(fid);
