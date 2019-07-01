function varargout = BEL1DSNMR(varargin)
% BEL1DSNMR is a graphical interface that helps the user to perform
% prediction-focused approach imaging of the subsurface using SNMR data.
% This code relies heavily on MRSmatlab codes (Mueller-Petke et al., 2016),
% and therefore REQUIRES the presence of all MRSmatlab functions in the
% MATLAB path (run the shortcut previously to any uses of the interface).
%
% In order to run properly:
%       - The functions of MRSmatlab (Mueller-Petke et al., 2016) are
%         available to MATLAB.
%       - The files formatting needs to correspond to the formatting of
%         MRSmatlab (Mueller-Petke et al., 2016).
%
% The different steps to use this interface are described in the User
% Manual provided with this code.
%
% This GUI is adapted for the use of multiple-loops datasets (Kremer et
% al., 2018). The multiple-loops datasets are suspected to improve accuracy
% on the results as shown by sensitivity analysis. However, the GUI work on
% simple configurations as well, given that either the loops are circular, 
% or that the user provides an accurate kernel in the MRSmatlab formatting.
%
% The GUI is organized in multiple panels that appear progressively and 
% enables the PFA imaging of the dataset in an intuitive way.
% 
% The panels are:
%
%   1) Data panel: This panel enables the user to load the data. The user
%                  is first required to select a directory where the 
%                  different files (the *.mrsd files with the data and 
%                  (optional) the *.mrsk files with the kernels) are 
%                  located. Then, the user must complete the table with the
%                  names (with the *.mrsd extension) of the data files 
%                  (e.g.: datatx50rx50.mrsd). Once the name is correctly
%                  written, the user can activate the file by checking the
%                  corresponding 'active' box.
%                  Once the file is loaded, the user can select a
%                  particular configuration (particularly suited for
%                  multiple-loops configurations. When the configuration is
%                  selected, click the 'Confirm configuration' pushbutton.
%                  This will lock the configuration (no coming back) an
%                  enable further steps.
%
%   2) The prior model space panel: This panel enables the fine tunning of
%                                   the prior model space. The user is
%                                   required to enter the distributions
%                                   parameters and to choose for a given
%                                   distribution type/sampling method. When
%                                   clicking the 'Generate models'
%                                   pushbutton, the models are sampled
%                                   accordingly. Once again, there is no
%                                   coming back.
%
%   3) The kernels panel: This panel enables to either select a previously
%                         computed kernel file, or to compute the kernel if
%                         the file does not exist. To load a computed file,
%                         enter its name (with the *.mrsk extension) in the
%                         file name  colum and select 'computed'. If the
%                         file is not yet computed, the user can compute
%                         the kernel through a call to MRSmatlab. Be aware
%                         that the only possibilities through this
%                         interface are to compute kernels under a
%                         resistive earth hypothesis or a 100 Ohm.m
%                         resistivity. For more options, the user is
%                         reffered to the MRSKernel GUI from MRSmatlab.
%   
%   4) The PFA imaging panel: This panel enables the computation of the
%                             posterior distributions using the PFA imaging
%                             process. To launch the computations, the user
%                             must click on the 'Run PFA imaging' button. A
%                             progress bar appears and shows at which step
%                             the computations are. When the computations
%                             are finished, the user can analyze the
%                             results through the parameters distributions
%                             or the model distributions. The user can also
%                             save the obtained results in a *.mrspfa file
%                             for further computations. The saved file will
%                             contain the simulated data for the posterior
%                             space only if the model distributions have
%                             been displayed previously.
%
% For further help, please contact the author:
%       Name: Hadrien MICHEL
%       University: University of Liège
%       Country: Belgium
%       e-mail: Hadrien.Michel@uliege.be
%
% © June 2018 Hadrien MICHEL, University of Liège and University of Ghent
%
% Modifications of version 2 (Septembre/October 2018):
%       - decoupling of datasets
%       - improvement of results visualisation

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BEL1DSNMR_OpeningFcn, ...
                   'gui_OutputFcn',  @BEL1DSNMR_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before BEL1DSNMR is made visible.
function BEL1DSNMR_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)
addpath([pwd '/GlobalFunctions']);

% Choose default command line output for BEL1DSNMR

handles.output = hObject;

% Custom variables that are contained in the GUI:
handles.Files = [];% Contains the files informations
handles.Files.path = 'Files path';
handles.Files.data.namesFiles = {'Name 1','Name 2','Name 3','Name 4'};
handles.Files.data.active = [false; false; false; false];
handles.Files.kernel.names = {'Name 1','Name 2','Name 3','Name 4', 'Name 5'...
    ,'Name 6','Name 7','Name 8', 'Name 9', 'Name 10', 'Name 11', 'Name 12',...
    'Name 13', 'Name 14', 'Name 15', ',Name 16'};
handles.Files.kernel.namesFiles = {'Name 1', 'Name 2', 'Name 3', 'Name 4',...
    'Name 5', 'Name 6', 'Name 7', 'Name 8', 'Name 9', 'Name 10', 'Name 11',...
    'Name 12', 'Name 13', 'Name 14', 'Name 15', ',Name 16'};
handles.Files.kernel.exist = false(1,16);% By default, all the kernels needs to be computed
handles.Data = [];% Contains the data-sets
handles.Data.proclog1 = [];
handles.Data.proclog2 = [];
handles.Data.proclog3 = [];
handles.Data.proclog4 = [];
handles.Data.used = [];
handles.Data.used.Rx = zeros(4,4);%Sizes of the receivers: handles.Data.used.Rx(i,j) = size of the jth receiver of the ith dataset (Tx) 
handles.Data.activated = false(4,4);%handles.Data.activated(i,j) = true --> The ith transmitter with the jth receiver is used!
handles.Data.used.Tx = [0 0 0 0];%Sizes of the receivers: handles.Data.used.Tx(i,j) = size of the transmitter of the ith dataset (Tx) 
handles.Data.used.config.OK = false(1,16);
handles.Data.proclogNamesCode = {'handles.Data.proclog1','handles.Data.proclog2','handles.Data.proclog3','handles.Data.proclog4'};
handles.Kernel = [];% Contains the kernel fucntions
handles.Models = [];% Contains the models and their characteristics 
handles.Models.param = [2.5 7.5 3.5 10 5 350;0 5 10 30 5 350;0 5 3.5 35 5 350;...
    0 5 3.5 35 5 350;0 5 3.5 35 5 350;0 5 3.5 35 5 350]; % Default prior model space
handles.Models.nbLayers = 2;% Default number of layers
set(handles.uitable_prior,'BackgroundColor',[repmat([1 1 1],handles.Models.nbLayers,1);repmat([0 0 0],6-handles.Models.nbLayers,1)]);
handles.Models.type = 1; % Absolute value data representation
handles.Models.N = 1000;
handles.Solution = [];

handles.uitable_data_names.Data(:,1) = handles.Files.data.namesFiles;

% Graphical tuning of the interface
i = 1;
while i <= 4,
    eval(['set(handles.CB_Rx' num2str(i) ',''Visible'',''off'');']);
    i = i + 1;
end
set(handles.listboxTx,'String',['unused     ' ; 'unused     ' ; 'unused     ' ; 'unused     ']);

set(handles.uipanel_model,'Visible','off');
set(handles.uipanel_kernel,'Visible','off');
set(handles.uipanel_PFA,'Visible','off');

set(handles.uibuttongroup_config,'Visible','off');
set(handles.pushbutton_Config,'Visible','off');

% Adding the logos of the Universtity of Liège and Ghent:
axes(handles.axes2)
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

% Status updates variables
handles.Status = [];
handles.Status.data = false;%data loaded?
handles.Status.config = false;%Configuration choosen?
handles.Status.models = false;%Model generated?
handles.Status.kernels = false;%Kernel computed?
handles.Status.PFA = false;%PFA performed?

% For reproductibility: controle random number generation:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                     /\                                  %
%                                    /  \                                 %
%                                   /  | \                                %
%                                  /   |  \                               %
%                                 /    .   \                              %
%                                /__________\                             %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%/!\ Comment to perform true random generation /!\%%%%%%%%%%%%%
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%%%%%%%%%%%%%/!\ Comment to perform true random generation /!\%%%%%%%%%%%%%
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = BEL1DSNMR_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.output;


% --- Executes when entered data in editable cell(s) in uitable_data_names.
function uitable_data_names_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable_data_names (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

% This function enables to get the names of the files for the PFA imaging. 
% Due to the adaptation to multiple-loops datasets, this function may
% appear quite complicated. However, the function performs:
%   - Registering the name of the file in the File structure
%   - Getting parameters of the different loops and registering those
%     caracteristics in the Data structure
%   - Producing the default names for the future panels

if eventdata.Indices(2) == 1,% Adding name of file
    handles.Files.data.namesFiles{eventdata.Indices(1)} = eventdata.EditData;
else % Enabling file
    handles.Files.data.active(eventdata.Indices(1)) = eventdata.NewData;% Activating the file
    if handles.Files.data.active(eventdata.Indices(1)),
        proclog = load([handles.Files.path '\' handles.Files.data.namesFiles{eventdata.Indices(1)}],'-mat');% Loading the activated file to gather the experimental configuration
        switch eventdata.Indices(1),% Switching between the possible 4 loaded files
            case 1
                handles.Data.proclog1 = proclog.proclog;% Assigning the data to the Data structure
                handles.Data.used.Tx(1) = handles.Data.proclog1.txinfo.loopsize;% Reading the transmitter diameter
                % Adding in the configuration menu
                tx_supp = ['Tx = ' num2str(handles.Data.used.Tx(1)) 'm'];
                while length(tx_supp) < 11,
                    tx_supp = [tx_supp ' '];
                end
                current_tx = get(handles.listboxTx,'String');
                current_tx(1,:) = tx_supp;
                set(handles.listboxTx,'String',current_tx);
                % Reading the receiver diameter(s)
                handles.Data.used.Rx(1,1:length([handles.Data.proclog1.rxinfo(:).loopsize])) = [handles.Data.proclog1.rxinfo(:).loopsize];
            case 2
                handles.Data.proclog2 = proclog.proclog;% Assigning the data to the Data structure
                handles.Data.used.Tx(2) = handles.Data.proclog2.txinfo.loopsize(end);% Reading the transmitter diameter
                % Adding in the configuration menu
                tx_supp = ['Tx = ' num2str(handles.Data.used.Tx(2)) 'm'];
                while length(tx_supp) < 11,
                    tx_supp = [tx_supp ' '];
                end
                current_tx = get(handles.listboxTx,'String');
                current_tx(2,:) = tx_supp;
                set(handles.listboxTx,'String',current_tx);
                % Reading the receiver diameter(s)
                handles.Data.used.Rx(2,1:length([handles.Data.proclog2.rxinfo(:).loopsize])) = [handles.Data.proclog2.rxinfo(:).loopsize];
            case 3
                handles.Data.proclog3 = proclog.proclog;% Assigning the data to the Data structure
                handles.Data.used.Tx(3) = handles.Data.proclog3.txinfo.loopsize(end);% Reading the transmitter diameter
                % Adding in the configuration menu
                tx_supp = ['Tx = ' num2str(handles.Data.used.Tx(3)) 'm'];
                while length(tx_supp) < 11,
                    tx_supp = [tx_supp ' '];
                end
                current_tx = get(handles.listboxTx,'String');
                current_tx(3,:) = tx_supp;
                set(handles.listboxTx,'String',current_tx);
                % Reading the receiver diameter(s)
                handles.Data.used.Rx(3,1:length([handles.Data.proclog3.rxinfo(:).loopsize])) = [handles.Data.proclog3.rxinfo(:).loopsize];
            case 4
                handles.Data.proclog4 = proclog.proclog;% Assigning the data to the Data structure
                handles.Data.used.Tx(4) = handles.Data.proclog4.txinfo.loopsize(end);% Reading the transmitter diameter
                % Adding in the configuration menu
                tx_supp = ['Tx = ' num2str(handles.Data.used.Tx(4)) 'm'];
                while length(tx_supp) < 11,
                    tx_supp = [tx_supp ' '];
                end
                current_tx = get(handles.listboxTx,'String');
                current_tx(4,:) = tx_supp;
                set(handles.listboxTx,'String',current_tx);
                % Reading the receiver diameter(s)
                handles.Data.used.Rx(4,1:length([handles.Data.proclog4.rxinfo(:).loopsize])) = [handles.Data.proclog4.rxinfo(:).loopsize];
        end
    else % Disabling files
        switch eventdata.Indices(1),
            case 1
                handles.Data.proclog1 = [];% Deleting the data from memory
                % Re-initializing the flags and parameters
                handles.Data.used.Rx(1,:) = [0 0 0 0];
                handles.Data.activated(1,:) = false;
                current_tx = get(handles.listboxTx,'String');
                current_tx(1,:) = 'unused     ';
                set(handles.listboxTx,'String',current_tx);
            case 2
                handles.Data.proclog2 = [];% Deleting the data from memory
                % Re-initializing the flags and parameters
                handles.Data.used.Rx(2,:) = [0 0 0 0];
                handles.Data.activated(2,:) = false;
                current_tx = get(handles.listboxTx,'String');
                current_tx(2,:) = 'unused     ';
                set(handles.listboxTx,'String',current_tx);
            case 3
                handles.Data.proclog3 = [];% Deleting the data from memory
                % Re-initializing the flags and parameters
                handles.Data.used.Rx(3,:) = [0 0 0 0];
                handles.Data.activated(3,:) = false;
                current_tx = get(handles.listboxTx,'String');
                current_tx(3,:) = 'unused     ';
                set(handles.listboxTx,'String',current_tx);
            case 4
                handles.Data.proclog4 = [];% Deleting the data from memory
                % Re-initializing the flags and parameters
                handles.Data.used.Rx(4,:) = [0 0 0 0];
                handles.Data.activated(4,:) = false;
                current_tx = get(handles.listboxTx,'String');
                current_tx(4,:) = 'unused     ';
                set(handles.listboxTx,'String',current_tx);
        end
    end
end

if not(isempty(handles.Files.data.active==true)),% If at least one file is activated: show the configuration panel
    set(handles.uibuttongroup_config,'Visible','on');
    set(handles.pushbutton_Config,'Visible','on');
    % Activating the listbox and the associated checkbox:
    item = 1;
    i = 1;
    % For the current transmitter, show the possible receivers
    while handles.Data.used.Rx(item,i) > 0,
        eval(['set(handles.CB_Rx' num2str(i) ',''String'',[''Rx = '' num2str(handles.Data.used.Rx(item,i)) '' m'']);']);
        eval(['set(handles.CB_Rx' num2str(i) ',''Visible'',''on'');']);
        i = i + 1;
    end
else% If no file is activated: hide the configuration panel
    set(handles.uibuttongroup_config,'Visible','off');
    set(handles.pushbutton_Config,'Visible','off');
end 
handles.Satuts.data = true;
guidata(hObject, handles);
    


function edit_data_path_Callback(hObject, eventdata, handles)
% hObject    handle to edit_data_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Used if the user manually enters the path to the files in the file path
% textbox
handles.Files.path = get(hObject,'String');
guidata(hObject, handles);


% --- Executes when selected cell(s) is changed in uitable_data_names.
function uitable_data_names_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable_data_names (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

% UNUSED --> For further developments (for example, using a gui to select
% the files)


% --- Executes on button press in pushbutton_Config.
function pushbutton_Config_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Config (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This function locks the configuration of the data space. It register 
% which datasets are used, with which transmitter/receiver couples and
% creates the default namles for the kernel files. It also searches for
% those default files in the folder reffered in the first step.

used_final = handles.Data.activated;% matrix containing the list of datasets used (diff. config)

j = 1;
k = 0;
l = 1;
handles.Data.used.config.OK = false(1,16);

% Initializing the names of the default kernel files, variables and headers
for i = 1 : 16,
    k = k + 1;
    if used_final(j,k),    
        handles.Data.used.config.OK(i) = true;
        NamesUsed{l} = ['Tx' num2str(handles.Data.used.Tx(j)) '/Rx' num2str(handles.Data.used.Rx(j,k))];
        KFiles_default{l} = ['Kernel_Tx' num2str(handles.Data.used.Tx(j)) '_Rx' num2str(handles.Data.used.Rx(j,k)) '.mrsk'];
        names_kdata{l} = ['kdata' num2str(j) num2str(k)];
        channeltx(l) = j;
        channelrx(l) = k;
        l = l + 1;
    end
    if mod(i,4) == 0,
        j = j + 1;
        k = 0;
    end
end

% Tunning the GUI for the particular configuration:
handles.Data.used.config.KNamesCode = names_kdata;
handles.Data.used.config.channelRx = channelrx;
handles.Data.used.config.channelTx = channeltx;
n = nnz(handles.Data.used.config.OK);
handles.Kernel.nb = n;
set(handles.uitable_kernels,'BackgroundColor',[repmat([1 1 1],handles.Kernel.nb,1);repmat([0 0 0],16-handles.Kernel.nb,1)]);
rowNames = get(handles.uitable_kernels,'RowName');
rowNames(1:n) = NamesUsed;
set(handles.uitable_kernels,'RowName',rowNames);
set(handles.uipanel_model,'Visible','on');

% Checking if the default kernel file already exists:
handles.Files.kernels.default_names = KFiles_default;
Data_current = get(handles.uitable_kernels,'Data');
Data_current(1:length(KFiles_default),1) = KFiles_default;
files = dir(handles.Files.path);
files = {files(:).name};
myindices = ~cellfun(@isempty,regexp(files,'.mrsk$'));
if ~isempty(find(myindices==1,1)),
    for i = 1 : length(KFiles_default),
        indices = find(myindices==1);
        for j = 1 : length(indices),
            tf = strcmp(KFiles_default{i},files{indices(j)});
            if tf == 1,
                Data_current{i,2}=true;
            end
        end
    end
end
set(handles.uitable_kernels,'Data',Data_current);
set(handles.uipanel_kernel,'Visible','on');
handles.uitable_prior.Data = handles.Models.param;
handles.Status.config = true;
guidata(hObject, handles);


% --- Executes on button press in pushbutton_addLayer.
function pushbutton_addLayer_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_addLayer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Adding layer to the model space
if handles.Models.nbLayers < 6,
    handles.Models.nbLayers = handles.Models.nbLayers + 1;
    set(handles.uitable_prior,'BackgroundColor',[repmat([1 1 1],handles.Models.nbLayers,1);repmat([0 0 0],6-handles.Models.nbLayers,1)]);
    guidata(hObject, handles);
else
    error('Too many layers in the model (max. is 6)!');
end


% --- Executes on button press in pushbutton_removeLayer.
function pushbutton_removeLayer_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_removeLayer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Removing layers from the model space
if handles.Models.nbLayers > 2,
    handles.Models.nbLayers = handles.Models.nbLayers - 1;
    set(handles.uitable_prior,'BackgroundColor',[repmat([1 1 1],handles.Models.nbLayers,1);repmat([0 0 0],6-handles.Models.nbLayers,1)]);
    guidata(hObject, handles);
else
    error('Too few layers in the model (min. is 2)!');
end


% --- Executes on selection change in popupmenu_TypeMod.
function popupmenu_TypeMod_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_TypeMod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Update of the type of sampling/sampler
handles.Models.type = get(handles.popupmenu_TypeMod, 'value');
guidata(hObject, handles);



function edit_nbLayer_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nbLayer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Updatting the number of models in the prior sample
handles.Models.N = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes on button press in pushbutton_generateModels.
function pushbutton_generateModels_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_generateModels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
addpath([pwd '/SNMR']);
% Function for model generation:
% [models] = ModelGenerator(type, N, parameters, nb_layer)
parameters = get(handles.uitable_prior,'Data');
parameters = parameters(1:handles.Models.nbLayers,:);
parameters(:,3:4) = parameters(:,3:4)./100;% Changes the values in percentage to dimensionless variables.
handles.Models.model = ModelGenerator(handles.Models.type, handles.Models.N, parameters, handles.Models.nbLayers);
handles.Status.models = true;
if handles.Status.kernels && handles.Status.models,
    set(handles.uipanel_PFA,'Visible','on');
end
guidata(hObject, handles);


% --- Executes on button press in pushbutton_loaKernels.
function pushbutton_loaKernels_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_loaKernels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Loading existing kernels:
for i = 1 : length(find(handles.Data.used.config.OK==true)),
    if handles.uitable_kernels.Data{i,2} == true,
        kdata = load([handles.Files.path '\' handles.uitable_kernels.Data{i,1}],'-mat');
        eval(['handles.Kernels.' handles.Data.used.config.KNamesCode{i} ' = kdata.kdata;']);
    end
end
if length(find([handles.uitable_kernels.Data{:,2}] == true)) == length(find(handles.Data.used.config.OK==true)),
    handles.Status.kernels = true;
end
if handles.Status.kernels && handles.Status.models,
    set(handles.uipanel_PFA,'Visible','on');
end
guidata(hObject, handles);


% --- Executes on button press in pushbutton_computeKernels.
function pushbutton_computeKernels_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_computeKernels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This is the only function that needs the presence of MRSmatlab in the
% accessible functions. It computes the missing kernel files and saves them
% in the default file in the selected folder.
addpath([pwd '/SNMR']);
Data_current = get(handles.uitable_kernels,'Data');
for i = 1 : length(find(handles.Data.used.config.OK==true)),
    if handles.uitable_kernels.Data{i,2} == false,
        kdata = ComputeKernel(eval(handles.Data.proclogNamesCode{handles.Data.used.config.channelTx(i)}), handles.Data.used.config.channelRx(i) , handles.checkbox_resInf.Value);
        eval(['handles.Kernels.' handles.Data.used.config.KNamesCode{i} ' = kdata;']);
        Data_current{i,2} = true;
        save([handles.Files.path '\' handles.uitable_kernels.Data{i,1}], 'kdata');
    end
end
set(handles.uitable_kernels,'Data',Data_current);
if length(find([handles.uitable_kernels.Data{:,2}] == true)) == length(find(handles.Data.used.config.OK==true)),
    handles.Status.kernels = true;
end
if handles.Status.kernels && handles.Status.models,
    set(handles.uipanel_PFA,'Visible','on');
end
guidata(hObject, handles);


% --- Executes on button press in checkbox_resInf.
function checkbox_resInf_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_resInf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_resInf
if get(hObject,'Value') == false,
    set(handles.text_resEarth,'Visible','on');
else
    set(handles.text_resEarth,'Visible','off');
end


% --- Executes on button press in pushbutton_Run.
function pushbutton_Run_Callback(hObject, eventdata, handles)
% This function performs the main object of the PFA imaging. It is
% decomposed in several steps:
%   1) Forward modelling
%       This steps consits on the computation of the simulated data for
%       each sampled model from the prior model space. It uses the function
%       'ForwardModelling_FID.m' (a modified version of the MRSmatlab
%       'ForwardModelling.m' function : Müller-Petke et al., 2016).
%   2) Arranging the data for further computations
%       The results obtained from the previous step are arranged into
%       another form, according to the type of data that is used (by
%       default, the absolute value is used).
%   3) Reducing dimensionality
%       a\ PCA
%       The Principle Component Annalysis (PCA) is applied to the different
%       dimensions of the problem. It enables a reduction of the
%       dimensionality while keeping at least 99% of the variability in the
%       dataset. 
%       The dimensions on the model parameters is notr reduced, as it is
%       unnecessary with uncorrelated variables.
%       b\ CCA
%       Establishing cannonical correlation between the models parameters
%       and the data. Main stage of the PFA imaging process.
%       Two figures are presented in this step: 
%           - The cannoncially correlated space
%           - The linear combinations visualisation
%   4) Kernel density estimation
%       This transforms the set of points in the CCA space into
%       probabilistic distributions that can be sampled for the next step.
%   5) Sampling and back-transformation (rejection of samples outside the
%      limits of the model space).
%
% All the lines with waitbar2a (Hatton, R. (2015). Waitbar2a(x,whichbar,
% varargin): Progress bar that can be embedded into a GUI uipanel element 
% and/or decremented (https://nl.mathworks.com/matlabcentral/fileexchange/
% 23169-waitbar2a-x-whichbar-varargin), MATLAB Central File Exchange) are 
% used to update the progress bar displayed in the GUI while the PFA
% imaging is running.
%
% The forward modelling is the most time-demanding step. It is possible to
% run this step with the parrallel computing toolbox to fasten
% computations. To do so, please, start the parallel pool *BEFORE* clicking
% on 'Run PFA Imaging' pushbutton.
%
addpath([pwd '/GlobalFunctions']);
addpath([pwd '/SNMR']);

if handles.Status.PFA ~= true,% Only possible to run if not previously run
    set(handles.uipanel_waitbar,'Visible','on');
    drawnow;
    
    % Adding the possibility to iterate (not on GUI)
    % The ability to iterate enables to reduce in a better way the
    % uncertainty on the models parameters. However, using already reduced
    % distributions may lead to errors in the posterior distributions,
    % especially if error is present in the data.
    iterate = true;% true = iterate, false = no iteration
    if iterate,
        iteration = 3; % Number of iterations
    else, 
        iteration = 1; 
    end
    iteration_init = iteration;
    passed = 0;
    % End iteration parameters
    while iteration >= 1,
        %% 1) Forward modelling:
        % Iteration initialization:
        if passed ~= 1,
            models = handles.Models.model;
        else
            if (iteration_init-iteration)==1,
                handles.Models.model_save = handles.Models.model;
            end
            handles.Models.model = handles.Solution.model;
            models = handles.Models.model;
        end 
        i = 1;
        type_mod = 1;% 1 = Absolute value, 2 = Complexe value, 3 = Real only, any other = Absolute value
        waitbar2a(0,handles.uipanel_waitbar,'waitbartext',sprintf('Please wait . . . \n Computing FID'));
        while i <= length(find(handles.Data.used.config.OK==true)),% For multiple loops datasets
            eval(['kdata = handles.Kernels.' handles.Data.used.config.KNamesCode{i} ';']);% Charging the corresponding kernel
            eval(['proclog = ' handles.Data.proclogNamesCode{handles.Data.used.config.channelTx(i)} ';']);% Charging the corresponding data
            proclog = proclog;
            kdata = kdata;
            t_register = [proclog.Q(1).rx(1).sig(2).t];% Timing of datapoints
            fid = zeros(handles.Models.N,size(kdata.K,1),length(t_register + proclog.Q(1).timing.tau_dead1));% Initialization of the matrix containing the data
            if type_mod == 2,
                fid = complex(fid);% If the values to store are complex, declare the matrix as complex.
            end
            waitbar2a(0,handles.uipanel_waitbar,sprintf('Please wait . . . \n Computing FID for model %d',(i)));
            % Forward modelling itself:
            if isempty(gcp('nocreate')),% If no parallel pool are existing, run on a simple loop
                for j = 1 : handles.Models.N,
                    if mod(j,50)==0,
                        waitbar2a(j/handles.Models.N,handles.uipanel_waitbar);
                    end
                    % Initializing the parameters:
                    mdata    = struct();
                    mdata.mod.nTau = 2;
                    mdata.mod.tau  = [0.1 0.3];

                    mdata.mod.Nlayer = handles.Models.nbLayers;
                    if models.depth(j,end) < kdata.model.zmax
                        mdata.mod.zlayer = [models.depth(j,:) kdata.model.zmax];
                    else
                        mdata.mod.zlayer = [models.depth(j,1:end-1) kdata.model.zmax-1 kdata.model.zmax];
                    end
                    mdata.mod.f      = models.water(j,:);
                    mdata.mod.T2s    = models.T2(j,:)./1e3;
                    mdata.mod.T2     = zeros(j, length(mdata.mod.zlayer)-1);
                    mdata.mod.T1     = zeros(j, length(mdata.mod.zlayer)-1);

                    mdata.mod.tfid1  = t_register + proclog.Q(1).timing.tau_dead1;
                    mdata.mod.noise  = 0;

                    %kdata.KT1 = [];

                    mdata = ForwardModelling_FID(mdata,kdata);% Call to the forward model
                    if type_mod == 1,%% Absolute values
                        fid(j,:,:) = abs(mdata.dat.fid1);
                    elseif type_mod == 2%% Complex values
                        fid(j,:,:) = mdata.dat.fid1;
                    elseif type_mod == 3,%% Real only
                        fid(j,:,:) = real(mdata.dat.fid1);
                    else
                        fid(j,:,:) = abs(mdata.dat.fid1);
                    end
                end
            else % If a parallel pool exists, run multiple times the parfor loop
                tested = [1 100 : 100 : handles.Models.N];
                for jklm = 1 : length(tested)-1,
                    waitbar2a(tested(jklm)/handles.Models.N,handles.uipanel_waitbar);
                    parfor j =  tested(jklm): tested(jklm+1),
                        % Initializing the parameters
                        mdata    = struct();
                        mdata.mod.nTau = 2;
                        mdata.mod.tau  = [0.1 0.3];

                        mdata.mod.Nlayer = handles.Models.nbLayers;
                        if models.depth(j,end) < kdata.model.zmax
                            mdata.mod.zlayer = [models.depth(j,:) kdata.model.zmax];
                        else
                            mdata.mod.zlayer = [models.depth(j,1:end-1) kdata.model.zmax-1 kdata.model.zmax];
                        end
                        mdata.mod.f      = models.water(j,:);
                        mdata.mod.T2s    = models.T2(j,:)./1e3;
                        mdata.mod.T2     = zeros(j, length(mdata.mod.zlayer)-1);
                        mdata.mod.T1     = zeros(j, length(mdata.mod.zlayer)-1);

                        mdata.mod.tfid1  = t_register + proclog.Q(1).timing.tau_dead1;
                        mdata.mod.noise  = 0;

                        mdata = ForwardModelling_FID(mdata,kdata);% Call to the forward model
                        if type_mod == 1,%% Absolute values
                            fid(j,:,:) = abs(mdata.dat.fid1);
                        elseif type_mod == 2%% Complex values
                            fid(j,:,:) = mdata.dat.fid1;
                        elseif type_mod == 3,%% Real only
                            fid(j,:,:) = real(mdata.dat.fid1);
                        else
                            fid(j,:,:) = abs(mdata.dat.fid1);
                        end
                    end
                end
            end
            eval(['handles.Models.model.results' num2str(i) ' = fid;']);% Assining the result to the right place
            i = i + 1;% Getting to the next receiver loop in the configuration
        end
        clearvars('kdata','proclog');% Useless in the next steps -> delete
        guidata(hObject, handles);% Updating the data in the GUI

        waitbar2a(0.5,handles.uipanel_waitbar,sprintf('Formatting the datasets . . . '));

        %% 2) Grouping real and simulated data:
        % 2.a) Simulated data:
        ismulti = length(find(handles.Data.used.config.OK==true));
        models = handles.Models.model;
        if ismulti > 1,
            i = 1;
            str_results = [''];
            while i <= ismulti,
                if type_mod == 2,% Complex
                    eval(['tmp1 = real(models.results' num2str(i) '); tmp2 = imag(models.results' num2str(i) ');']);
                    tmp_global = zeros(N, size(tmp1,2)*2,size(tmp1,3));
                    tmp_global(:,1:end/2,:)=tmp1;
                    tmp_global(:,end/2+1:end,:)=tmp2;
                    eval(['models.results' num2str(i) ' = tmp_global;']);
                    if i < ismulti,
                        str_results = [str_results [' models.results' num2str(i) ',']];
                    else
                        str_results = [str_results [' models.results' num2str(i)]];
                    end
                    i = i + 1;
                elseif type_mod == 3,% Real only
                    if i < ismulti,
                        str_results = [str_results [' real(models.results' num2str(i) '),']];
                    else
                        str_results = [str_results [' real(models.results' num2str(i) ')']];
                    end
                    i = i + 1;
                else % Absolute value
                    if i < ismulti,
                        str_results = [str_results [' abs(models.results' num2str(i) '),']];
                    else
                        str_results = [str_results [' abs(models.results' num2str(i) ')']];
                    end
                    i = i + 1;
                end
            end
            eval(['models.results = [' (str_results) '];']);
        else
            if type_mod == 2,
                tmp1 = real(models.results1);
                tmp2 = imag(models.results1);

                tmp_global = zeros(N, size(tmp1,2)*2,size(tmp1,3));
                tmp_global(:,1:end/2,:)=tmp1;
                tmp_global(:,end/2+1:end,:)=tmp2;

                models.results = tmp_global;
                clearvars tmp1 tmp2 tmp_global;
            elseif type_mod == 3,
                models.results = real(models.results1);
            else
                models.results = abs(models.results1);
            end
        end
        
        handles.Models.model.results = models.results;
        for i = 1 : ismulti,
            % For memory gains !
            handles.Models.model = rmfield(handles.Models.model,['results' num2str(i)]);
        end
        guidata(hObject, handles);

        % 2.b) Real data:
        d_real_obs = [];
        if ismulti > 1,
            i = 1;
            while i <= ismulti,
                eval(['proclog = ' handles.Data.proclogNamesCode{handles.Data.used.config.channelTx(i)} ';']);
                taille = size([proclog.Q(:).q],2);
                for j = 1 : taille,
                    eval(['d_real_obs = [d_real_obs proclog.Q(j).rx(' num2str(handles.Data.used.config.channelRx(i)) ').sig(2).V];']);
                end
                i = i + 1;
            end
        else
            eval(['proclog = ' handles.Data.proclogNamesCode{handles.Data.used.config.channelTx(1)} ';']);
            t_register = [proclog.Q(1).rx(1).sig(2).t];
            d_real_obs = zeros(length([proclog.Q(:).q]),length(t_register));
            taille = size([proclog.Q(1).rx(handles.Data.used.config.channelRx(1)).sig(2).V],2);
            for i = 1 : size(d_real_obs,1),
                d_real_obs(i,1:taille) = proclog.Q(i).rx(handles.Data.used.config.channelRx(1)).sig(2).V;
            end
        end
        if type_mod == 2,
            d_real_obs = [real(d_real_obs); imag(d_real_obs)];
        elseif type_mod == 3,
            d_real_obs = real(d_real_obs);
        else
            d_real_obs = abs(d_real_obs);
        end
        tmp_real = [];
        for i = 1 : size(d_real_obs,1),
            tmp_real = [tmp_real d_real_obs(i,:)];
        end

        d_real_obs = tmp_real;
        handles.Data.d_real_obs = d_real_obs;
        guidata(hObject, handles);

        %% 3) Dimension reduction

        % If a parallel pool does exist, it will take a large part of the
        % memory available... It is therefor shut down! If the user wants to
        % use parallel computing later, they will need to restart the pool.
        p = gcp('nocreate');
        if ~isempty(p),
            poolsize = p.NumWorkers;
            delete(gcp('nocreate'));
            fprintf('The parallel pool (%d workers) was shut down for memory gains! \n Please restart the parallel pool later if needed.',poolsize);
        end

        waitbar2a(0.5,handles.uipanel_waitbar,sprintf('Reducing dimensionality \n on FID. . . '));
        models = handles.Models.model;
        dobs_f = handles.Data.d_real_obs;

        Prior_h = [models.thick, models.water, (models.T2)];
        Prior_d = [];
        for i = 1 : size(models.results,2),
            Prior_d = [Prior_d squeeze(models.results(:,i,:))];
        end

        % 3.a) PCA on model parameters --> not done (useless)

        PCA_h = Prior_h - repmat(mean(Prior_h,1),size(Prior_h,1),1);
        coeff_h = eye(size(Prior_h,2));
        dimh = size(PCA_h,2);

        % 3.b) PCA on data

        level_var = 0.1;
        warning('off','all');
        [coeff_d,score,~,~,explained_F,~]=pca(Prior_d);
        warning('on','all');
        nb_PC = max([find(explained_F>level_var,1,'last'),dimh]);% To be sure that there is at least more or equal d than h dim
        PCA_d = score(:,1:nb_PC);
        dimd = nb_PC;

        dobs_f = (d_real_obs-(mean(Prior_d)))*coeff_d(:,1:nb_PC);


        % 3.c) CCA

        % set(handles.text_progress,'String','Canonical correlation analysis');
        waitbar2a(0.5,handles.uipanel_waitbar,sprintf('Canonical correlation \n analysis . . . '));
        [A,B,~,Dc,Hc] = canoncorr(PCA_d,PCA_h);
        dobs_c = (dobs_f-mean(PCA_d))*A;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %
%  %                                                                   %  %
%  %           /\        Printing the graph on the influence of each   %  %
%  %          /  \       parameter on the different CCA dimensions     %  %
%  %         /  | \                                                    %  %
%  %        /   |  \     COMMENT if unwanted!!!                        %  %
%  %       /    |   \                                                  %  %
%  %      /     .    \                                                 %  %
%  %     /____________\                                                %  %
%  %                                                                   %  %
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure('Units','normalized','Position',[0.05 0.05 0.8 0.8])
        subplots = zeros(size(B,1),2);
        for jkl = 1 : size(B,1),
            if jkl <= ceil(size(B,1)/2),
                subplots(jkl,:) = 2*(jkl-1)+1:2*(jkl-1)+2;
            else
                subplots(jkl,:) = 2*(jkl-1)+1:2*(jkl-1)+2;
            end
        end
        for jkl = 1 : size(B,1),
            subplot(2,ceil(size(B,1)/2)*2,subplots(jkl,:));
            pos = get(gca, 'Position');
            pos(3:4) = 0.8*pos(3:4);
            set(gca, 'Position', pos)
            scatter(Dc(:,jkl),Hc(:,jkl));
            hold on;
            plot([dobs_c(jkl) dobs_c(jkl)],[-3 3],'linewidth',3);
            xlabel(['d^c_' num2str(jkl)])
            ylabel(['m^c_' num2str(jkl)])
            set(gca,'FontSize',16)
        end
        h = figure;
        nb_layer = handles.Models.nbLayers;
        for i = 1 : (3*nb_layer-1),
           c{i} = ['M^c_' num2str(i)];
           if i < nb_layer,
               leg{i} = ['e_' num2str(i)];
           elseif i< 2*nb_layer,
               leg{i} = ['W_' num2str(i-nb_layer+1)];
           else
               leg{i} = ['T^*_{2,' num2str(1+i-nb_layer*2) '}'];
           end
        end   
        bar((abs(B)'./repmat(sum(abs(B)',1),size(B,1),1)),'stacked');
        set(gca,'xticklabel',c);
        set(gca,'yticklabel',[]);
        legend(leg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------- End: graph display ---------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Verification of noise impact on a sample of prior models
        % Default bandwidth for the kernel density estimator:
        default_bandwidth = 0.01;
        answer = questdlg('Would you like to test the impact of noise?','Noise impact','Yes','No','No');
        if strcmp(answer,'Yes'),
            [~, error] = testNoise(Prior_d, PCA_d, A, default_bandwidth);
        else
            error = ones(dimh,1).*default_bandwidth;
        end
        %% 4) Kernel density 

        disc = 10000;% Number of points in the pdf
        values = -10 : 20/(disc-1) : 10;
        for i = 1 : size(Dc,2)
            waitbar2a(i/(size(Dc,2)*2),handles.uipanel_waitbar,sprintf('Constructing posterior \n in reduced space . . .'));
            [fi(:,i),out] = KernelDensity([Dc(:,i) Hc(:,i)],dobs_c(i),values,[error(i) 0.01]); 
            if out,
                if iterate,
                    break
                end
            end
        end
        
        if ~out | ~iterate,
            cdf = fi.*0;
            for i = 2 : disc,
                cdf(i,:) = trapz(values(1:i),fi(1:i,:));
            end

            %% 5) Sampling and back-transformation

            nb_sample = handles.Models.N;
            nb_layer = handles.Models.nbLayers;
            CCAi = zeros(nb_sample,size(fi,2));
            HpostCoeff = zeros(nb_sample, size(fi,2));
            nb_trial = 0;
            i = 0;
            % For the rejection sampler:
            parameters = get(handles.uitable_prior,'Data');
            parameters = parameters(1:nb_layer,:);
            parameters(:,3:4) = parameters(:,3:4)./100;
            param_min = parameters(:,1:2:end);
            param_max = parameters(:,2:2:end);
            tmp = 1;
            for k = 1 : size(fi,2),
                if k == nb_layer,
                    tmp = tmp + 1;
                end
                Pmin(k) = param_min(tmp);
                Pmax(k) = param_max(tmp);
                tmp = tmp + 1;
            end
            while i < nb_sample,
                i = i + 1;
                waitbar2a(i/nb_sample,handles.uipanel_waitbar,'Sampling and back-transformation . . .');
                for j = 1 : size(cdf,2),
                    ni = rand(1);
                    idx = find(cdf(:,j)<=ni,1,'last');
                    CCAi(i,j) = values(idx);
                end
                HpostCoeff(i,:) = CCAi(i,:)*pinv(B)+repmat(mean(PCA_h,1)',1,1)';% Return to PCA space
                HpostCoeff(i,:) = HpostCoeff(i,:)*pinv(coeff_h)+repmat(mean(Prior_h,1)',1,1)';% Return to real space
                if not(isempty(find(HpostCoeff(i,:)<0, 1))) || not(isempty(find(HpostCoeff(i,:)<Pmin,1))) || not(isempty(find(HpostCoeff(i,:)>Pmax,1))),
                    i = i - 1;
                    nb_trial = nb_trial + 1;
                    if nb_trial >= 100000,
                        warning('Unable to sample correctly! \n \t The data are outside the physical boundaries of the domain!\n');
                        break
                    end
                end
            end

            depth_tmp = zeros(size(HpostCoeff,1),2);
            for i = 2 : nb_layer,
                depth_tmp(:,i) = depth_tmp(:,i-1) + HpostCoeff(:,i-1);
            end

            models = struct('thick',HpostCoeff(:,1:nb_layer-1),'depth',depth_tmp,...
                'water',HpostCoeff(:,nb_layer:2*nb_layer-1),'T2',HpostCoeff(:,2*nb_layer:end),'res',ones(size(HpostCoeff,1),1)*100);
            clear depth_tmp;

            handles.Solution.model = models;


            waitbar2a(1,handles.uipanel_waitbar,'Finalizing!');

            passed = 1;
            iteration = iteration - 1;
        else
            waitbar2a(0.5,handles.uipanel_waitbar,sprintf('Comming back to previous\n iteration...'));
            pause(1)
            iteration = 0;
        end
    end
    if iteration_init>1,
        handles.Models.model = handles.Models.model_save;
    end
    handles.Status.PFA = true;
    guidata(hObject, handles);
end

% GUI adaptation:
set(handles.uipanel_waitbar,'Visible','off');

set(handles.uibuttongroup_analyze,'Visible','on');
set(handles.pushbutton_save,'Visible','on');
% End computing PFA


% --- Executes on button press in pushbutton_Show_post.
function pushbutton_Show_post_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Show_post (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

models = handles.Solution.model;
models_prior = handles.Models.model;
nb_layer = handles.Models.nbLayers;
figure('Units','normalized','OuterPosition',[0.1 0.1 0.5 0.6]);
subplot(nb_layer,6, 6*nb_layer-1:6*nb_layer);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
position = get(gca,'Position');
for j = 1 : nb_layer,
    subplot(nb_layer,6,(j-1)*6+1:(j-1)*6+2);% Water content
    if j == 1,
        title('Water content [%]','FontSize',16);
    end
    hold on;
    histogram(models_prior.water(:,j).*100,10,'Normalization','pdf');
    histogram(models.water(:,j).*100,'Normalization','pdf');
    ylabel(['Layer ' num2str(j)],'FontSize',14);
    subplot(nb_layer,6,(j-1)*6+3:(j-1)*6+4);% T2
    if j == 1,
        title('Relaxation time [ms]','FontSize',16);
    end
    hold on;
    histogram(models_prior.T2(:,j),10,'Normalization','pdf');
    histogram(models.T2(:,j),'Normalization','pdf');
    subplot(nb_layer,6,(j-1)*6+5:(j-1)*6+6);
    if j == 1,
        title('Thickness [m]','FontSize',16);
    end
    hold on;
    if j ~= nb_layer,
        histogram(models_prior.thick(:,j),10,'Normalization','pdf');
        histogram(models.thick(:,j),'Normalization','pdf');
    end
    if j == nb_layer-1,
        w = legend('Prior','Posterior');
        w.Position = position;
        w.FontSize = 12;
    end
end

param = [models.water models.T2 models.thick];
param_prior = [models_prior.water models_prior.T2, models_prior.thick];
nb_param = 3;
param_pre = {'W_{','T^*_{2,','e_{'};
units = {'[m^3/m^3]','[ms]','[m]'};
param_names = {};
for i = 1 : nb_param,
    for j = 1 : nb_layer,
        param_names{(i-1)*nb_layer+j} = [param_pre{i} num2str(j) '} ' units{i}];
    end
end
nb_param = nb_param*nb_layer-1;
figure('Units','centimeters','OuterPosition',[1.5 1.5 25 25]);
subplot(nb_param,nb_param,nb_param);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
position = get(gca,'Position');
for i = 1 : nb_param,
    for j = 1 : nb_param,
        if j < i,
            subplot(nb_param,nb_param,(i-1)*nb_param+j);
            hold on;
            scatter(param_prior(:,j),param_prior(:,i),'.');
            scatter(param(:,j),param(:,i),'.');
            if j == 1,
                ylabel(param_names{i});
            end
            if i == nb_param,
                xlabel(param_names{j});
            end
        elseif i == j,
            subplot(nb_param,nb_param,(i-1)*nb_param+j);
            hold on;
            histogram(param_prior(:,j),10,'DisplayStyle','stairs','Normalization','pdf');
            histogram(param(:,j),'DisplayStyle','stairs','Normalization','pdf');
            if i == nb_param,
                xlabel(param_names{j});
                w = legend('Prior','Posterior');
                w.Position = position;
                w.FontSize = 12;
            elseif j == 1,
                ylabel(param_names{j});
            end
        end
    end
end



% --- Executes on button press in pushbutton_RMS.
function pushbutton_RMS_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_RMS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
addpath([pwd '/SNMR']);
%% 1) Forward modelling:
models = handles.Solution.model;
i = 1;
type_mod = 1;% Absolute value
while i <= length(find(handles.Data.used.config.OK==true)),
    eval(['kdata = handles.Kernels.' handles.Data.used.config.KNamesCode{i} ';']);
    eval(['proclog = ' handles.Data.proclogNamesCode{handles.Data.used.config.channelTx(i)} ';']);
    t_register = [proclog.Q(1).rx(1).sig(2).t];    
    for j = 1 : length(models.thick(:,1)),
        mdata    = struct();% Initializing
        mdata.mod.nTau = 2;
        mdata.mod.tau  = [0.1 0.3];
        mdata.mod.Nlayer = handles.Models.nbLayers;
        if models.depth(i,end) < kdata.model.zmax
            mdata.mod.zlayer = [models.depth(j,:) kdata.model.zmax];
        else
            mdata.mod.zlayer = [models.depth(j,1:end-1) kdata.model.zmax-1 kdata.model.zmax];
        end
        mdata.mod.f      = models.water(j,:);
        mdata.mod.T2s    = models.T2(j,:)./1e3;
        mdata.mod.T2     = zeros(1, length(mdata.mod.zlayer)-1);
        mdata.mod.T1     = zeros(1, length(mdata.mod.zlayer)-1);

        mdata.mod.tfid1  = t_register + proclog.Q(1).timing.tau_dead1;
        mdata.mod.noise  = 0;

        kdata.KT1 = [];

        mdata = ForwardModelling_FID(mdata,kdata);
        if type_mod == 1,%% Absolute values
            data_model.fid(j,:,:) = abs(mdata.dat.fid1);
        elseif type_mod == 2%% Complex values
            data_model.fid(j,:,:) = mdata.dat.fid1;
        elseif type_mod == 3,%% Real only
            data_model.fid(j,:,:) = real(mdata.dat.fid1);
        else
            data_model.fid(j,:,:) = abs(mdata.dat.fid1);
        end
    end
    eval(['handles.Solution.model.results' num2str(i) ' = data_model.fid;']);
    i = i + 1;
end
clearvars('kdata','proclog');
guidata(hObject, handles);

%% 2) Grouping real and simulated data:
% 2.a) Simulated data:
ismulti = length(find(handles.Data.used.config.OK==true));
models = handles.Solution.model;

if ismulti > 1,
    i = 1;
    str_results = [''];
    while i <= ismulti,
        if type_mod == 2,% Complex
            eval(['tmp1 = real(models.results' num2str(i) '); tmp2 = imag(models.results' num2str(i) ');']);
            tmp_global = zeros(N, size(tmp1,2)*2,size(tmp1,3));
            tmp_global(:,1:end/2,:)=tmp1;
            tmp_global(:,end/2+1:end,:)=tmp2;
            eval(['models.results' num2str(i) ' = tmp_global;']);
            if i < ismulti,
                str_results = [str_results [' models.results' num2str(i) ',']];
            else
                str_results = [str_results [' models.results' num2str(i)]];
            end
            i = i + 1;
        elseif type_mod == 3,% Real only
            if i < ismulti,
                str_results = [str_results [' real(models.results' num2str(i) '),']];
            else
                str_results = [str_results [' real(models.results' num2str(i) ')']];
            end
            i = i + 1;
        else % Absolute value
            if i < ismulti,
                str_results = [str_results [' abs(models.results' num2str(i) '),']];
            else
                str_results = [str_results [' abs(models.results' num2str(i) ')']];
            end
            i = i + 1;
        end
    end
    eval(['models.results = [' (str_results) '];']);
else
    if type_mod == 2,
        tmp1 = real(models.results1);
        tmp2 = imag(models.results1);

        tmp_global = zeros(N, size(tmp1,2)*2,size(tmp1,3));
        tmp_global(:,1:end/2,:)=tmp1;
        tmp_global(:,end/2+1:end,:)=tmp2;

        models.results = tmp_global;
        clearvars tmp1 tmp2 tmp_global;
    elseif type_mod == 3,
        models.results = real(models.results1);
    else
        models.results = abs(models.results1);
    end
end

clearvars tmp1 tmp2 tmp_global;

handles.Solution.model.results = models.results;
for i = 1 : ismulti,
    % For memory gains !
    handles.Solution.model = rmfield(handles.Solution.model,['results' num2str(i)]);
end
clearvars tmp1 tmp2 tmp_global;
guidata(hObject, handles);

data_true = handles.Data.d_real_obs;
data_post = handles.Solution.model.results;
tmp = [];
for i = 1 : size(data_post,2),
    tmp = [tmp squeeze(data_post(:,i,:))];
end
data_post = tmp;
clear tmp;

%% 3) Computing RMS
RMS = rms(data_post.*1e9-repmat(data_true.*1e9,size(data_post,1),1),2);

[RMS,I] = sort(RMS,'descend');

RMS_fact = (RMS)./(max(RMS)+0.005);

HpostCoeff = [models.thick, models.water, (models.T2)];
HpostCoeff = HpostCoeff(I,:);
%% 4) Presenting the model in 1D profil

map = colormap(jet);
% Adding feature for the same number of models per share of color
nb_map = length(map);
Q_used = (0 : nb_map-1)./(nb_map-1);
RMS_Interval = quantile(RMS,Q_used);
RMS_fact_map = ones(1,length(RMS));
Q_act = 1;
for i = length(RMS) : -1 : 2,
    RMS_fact_map(i) = Q_act;
    if RMS(i-1) >= RMS_Interval(Q_act),
        Q_act = Q_act + 1;
    end
end
RMS_fact_map(1) = nb_map;
% End addition!
% RMS_fact_map = ceil((RMS_fact).*(length(map)/(max((RMS_fact))-min(RMS_fact))));
% RMS_fact_map = RMS_fact_map - min(RMS_fact_map);
% RMS_fact_map(RMS_fact_map == 0) = 1;

figure;
subplot(10,2,[1:2:13]);% Water content
hold on;
for i = 1 : 10 : size(HpostCoeff,1),
    stairs([models.water(I(i),1) models.water(I(i),:)],[models.depth(I(i),:) 100],'Color',map(RMS_fact_map(i),:),'LineWidth',(1-RMS_fact(i))*3);
end
set (gca,'Ydir','reverse');
xlabel('Water content (\phi) [/]');
xlim([0 0.5]);
ylim([0 30]);
ylabel('Depth [m]');
grid on;
hold off;

subplot(10,2,[2:2:14]);% Relaxation time
hold on;
for i = 1 : 10 : size(HpostCoeff,1),
    stairs([models.T2(I(i),1) models.T2(I(i),:)],[models.depth(I(i),:) 100],'Color',map(RMS_fact_map(i),:),'LineWidth',(1-RMS_fact(i))*3);
end
set (gca,'Ydir','reverse');
xlabel('Relaxation time (T_2^*) [ms]');
xlim([0 500]);
ylim([0 30]);
ylabel('Depth [m]');
grid on;
hold off;

subplot(10,2,19:20);
%Adding feature for same number of models per color
tmp = min(RMS) : (max(RMS)-min(RMS))/999 : max(RMS);
for i = 1 : length(tmp),
    tmp(i) = find(RMS_Interval<=tmp(i),1,'last');
end
imagesc(RMS_Interval,0:1,[(tmp)' (tmp)']');
% End addition
% imagesc([min(RMS) : (max(RMS)-min(RMS))/(size(map,1)-1) : max(RMS)],0:1,[(1:size(map,1))' (1:size(map,1))']');
colormap(map);
set(gca,'YTick',[]);
xlabel('RMS error [nV]');


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SNMR.Data = handles.Data;
SNMR.Models = handles.Models;
SNMR.Kernels = handles.Kernel;
SNMR.Solution = handles.Solution;
save([handles.Files.path '\Results_SNMR.mrsbel'],'SNMR','-v7.3');
[file,path] = uiputfile('*.mrsbel');
resp = movefile([handles.Files.path '\Results_SNMR.mrsbel'],[path file]);
if resp,
    fprintf('File save at : \n \t %s%s \n',path,file);
else
    fprintf('Error while saving the file! Try again . . .');
end


% --- Executes on selection change in listboxTx.
function listboxTx_Callback(hObject, eventdata, handles)
% hObject    handle to listboxTx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxTx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxTx
for i = 1 : 4,
    eval(['set(handles.CB_Rx' num2str(i) ',''Visible'',''off'');']);
end
item = get(hObject,'Value');
i = 1;
while handles.Data.used.Rx(item,i) > 0,
    eval(['set(handles.CB_Rx' num2str(i) ',''String'',[''Rx = '' num2str(handles.Data.used.Rx(item,i)) '' m'']);']);
    eval(['set(handles.CB_Rx' num2str(i) ',''Value'',handles.Data.activated(item,i));']);
    eval(['set(handles.CB_Rx' num2str(i) ',''Visible'',''on'');']);
    i = i + 1;
end


% --- Executes during object creation, after setting all properties.
function listboxTx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxTx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CB_Rx1.
function CB_Rx1_Callback(hObject, eventdata, handles)
% hObject    handle to CB_Rx1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CB_Rx1

i = get(handles.listboxTx,'Value');
handles.Data.activated(i,1) = get(hObject,'Value');
guidata(hObject, handles);


% --- Executes on button press in CB_Rx2.
function CB_Rx2_Callback(hObject, eventdata, handles)
% hObject    handle to CB_Rx2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CB_Rx2
i = get(handles.listboxTx,'Value');
handles.Data.activated(i,2) = get(hObject,'Value');
guidata(hObject, handles);


% --- Executes on button press in CB_Rx3.
function CB_Rx3_Callback(hObject, eventdata, handles)
% hObject    handle to CB_Rx3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CB_Rx3
i = get(handles.listboxTx,'Value');
handles.Data.activated(i,3) = get(hObject,'Value');
guidata(hObject, handles);


% --- Executes on button press in CB_Rx4.
function CB_Rx4_Callback(hObject, eventdata, handles)
% hObject    handle to CB_Rx4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CB_Rx4
i = get(handles.listboxTx,'Value');
handles.Data.activated(i,4) = get(hObject,'Value');
guidata(hObject, handles);


% --- Executes on button press in pushbutton_file.
function pushbutton_file_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Files.path = uigetdir;
handles.edit_data_path.String = handles.Files.path;
guidata(hObject, handles);
