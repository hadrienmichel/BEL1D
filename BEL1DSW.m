function varargout = BEL1DSW(varargin)
% BEL1DSW is a graphical interface that helps the user to perform
% prediction-focused approach imaging of the subsurface using Analysis of 
% Surface Waves data through the fundamental mode of the dispersion curve.
% 
%
% This code relies heavily on the gpdc code from GEOPSY 
% [http://www.geopsy.org](Wathlet, 2008) that must be accessible to matlab 
% (add the Geopsy folder to the general path in windows).
% The GUI is organized in multiple panels that appears progressively and
% enables the BEL modeling of the dataset in an intuitive way.
%
% The panels are:
%
%   1) Data panel: This panel enables the user to load the data, either in
%                  the '*.maswd' format (see example below) or in the '*'
%                  format issued from the f-k analysis of Geopsy.
%
%                  o The '*.maswd' file format is a '-mat' file that 
%                    contains a structure called 'dispersion' and contains 
%                    the next elements:
%                   >> dispersion.Curve -> the fundamental mode dispersion
%                                          curve.
%                           dispersion.Curve.f ->  the frequencies (vector)
%                           dispersion.Curve.ct -> the corresponding phase
%                                                  velocities (vector)
%                           dispersion.Curve.lambda -> the wavenumbers
%                                                      (vector)
%                   >> dispersion.FK (optional) -> the full fk image
%                           dispersion.FK.f -> the frequencies (vector)
%                           dispersion.FK.ct -> the phase velocities
%                                               (vector)
%                           dispersion.FK.normed -> the normed amplitude 
%                                                   (matrix)
%                   >> disperstion.true (optional) -> the true model
%                           dispersion.true.thick -> the thickness of the
%                                                    N-1 first layers
%                                                    (vector)
%                           dispersion.true.Vp -> the P-wave velocity of
%                                                 the N layers (vector)
%                           dispersion.true.Vs -> the S-wave velocity of
%                                                 the N layers (vector)
%                           dispersion.true.rho -> the density of the N
%                                                  layers (vector)
%
%                   o The '*' (no extension) file format (Geopsy) is a 
%                     *space-separated* multi-column text file organized as
%                     described below.
%                       # 
%                       # Some
%                       # commentary
%                       # lines
%                       #
%                       # | Frequency (Hz) | Slowness (s/m) | Stddev (s/m) | Weight | 
%                       10 0.001 0 1 
%                       10.1638981180644 0.001 0 1 
%                       10.3304824954393 0.00107214428857715 0 1 
%                       10.4997971594093 0.00123446893787575 0 1 
%                       10.6718868588578 0.00128857715430862 0 1 
%                       10.8467970760941 0.00132464929859719 0 1 
%                       11.0245740388739 0.00145090180360721 0 1 
%                       11.2052647326172 0.00177555110220441 0 1 
%                       11.3889169128261 0.00215430861723447 0 1 
%                       11.5755791177065 0.00237074148296593 0 1 
%                       11.7653006809963 0.00251503006012024 0 1 
%                       . . . etc.
%
%                  Once the data is loaded properly, the next panel,
%                  enabling a tuning of the prior model space, appears.
%
%   2) Prior model space panel: In this panel, the user is invited to
%                               describe the prior model space through the
%                               different parameters that are proposed. The
%                               table sets the limits of the domain, the 
%                               number of model generated controls the
%                               number of sampled models, the type of
%                               distribution/sampler is set by selecting
%                               the type in the pop-up menu.
%                               Once all the parameters are set optimaly,
%                               the user is invited to click on 'Generate
%                               models'. This will enable the next panel to
%                               appear.
%
%   3) BEL modeling panel: In this panel, the BEL modeling is taking place.
%                         By clicking on the 'Run BEL modeling' pushbutton,
%                         the user launches the computations. A progressbar
%                         will appear and the different steps advancement
%                         described. 
%                         Once all the computations are done, the user can 
%                         analyze the results through the parameters 
%                         distributions or the model distributions. The 
%                         user can also save the obtained results in a 
%                         '*.maswbel' file for further computations. The 
%                         saved file will contain the simulated data for 
%                         the posterior space only if the model 
%                         distributions have been displayed previously.
%
% For further help, please contact the author:
%       Name: Hadrien MICHEL
%       University: University of Liège and University of Ghent
%       Country: Belgium
%       e-mail: Hadrien.Michel@uliege.be
%
% © 2019, Hadrien MICHEL, University of Liège, University of Ghent
% and F.R.S.-FNRS
%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BEL1DSW_OpeningFcn, ...
                   'gui_OutputFcn',  @BEL1DSW_OutputFcn, ...
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


% --- Executes just before BEL1DSW is made visible.
function BEL1DSW_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)
addpath([pwd '/GlobalFunctions']);

% Choose default command line output for BEL1DSW

handles.output = hObject;

% Adding the path to the MASWave functions (must be contained in the same
% directory):
%addpath(genpath(pwd),'-begin');

% Custom variables that are contained in the GUI:
handles.Files = [];% Contains the files informations
handles.Files.path = 'Files path';
handles.Files.name = 'Name File';
handles.Data = [];% Contains the data-sets
handles.Models = [];
handles.Models.param = [2.5 7.5 300 1200 NaN NaN 100 500 1500 2000;0 5 1200 2500 NaN NaN 400 1200 1900 2400;0 5 100 2000 NaN NaN 100 1500 1000 3000;...
    0 5 100 2500 NaN NaN 100 2000 1000 3000;0 5 100 3000 NaN NaN 100 2500 1000 3000;0 5 100 3500 NaN NaN 100 3000 1000 3000];
handles.Models.nbLayers = 2;
set(handles.uitable_prior,'BackgroundColor',[repmat([1 1 1],handles.Models.nbLayers,1);repmat([0 0 0],6-handles.Models.nbLayers,1)]);
handles.Models.type = 1;
handles.Models.N = 1000;
handles.Solution = [];

% Graphical tuning of the interface
set(handles.uipanel_model,'Visible','off');
set(handles.uipanel_PFA,'Visible','off');

set(handles.axes2,'Visible','off');

% Adding the logos of the Universtity of Liège and Ghent:
axes(handles.axes3)
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

axes(handles.axes4)
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
handles.Status.models = false;%Model generated?
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
function varargout = BEL1DSW_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;  

function edit_data_path_Callback(hObject, eventdata, handles)
% hObject    handle to edit_data_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Enables the loading of the data if the user has entered manually the path
% to the data file. The function performs:
%   - Loading of the data
%   - Graphical display of the data

%
% Loading of the data:
%
tmp = get(hObject,'String');
[filepath,name,ext] = fileparts(tmp);
% Need to differenciate between the geopsy files and the mat files:
handles.Files.path = filepath;
handles.Files.name = [name ext];

if strcmp(ext,'.maswd'), % Matlab generated file
    dispersion = load([handles.Files.path handles.Files.name],'-mat');
    handles.Data = dispersion.dispersion;
elseif strcmp(ext,''), % Geopsy generated file
    fileID = fopen(handles.edit_data_path.String,'r');
    c = textscan(fileID, '%f %f %f %f','Delimiter',' ','CommentStyle','#');
    fclose(fileID);
    dispersion.Curve.f = c{1};
    dispersion.Curve.ct = [c{2}].^(-1);
    dispersion.Curve.lambda = dispersion.Curve.ct.*dispersion.Curve.f;
    handles.Data = dispersion;
else
    error('The file you selected is not readable by the GUI!');
end

%
% Displaying the data
%
axes(handles.axes2);
if isfield(handles.Data,'FK'),% If the full f-k analysis is available it is displayed in background
    h=pcolor(handles.Data.FK.f,handles.Data.FK.ct,handles.Data.FK.normed);
    set(h, 'EdgeColor', 'none');
    colormap(jet)
end
hold on;
plot(handles.Data.Curve.f,handles.Data.Curve.ct,'k','linewidth',2');
set(handles.axes2,'Visible','on');
xlabel('Frequency [Hz]');
ylabel('Phase Velocity [m/s]');

set(handles.uitable_prior,'Data',handles.Models.param); 
set(handles.uipanel_model,'Visible','on');
handles.Status.data = true;
handles.Models.increasing = false;
guidata(hObject, handles);

% --- Executes on button press in pushbutton_file.
function pushbutton_file_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Enables the loading of the data if the user has entered the path to the 
% data file through the GUI. The function performs:
%   - Loading of the data
%   - Graphical display of the data

%
% Loading of the data:
%
[handles.Files.name ,handles.Files.path] = uigetfile({'*.maswd','MATLAB MASW data file (*.maswd)';'*','Geopsy file (no extension)'},'Select a data file');
[~,~,ext] = fileparts([handles.Files.path handles.Files.name]);
handles.edit_data_path.String = [handles.Files.path handles.Files.name];
if strcmp(ext,'.maswd'), % Matlab generated file
    dispersion = load([handles.Files.path handles.Files.name],'-mat');
    handles.Data = dispersion.dispersion;
elseif strcmp(ext,''), % Geopsy generated file
    fileID = fopen(handles.edit_data_path.String,'r');
    c = textscan(fileID, '%f %f %f %f','Delimiter',' ','CommentStyle','#');
    fclose(fileID);
    dispersion.Curve.f = c{1};
    dispersion.Curve.ct = [c{2}].^(-1);
    dispersion.Curve.lambda = dispersion.Curve.ct./dispersion.Curve.f;
    handles.Data = dispersion;
else
    error('The file you selected is not readable by the GUI!');
end

%
% Displaying the data
%
axes(handles.axes2);
if isfield(handles.Data,'FK');
    h=pcolor(handles.Data.FK.f,handles.Data.FK.ct,handles.Data.FK.normed);
    set(h, 'EdgeColor', 'none');
    colormap(jet)
end
hold on;
plot(handles.Data.Curve.f,handles.Data.Curve.ct,'k','linewidth',2');
set(handles.axes2,'Visible','on');
xlabel('Frequency [Hz]');
ylabel('Phase Velocity [m/s]');

set(handles.uitable_prior,'Data',handles.Models.param); 
set(handles.uipanel_model,'Visible','on');
handles.Status.data = true;
handles.Models.increasing = false;
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
addpath([pwd '/SurfaceWave']);
% Function for model generation:
% [models] = ModelGenerator(type, N, parameters, nb_layer)
parameters = get(handles.uitable_prior,'Data');
parameters = parameters(1:handles.Models.nbLayers,:);
handles.Models.model = ModelGeneratorMASW(handles.Models.type, handles.Models.N, parameters, handles.Models.nbLayers,handles.Models.increasing);
handles.Status.models = true;
if handles.Status.models,
    set(handles.uipanel_PFA,'Visible','on');
end
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Run.
function pushbutton_Run_Callback(hObject, eventdata, handles)
% This function performs the main object of the PFA imaging. It is
% decomposed in several steps:
%   1) Forward modelling
%       This steps consits on the computation of the simulated data for
%       each sampled model from the prior model space. It uses the function
%       'MASWaves_theoretical_dispersion_curve_ModHadrien.m' (a modified 
%       version of the MASWave 'MASWaves_theoretical_dispersion_curve.m'
%       function : Olafsdottir et al., 2018).
%   2) Reducing dimensionality
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
%   3) Kernel density estimation
%       This transforms the set of points in the CCA space into
%       probabilistic distributions that can be sampled for the next step.
%   4) Sampling and back-transformation (rejection of samples outside the
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
addpath([pwd '/SurfaceWave']);
addpath([pwd '/GlobalFunctions']);
if handles.Status.PFA ~= true,
    set(handles.uipanel_waitbar,'Visible','on');
    drawnow;
    % Computing PFA

    %% 1) Forward modelling:
    
    waitbar2a(0,handles.uipanel_waitbar,'waitbartext',sprintf('Please wait . . . \n Computing Dispersion Curves'));

%     waitbar2a(0,handles.uipanel_waitbar,sprintf('Please wait . . . \n Computing forward models'));
% 
%     lambda_t = zeros(handles.Models.N,length(handles.Data.Curve.lambda));
%     c_t = lambda_t;
%     if isempty(gcp('nocreate')),
%         for j = 1 : handles.Models.N,
%             waitbar2a(j/handles.Models.N,handles.uipanel_waitbar);
%             [c_t(j,:),lambda_t(j,:)] = MASWaves_theoretical_dispersion_curve_ModHadrien(10 : 1 : 1500,handles.Data.Curve.lambda,handles.Models.model.thick(j,:),handles.Models.model.Vp(j,:),handles.Models.model.Vs(j,:),handles.Models.model.rho(j,:),handles.Models.nbLayers-1);
%         end
%     else
%         tested = [1 100 : 100 :handles.Models.N];
%         hjkl = 1;
%         lambda = handles.Data.Curve.lambda;
%         thick = handles.Models.model.thick;
%         Vp = handles.Models.model.Vp;
%         Vs = handles.Models.model.Vs;
%         rho = handles.Models.model.rho;
%         nbLayers = handles.Models.nbLayers;
%         for jklm = tested(2 : end),
%             waitbar2a(jklm/handles.Models.N,handles.uipanel_waitbar);
%             parfor j = tested(hjkl) : jklm,
%                 [c_t(j,:),~] = MASWaves_theoretical_dispersion_curve_ModHadrien(10 : 1 : 1500,lambda,thick(j,:),Vp(j,:),Vs(j,:),rho(j,:),nbLayers-1);
%             end
%             hjkl = hjkl + 1;
%         end
%     end
%     handles.Models.model.results = c_t;

    [handles.Models,out_error] = gpdcCall_bisMEX(handles.Models,handles.Data.Curve.f);
    handles.Models.N = handles.Models.N - sum(out_error);
    handles.Models.model.thick = handles.Models.model.thick(~out_error,:);
    handles.Models.model.Vp = handles.Models.model.Vp(~out_error,:);
    handles.Models.model.Vs = handles.Models.model.Vs(~out_error,:);
    handles.Models.model.rho = handles.Models.model.rho(~out_error,:);
    handles.Models.model.results = handles.Models.model.results(~out_error,:);
    
    guidata(hObject, handles);

    %% 2) Dimension reduction

    waitbar2a(0.5,handles.uipanel_waitbar,sprintf('Reducing dimensionality \n on data. . . '));
    models = handles.Models.model;
    d_real_obs = handles.Data.Curve.ct;

    Prior_h = [models.thick, models.Vp, models.Vs, models.rho];
    Prior_d = handles.Models.model.results;

    % 2.a) PCA on model parameters --> not done (useless)
    
%     level_var = 5;
%     warning('off','all');
%     [coeff_h, score,~,~,explained_F,~] = pca(Prior_h);
%     nb_PC = find(explained_F>level_var,1,'last');
%     PCA_h = score(:,1:nb_PC);
%     coeff_h = coeff_h(:,1:nb_PC);
%     dimh = nb_PC;
%     figure(100);
%     nb_layer = handles.Models.nbLayers;
%     for i = 1 : (4*nb_layer-1),
%        c{i} = ['M^c_' num2str(i)];
%        if i < nb_layer,
%            leg{i} = ['e_' num2str(i)];
%        elseif i< 2*nb_layer,
%            leg{i} = ['V_{P,' num2str(i-nb_layer+1) '}'];
%        elseif i < 3*nb_layer,
%            leg{i} = ['V_{S,' num2str(1+i-nb_layer*2) '}'];
%        else
%            leg{i} = ['\rho_{' num2str(1+i-nb_layer*3) '}'];
%        end
%     end   
%     bar((abs(coeff_h)'./repmat(sum(abs(coeff_h)',1),size(coeff_h,2),1)),'stacked');
%     set(gca,'xticklabel',c);
%     set(gca,'yticklabel',[]);
%     legend(leg);
    
    
    PCA_h = Prior_h - repmat(mean(Prior_h,1),size(Prior_h,1),1); 
    coeff_h = eye(size(Prior_h,2)); 
    dimh = size(PCA_h,2);

    % 2.b) PCA on data

    level_var = 0.1;
    warning('off','all');
    [coeff_d,score,~,~,explained_F,~]=pca(Prior_d);
    warning('on','all');
    nb_PC = max([find(explained_F>level_var,1,'last'),dimh]);% To be sure that there is at least more or equal d than h dim
    PCA_d = score(:,1:nb_PC);
    dimd = nb_PC;
    
    dobs_f = (d_real_obs'-(mean(Prior_d)))*coeff_d(:,1:nb_PC);


    % 2.c) CCA

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
            subplots(jkl,:) = 2*(jkl-1)+1:2*(jkl-1)+2;%2*(jkl-1)+2:2*(jkl-1)+3;
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
    figure(101);
    nb_layer = handles.Models.nbLayers;
    for i = 1 : (4*nb_layer-1),
       c{i} = ['M^c_' num2str(i)];
       if i < nb_layer,
           leg{i} = ['e_' num2str(i)];
       elseif i< 2*nb_layer,
           leg{i} = ['V_{P,' num2str(i-nb_layer+1) '}'];
       elseif i < 3*nb_layer,
           leg{i} = ['V_{S,' num2str(1+i-nb_layer*2) '}'];
       else
           leg{i} = ['\rho_{' num2str(1+i-nb_layer*3) '}'];
       end
    end   
    bar((abs(B)'./repmat(sum(abs(B)',1),size(B,1),1)),'stacked');
    set(gca,'xticklabel',c);
    set(gca,'yticklabel',[]);
    legend(leg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------- End: graph display ---------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% 3) Kernel density 

    disc = 10000;% Number of points in the pdf
    values = -10 : 20/(disc-1) : 10;
    for i = 1 : size(Dc,2)
        waitbar2a(i/(size(Dc,2)*2),handles.uipanel_waitbar,sprintf('Constructing posterior \n in reduced space . . .'));
        fi(:,i) = KernelDensity([Dc(:,i) Hc(:,i)],dobs_c(i),values,[0.01 0.01]); 
    end

    cdf = fi.*0;
    for i = 2 : disc,
        cdf(i,:) = trapz(values(1:i),fi(1:i,:));
    end

    %% 4) Sampling and back-transformation

    nb_sample = 1000;
    nb_layer = handles.Models.nbLayers;
    CCAi = zeros(nb_sample,size(fi,2));
    HpostCoeff = zeros(nb_sample, size(fi,2));
    Hpost = zeros(nb_sample, size(Prior_h,2));
    nb_trial = 0;
    i = 0;
    % For the rejection sampler:
    parameters = get(handles.uitable_prior,'Data');
    parameters = parameters(1:nb_layer,:);
    param_min = parameters(:,[1 3 7 9]);
    param_max = parameters(:,[2 4 8 10]);
    if ~isnan(parameters(1,5)),
        nu_min = parameters(:,5);
        nu_max = parameters(:,6);
    end
    tmp = 1;
    for k = 1 : size(Prior_h,2),
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
        Hpost(i,:) = HpostCoeff(i,:)*pinv(coeff_h)+repmat(mean(Prior_h,1)',1,1)';% Return to real space
        if ~isnan(parameters(1,5)) && handles.Models.increasing,% Nu and increasing
            Vp = Hpost(i,nb_layer:2*nb_layer-1);
            Vs = Hpost(i,2*nb_layer:3*nb_layer-1);
            rho = Hpost(i,3*nb_layer:end);
            Nu = (Vp.^2-(Vs.^2).*2)./((Vp.^2-Vs.^2).*2);
            if not(isempty(find(Hpost(i,:)<0, 1))) || not(isempty(find(Hpost(i,:)<Pmin,1))) || not(isempty(find(Hpost(i,:)>Pmax,1))) || not(isempty(find(Nu<nu_min,1))) || not(isempty(find(Nu>nu_max,1))) || not(issorted(Vp)) || not(issorted(Vs)) || not(issorted(rho)),
                i = i - 1;
                nb_trial = nb_trial + 1;
                if nb_trial >= 10000,
                    warning('Unable to sample correctly! \n \t The data are outside the physical boundaries of the domain!');
                    break
                end
            else
                nb_trial = 0;
            end
        elseif ~isnan(parameters(1,5)),% Nu but not increasing
            Vp = Hpost(i,nb_layer:2*nb_layer-1);
            Vs = Hpost(i,2*nb_layer:3*nb_layer-1);
            Nu = (Vp.^2-(Vs.^2).*2)./((Vp.^2-Vs.^2).*2);
            if not(isempty(find(Hpost(i,:)<0, 1))) || not(isempty(find(Hpost(i,:)<Pmin,1))) || not(isempty(find(Hpost(i,:)>Pmax,1))) || not(isempty(find(Nu<nu_min,1))) || not(isempty(find(Nu>nu_max,1))),
                i = i - 1;
                nb_trial = nb_trial + 1;
                if nb_trial >= 100000,
                    warning('Unable to sample correctly! \n \t The data are outside the physical boundaries of the domain!');
                    break
                end
            end
        elseif handles.Models.increasing,% No nu and increasing
            Vp = Hpost(i,nb_layer:2*nb_layer-1);
            Vs = Hpost(i,2*nb_layer:3*nb_layer-1);
            rho = Hpost(i,3*nb_layer:end);
            if not(isempty(find(Hpost(i,:)<0, 1))) || not(isempty(find(Hpost(i,:)<Pmin,1))) || not(isempty(find(Hpost(i,:)>Pmax,1))) || not(issorted(Vp)) || not(issorted(Vs)) || not(issorted(rho)),
                i = i - 1;
                nb_trial = nb_trial + 1;
                if nb_trial >= 100000,
                    warning('Unable to sample correctly! \n \t The data are outside the physical boundaries of the domain!');
                    break
                end
            end
        else % Only prior Vp Vs and rho
            if not(isempty(find(Hpost(i,:)<0, 1))) || not(isempty(find(Hpost(i,:)<Pmin,1))) || not(isempty(find(Hpost(i,:)>Pmax,1))),
                i = i - 1;
                nb_trial = nb_trial + 1;
                if nb_trial >= 100000,
                    warning('Unable to sample correctly! \n \t The data are outside the physical boundaries of the domain!');
                    break
                end
            end
        end
    end

    models = struct('thick',Hpost(:,1:nb_layer-1),'Vp',Hpost(:,nb_layer:2*nb_layer-1),...
        'Vs',Hpost(:,2*nb_layer:3*nb_layer-1),'rho',Hpost(:,3*nb_layer:end));

    handles.Solution.model = models;
    handles.Solution.N = nb_sample;
    handles.Solution.nbLayers = handles.Models.nbLayers;

    waitbar2a(1,handles.uipanel_waitbar,'Finalizing!');
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
subplot(nb_layer,8, 8*nb_layer-1:8*nb_layer);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
position = get(gca,'Position');
for j = 1 : nb_layer,
    subplot(nb_layer,8,(j-1)*8+1:(j-1)*8+2);% Water content
    if j == 1,
        title('V_p [m/s]','FontSize',16);
    end
    hold on;
    histogram(models_prior.Vp(:,j),10,'Normalization','pdf');
    histogram(models.Vp(:,j),'Normalization','pdf');
    ylabel(['Layer ' num2str(j)],'FontSize',14);
    subplot(nb_layer,8,(j-1)*8+3:(j-1)*8+4);% T2
    if j == 1,
        title('V_s [m/s]','FontSize',16);
    end
    hold on;
    histogram(models_prior.Vs(:,j),10,'Normalization','pdf');
    histogram(models.Vs(:,j),'Normalization','pdf');
    subplot(nb_layer,8,(j-1)*8+5:(j-1)*8+6);
    if j == 1,
        title('\rho [kg/m^3]','FontSize',16);
    end
    hold on;
    histogram(models_prior.rho(:,j),10,'Normalization','pdf');
    histogram(models.rho(:,j),'Normalization','pdf');
    subplot(nb_layer,8,(j-1)*8+7:(j-1)*8+8);
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

param = [models.Vp models.Vs models.rho models.thick];
param_prior = [models_prior.Vp models_prior.Vs, models_prior.rho, models_prior.thick];
nb_param = 4;
param_pre = {'V_{P,','V_{S,','\rho_{','e_{'};
units = {'[m/s]','[m/s]','[kg/m^3]','[m]'};
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
addpath([pwd '/SurfaceWave']);

%% 1) Forward modelling:
% models = handles.Solution.model;
% lambda_t = zeros(length(models.thick),length(handles.Data.Curve.lambda));
% c_t = lambda_t;

% if isempty(gcp('nocreate')),
%     for j = 1 : 1000,
%         [c_t(j,:),lambda_t(j,:)] = MASWaves_theoretical_dispersion_curve_ModHadrien(10 : 1 : 1500,handles.Data.Curve.lambda,handles.Solution.model.thick(j,:),handles.Solution.model.Vp(j,:),handles.Solution.model.Vs(j,:),handles.Solution.model.rho(j,:),handles.Models.nbLayers-1);
%     end
% else
%     thick = models.thick;
%     Vp = models.Vp;
%     Vs = models.Vs;
%     rho = models.rho;
%     nbLayers = handles.Models.nbLayers;
%     lambda = handles.Data.Curve.lambda;
%     parfor j = 1 : 1000,
%         [c_t(j,:),~] = MASWaves_theoretical_dispersion_curve_ModHadrien(10 : 1 : 1500,lambda,thick(j,:),Vp(j,:),Vs(j,:),rho(j,:),nbLayers-1);
%     end
% end
% clear thick Vp Vs rho nbLayers lambda;
% handles.Solution.model.results = c_t;

[handles.Solution,out_error] = gpdcCall_bisMEX(handles.Solution,handles.Data.Curve.f);
handles.Solution.N = handles.Solution.N - sum(out_error);
handles.Solution.model.thick = handles.Solution.model.thick(~out_error,:);
handles.Solution.model.Vp = handles.Solution.model.Vp(~out_error,:);
handles.Solution.model.Vs = handles.Solution.model.Vs(~out_error,:);
handles.Solution.model.rho = handles.Solution.model.rho(~out_error,:);
handles.Solution.model.results = handles.Solution.model.results(~out_error,:);

guidata(hObject, handles);


data_true = handles.Data.Curve.ct';
data_post = handles.Solution.model.results;

models = handles.Solution.model;

%% 2) Computing RMS
RMS = rms(data_post-repmat(data_true,size(data_post,1),1),2);

[RMS,I] = sort(RMS,'descend');

RMS_fact = (RMS)./(max(RMS)+0.005);

HpostCoeff = [models.thick];
HpostCoeff = HpostCoeff(I,:);
%% 3) Presenting the model in 1D profil

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
        if Q_act == length(RMS_Interval),
            Q_act = Q_act-1;
        end
    end
end
RMS_fact_map(1) = nb_map;
% End addition!
% RMS_fact_map = ceil((RMS_fact).*(length(map)/(max((RMS_fact))-min(RMS_fact))));
% RMS_fact_map = RMS_fact_map - min(RMS_fact_map);
% RMS_fact_map(RMS_fact_map == 0) = 1;

figure('Units','normalized','OuterPosition',[0.1 0.1 0.5 0.8]);
subplot(10,3,[1:3:7*3-1]);% V_p
hold on;
for i = 1 : 10 : size(HpostCoeff,1),
    stairs([models.Vp(I(i),1) models.Vp(I(i),:)],[0 cumsum(models.thick(I(i),1:end)) 1000],'Color',map(RMS_fact_map(i),:),'LineWidth',(1-RMS_fact(i))*3);
end
% Adding the true model if existing:
if isfield(handles.Data,'true');
    stairs([handles.Data.true.Vp(1) handles.Data.true.Vp],[0 cumsum(handles.Data.true.thick) 1000],':w','LineWidth',2);
end
set (gca,'Ydir','reverse');
xlabel('V_p [m/s]');
xlim([0 5000]);
ylim([0 50]);
ylabel('Depth [m]');
hold off;
grid on;
subplot(10,3,[2:3:7*3]);% V_s
hold on;
for i = 1 : 10 : size(HpostCoeff,1),
    stairs([models.Vs(I(i),1) models.Vs(I(i),:)],[0 cumsum(models.thick(I(i),1:end)) 1000],'Color',map(RMS_fact_map(i),:),'LineWidth',(1-RMS_fact(i))*3);
end
% Adding the true model if existing:
if isfield(handles.Data,'true');
    stairs([handles.Data.true.Vs(1) handles.Data.true.Vs],[0 cumsum(handles.Data.true.thick) 1000],':w','LineWidth',2);
end
set (gca,'Ydir','reverse');
xlabel('V_s [m/s]');
xlim([0 5000]);
ylim([0 50]);
ylabel('Depth [m]');
hold off;
grid on;
subplot(10,3,[3:3:7*3+1]);% Rho
hold on;
for i = 1 : 10 : size(HpostCoeff,1),
    stairs([models.rho(I(i),1) models.rho(I(i),:)],[0 cumsum(models.thick(I(i),1:end)) 1000],'Color',map(RMS_fact_map(i),:),'LineWidth',(1-RMS_fact(i))*3);
end
% Adding the true model if existing:
if isfield(handles.Data,'true'),
    stairs([handles.Data.true.rho(1) handles.Data.true.rho],[0 cumsum(handles.Data.true.thick) 1000],':w','LineWidth',2);
end
set (gca,'Ydir','reverse');
xlabel('\rho [kg/m^3]');
xlim([0 3000]);
ylim([0 50]);
ylabel('Depth [m]');
hold off;
grid on;

subplot(10,3,28:30);
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
xlabel('RMS error (on phase velocity) [m/s]');


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MASW.Data = handles.Data;
MASW.Models = handles.Models;
MASW.Solution = handles.Solution;
save([handles.Files.path '\Results_MASW.maswbel'],'MASW','-v7.3');
[file,path] = uiputfile('*.maswbel');
resp = movefile([handles.Files.path '\Results_MASW.maswbel'],[path file]);
if resp,
    fprintf('File save at : \n \t %s\\ %s \n',path,file);
else
    fprintf('Error while saving the file! Try again . . .');
end


% --- Executes on button press in checkbox11.
function checkbox11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Models.increasing = get(hObject,'Value');
guidata(hObject, handles);
