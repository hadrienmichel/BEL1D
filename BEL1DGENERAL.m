function varargout = BEL1DGENERAL(varargin)
% BEL1DGENERAL is a graphical interface that helps the user to perform
% BEL modeling of the subsurface using WATHEVER TYPE OF DATA YOU WANT. 
%
% The function requiers the use of a forward model according to the 
% formatting of the function 'forwardModel.m'. This later must correspond 
% to the formatting: 
%   function Y = forwardModel(X,param)
%   % Inputs: - X: the coordonates for wich the model must be calculated
%   %         - param: the parameters of the model in a vector form [1 x
%   %                  ((nxm)-1)]
%   %                 o n is the number of layers
%   %                 o m is the number of parameters
%   % Outputs: - Y: the values of the model at coordinates X
%   Y = some function of X and param
%   end
%
% For further help, please contact the author:
%       Name: Hadrien MICHEL
%       University: University of Liège
%       Country: Belgium
%       e-mail: Hadrien.Michel@uliege.be
%
% © 2019 Hadrien MICHEL, University of Liège, University of Ghent and
% F.R.S.-FNRS
%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BEL1DGENERAL_OpeningFcn, ...
                   'gui_OutputFcn',  @BEL1DGENERAL_OutputFcn, ...
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


% --- Executes just before BEL1DGENERAL is made visible.
function BEL1DGENERAL_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)
addpath([pwd '/GlobalFunctions']);
% Choose default command line output for BEL1DGENERAL
addpath(genpath(pwd),'-begin');

handles.output = hObject;

handles.Files = [];% Contains the files informations
handles.Files.path = 'Files path';
handles.Files.name = 'Name File';
handles.Data = [];% Contains the data-sets
handles.Models = [];
handles.Models.nbLayers = 2;
param(6,1) = struct;
for i = 1 : 6,
    param(i).names = {'Thickness [m]','nan','nan','nan'};
    param(i).min = [0 nan nan nan];
    param(i).max = [5 nan nan nan];
end
handles.Models.param = param;
tmp = {param(1).names' param(1).min' param(1).max'};
for i = 1 : 4,
    for j = 1 : 3,
        if j ~= 1,
            data{i,j} = tmp{j}(i);
        else
            data{i,j} = tmp{j}{i};
        end
    end
end
set(handles.uitable_prior,'Data',data);
handles.Models.type = 1;
handles.Models.N = 1000;
handles.Solution = [];

handles.ForwardModel = str2func('forwardModel');

handles.Status = [];
handles.Status.data = false;%data loaded?
handles.Status.models = false;%Model generated?
handles.Status.PFA = false;%PFA performed?

set(handles.uipanel_model,'Visible','off');
set(handles.uipanel_PFA,'Visible','off');

set(handles.axes2,'Visible','off');
set(handles.axes3,'Visible','off');
set(handles.editX,'Visible','off');
set(handles.editY,'Visible','off');
set(handles.text10,'Visible','off');
set(handles.text11,'Visible','off');
handles.Data_Xname = 'X';
handles.Data_Yname = 'Y';

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

% For reproductibility: controle random number generation:
rng(0);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BEL1DGENERAL wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BEL1DGENERAL_OutputFcn(hObject, eventdata, handles)
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
tmp = get(hObject,'String');
[filepath,name,ext] = fileparts(tmp);
handles.Files.path = filepath;
handles.Files.name = [name ext];

handles.Data = load([handles.Files.path handles.Files.name]);% The file MUST contain a 'data' matrix !!!
handles.Data = handles.Data.data;
axes(handles.axes2);

plot(handles.Data(:,1),handles.Data(:,2),'k','linewidth',2');
set(handles.axes2,'Visible','on');
xlabel('X');
ylabel('Y');
set(handles.editX,'Visible','on');
set(handles.editY,'Visible','on');
set(handles.text10,'Visible','on');
set(handles.text11,'Visible','on');


set(handles.uipanel_model,'Visible','on');
handles.Status.data = true;
guidata(hObject, handles);

% --- Executes on button press in pushbutton_file.
function pushbutton_file_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.Files.name ,handles.Files.path] = uigetfile({'*.mat','MATLAB data file (*.mat)'},'Select a data file');
set(handles.edit_data_path,'String',[handles.Files.path handles.Files.name]);
handles.Data = load([handles.Files.path handles.Files.name]);% The file MUST contain a 'data' matrix !!!
handles.Data = handles.Data.data;
axes(handles.axes2);

plot(handles.Data(:,1),handles.Data(:,2),'k','linewidth',2');
set(handles.axes2,'Visible','on');
xlabel('X');
ylabel('Y');
set(handles.editX,'Visible','on');
set(handles.editY,'Visible','on');
set(handles.text10,'Visible','on');
set(handles.text11,'Visible','on');

set(handles.uipanel_model,'Visible','on');
handles.Status.data = true;
guidata(hObject, handles);


% --- Executes on button press in pushbutton_addLayer.
function pushbutton_addLayer_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_addLayer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.Models.nbLayers < 6,
    handles.Models.nbLayers = handles.Models.nbLayers + 1;
    tmp = get(handles.listbox3,'String');
    tmp = char(tmp);
    contents = cellstr(get(handles.listbox3,'String'));% returns listbox3 contents as cell array
    num = char(contents(end));
    num = str2double(num(end));
    added = ['Layer ' num2str(num+1)];
    tmp = char(tmp,added);
    set(handles.listbox3,'String',cellstr(tmp));
    guidata(hObject, handles);
else
    error('Too many layers in the model (max. is 6)!');
end


% --- Executes on button press in pushbutton_removeLayer.
function pushbutton_removeLayer_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_removeLayer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.Models.nbLayers > 2,
    handles.Models.nbLayers = handles.Models.nbLayers - 1;
    tmp = get(handles.listbox3,'String');
    tmp = char(tmp);
    tmp = tmp(1:end-1,:);
    set(handles.listbox3,'String',cellstr(tmp));
    guidata(hObject, handles);
else
    error('Too few layers in the model (min. is 2)!');
end


% --- Executes on selection change in popupmenu_TypeMod.
function popupmenu_TypeMod_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_TypeMod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_TypeMod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_TypeMod
handles.Models.type = get(handles.popupmenu_TypeMod, 'value');
guidata(hObject, handles);



function edit_nbLayer_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nbLayer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nbLayer as text
%        str2double(get(hObject,'String')) returns contents of edit_nbLayer as a double
handles.Models.N = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes on button press in pushbutton_generateModels.
function pushbutton_generateModels_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_generateModels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
addpath([pwd '/GlobalFunctions']);
addpath([pwd '/General']);
% Function for model generation:
% [models] = ModelGenerator(type, N, parameters, nb_layer)

% Catching the number of parameters:
nb_param = sum(not(isnan(handles.Models.param(1).min)));
parameters = zeros(handles.Models.nbLayers,8);
for i = 1 : handles.Models.nbLayers,
    for j = 1 : nb_param,
        parameters(i,(j-1)*2+1) = handles.Models.param(i).min(j);
        parameters(i,(j-1)*2+2) = handles.Models.param(i).max(j);
    end
end
handles.Models.model = ModelGenerator(handles.Models.type, handles.Models.N, parameters, handles.Models.nbLayers);
if nb_param < 4,
    handles.Models.model = rmfield(handles.Models.model,'param4');
end
if nb_param <3,
    handles.Models.model = rmfield(handles.Models.model,'param3');
end
if nb_param < 2,
    handles.Models.model = rmfield(handles.Models.model,'param2');
end
handles.Status.models = true;
if handles.Status.models,
    set(handles.uipanel_PFA,'Visible','on');
end
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Run.
function pushbutton_Run_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
addpath([pwd '/GlobalFunctions']);
addpath([pwd '/General']);
% Main object of the PFA imaging --> runing dimension reductions and kernel
% density
if handles.Status.PFA ~= true,
    set(handles.uipanel_waitbar,'Visible','on');
    drawnow;
    % Computing PFA

    %% 1) Forward modelling:
    
    waitbar2a(0,handles.uipanel_waitbar,'waitbartext',sprintf('Please wait . . . \n Computing forward models'));
    function_used = handles.ForwardModel;
    X_true = handles.Data(:,1);
    models = handles.Models.model;
    Y = zeros(handles.Models.N,length(X_true));
    nb_param = sum(not(isnan(handles.Models.param(1).min)));
    if nb_param <2,
        param = [models.thick];
    elseif nb_param < 3,
        param = [models.thick, models.param2];
    elseif nb_param < 4,
        param = [models.thick, models.param2, models.param3];
    else
        param = [models.thick, models.param2, models.param3, models.param4];
    end
    if isempty(gcp('nocreate')),
        for j = 1 : handles.Models.N,
            if (mod(j,50)==0),
                waitbar2a(j/handles.Models.N,handles.uipanel_waitbar);
            end        
            Y(j,:) = function_used(X_true,param(j,:));
        end
    else
        tested = [1 100 : 100 : handles.Models.N];
        for jklm = 1 : length(tested)-1,
            waitbar2a(tested(jklm)/handles.Models.N,handles.uipanel_waitbar);
            parfor j = tested(jklm): tested(jklm+1),       
                Y(j,:) = function_used(X_true,param(j,:));
            end
        end
    end
    handles.Models.model.results = Y;
    guidata(hObject, handles);

    waitbar2a(0.5,handles.uipanel_waitbar,sprintf('Formatting the datasets . . . '));

    %% 2) Dimension reduction

    waitbar2a(0.5,handles.uipanel_waitbar,sprintf('Reducing dimensionality \n on data. . . '));
    d_real_obs = handles.Data(:,2);

    Prior_h = param;
    Prior_d = handles.Models.model.results;

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
    
    dobs_f = (d_real_obs'-(mean(Prior_d)))*coeff_d(:,1:nb_PC);


    % 3.c) CCA

    waitbar2a(0.5,handles.uipanel_waitbar,sprintf('Canonical correlation \n analysis . . . '));
    [A,B,~,Dc,Hc] = canoncorr(PCA_d,PCA_h);
    dobs_c = (dobs_f-mean(PCA_d))*A;
    
    %       /\          Printing the graph on the influence of each 
    %      /  \         parameter on the different CCA dimensions
    %     /  | \        
    %    /   |  \       COMMENT if unwanted!!!
    %   /    |   \
    %  /     .    \
    % /____________\
    figure('Units','normalized','Position',[0.05 0.05 0.8 0.8])
    subplots = zeros(size(B,1),2);
    for jkl = 1 : size(B,1),
        if jkl <= ceil(size(B,1)/2),
            subplots(jkl,:) = 2*(jkl-1)+1:2*(jkl-1)+2;
        else
            subplots(jkl,:) = 2*(jkl-1)+2:2*(jkl-1)+3;
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
    for i = 1 : (nb_param*nb_layer-1),
       c{i} = ['M^c_' num2str(i)];
       if i < nb_layer,
           leg{i} = ['e_' num2str(i)];
       elseif i< 2*nb_layer,
           leg{i} = ['Param_{2' num2str(i-nb_layer+1) '}'];
       elseif i < 3*nb_layer,
           leg{i} = ['Param_{3,' num2str(1+i-nb_layer*2) '}'];
       else
           leg{i} = ['Param_{4' num2str(1+i-nb_layer*3) '}'];
       end
    end   
    bar((abs(B)'./repmat(sum(abs(B)',1),size(B,1),1)),'stacked');
    set(gca,'xticklabel',c);
    set(gca,'yticklabel',[]);
    legend(leg);
    % end graph influance

    %% 3) Kernel density 

    disc = 10000;% Number of points in the pdf
    values = -10 : 20/(disc-1) : 10;
    for i = 1 : size(Dc,2)
        waitbar2a(i/(size(Dc,2)*2),handles.uipanel_waitbar,sprintf('Constructing posterior \n in reduced space . . .'));
        fi(:,i) = KernelDensity([Dc(:,i) Hc(:,i)],dobs_c(i),values,[0.1 0.1]);
%         [fi(:,i),~] = ksdensity([Dc(:,i) Hc(:,i)],[(ones(disc,1)*dobs_c(i)) (values)'],'Bandwidth',[0.1 0.1]);% Adding parameter for bandwidth
%         %trapz(values,fi(:,i))
%         % The returned pdf is not normalised --> normalisation :
%         fi(:,i) = fi(:,i)./trapz(values,fi(:,i));   
    end

    cdf = fi.*0;
    for i = 2 : disc,
        cdf(i,:) = trapz(values(1:i),fi(1:i,:));
    end

    % Checking that there exists a distribution for all the CCA dimensions
    is_outside = false(1,size(cdf,2));
    for i = 1 : size(cdf,2),
        waitbar2a((i+(size(Dc,2)*2))/(size(Dc,2)*2),handles.uipanel_waitbar);
        if isnan(cdf(1,i)),
            % warning(sprintf('The data you are trying to interpret is not part of the prior model!\n\tResults are probably non-representative.\n\n\tComputing the FULL Kernel density function is long . . . '))% If no information on data, fully spread
            is_outside(i) = true;
            %% New code 30-04-18 (Hadrien MICHEL)
            [X,Y] = meshgrid(-10:0.5:10,values);
            grid_eval = zeros(size(X,1)*size(X,2),2);
            for j = 1 : length(grid_eval),
                grid_eval(j,:) = [X(j),Y(j)];
            end
            [fi_tmp,xi] = ksdensity([Dc(:,i) Hc(:,i)],grid_eval);
            for j = 2 : length(values)-1,
                idx = find(xi(:,2) == values(j));
                fi(j,i) = trapz(xi(idx,1),fi_tmp(idx));                
            end
            fi(1,i) = 0;
            fi(end,i) = 0;
            cdf(1,i) = 0;
            for j = 2 : disc,
                cdf(j,i) = trapz(values(1:j),fi(1:j,i));
            end
        end
    end

    %% 4) Sampling and back-transformation

    % set(handles.text_progress,'String','Sampling and back-transformation . . .');
    nb_sample = 1000;
    nb_layer = handles.Models.nbLayers;
    CCAi = zeros(nb_sample,size(fi,2));
    HpostCoeff = zeros(nb_sample, size(fi,2));
    nb_trial = 0;
    % For the rejection sampler:
    nb_param = sum(not(isnan(handles.Models.param(1).min)));
    parameters = zeros(handles.Models.nbLayers,8);
    for i = 1 : handles.Models.nbLayers,
        for j = 1 : nb_param,
            parameters(i,(j-1)*2+1) = handles.Models.param(i).min(j);
            parameters(i,(j-1)*2+2) = handles.Models.param(i).max(j);
        end
    end
    i = 0;
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
        if not(isempty(find(HpostCoeff(i,:)<Pmin,1))) || not(isempty(find(HpostCoeff(i,:)>Pmax,1))),
            i = i - 1;
            nb_trial = nb_trial + 1;
            if nb_trial >= 100000,
                warning('Unable to sample correctly! \n \t The data are outside the physical boundaries of the domain!');
                break
            end
        end
        
    end

    if nb_param < 2,
        models = struct('thick',HpostCoeff(:,1:nb_layer-1));
    elseif nb_param < 3,
        models = struct('thick',HpostCoeff(:,1:nb_layer-1),'param2',HpostCoeff(:,nb_layer:2*nb_layer-1));
    elseif nb_param < 4,
        models = struct('thick',HpostCoeff(:,1:nb_layer-1),'param2',HpostCoeff(:,nb_layer:2*nb_layer-1),...
            'param3',HpostCoeff(:,2*nb_layer:3*nb_layer-1));
    else
        models = struct('thick',HpostCoeff(:,1:nb_layer-1),'param2',HpostCoeff(:,nb_layer:2*nb_layer-1),...
            'param3',HpostCoeff(:,2*nb_layer:3*nb_layer-1),'param4',HpostCoeff(:,3*nb_layer:end));
    end

    handles.Solution.model = models;

    waitbar2a(1,handles.uipanel_waitbar,'Finalizing!');
    handles.Status.PFA = true;
    guidata(hObject, handles);
end

% set(handles.axes_progressBar,'Visible','off');
% set(handles.text_progress,'Visible','off');
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
nb_param = sum(not(isnan(handles.Models.param(1).min)));
if nb_param <2,
    param = [models.thick];
elseif nb_param < 3,
    param = [models.param2, models.thick];
elseif nb_param < 4,
    param = [models.param2, models.param3, models.thick];
else
    param = [models.param2, models.param3, models.param4, models.thick];
end
if nb_param <2,
    param_prior = [models_prior.thick];
    param_pre = {handles.Models.param(1).names(1)};
elseif nb_param < 3,
    param_prior = [models_prior.param2, models_prior.thick];
    param_pre = {handles.Models.param(1).names(2),handles.Models.param(1).names(1)};
elseif nb_param < 4,
    param_prior = [models_prior.param2, models_prior.param3, models_prior.thick];
    param_pre = {handles.Models.param(1).names(2),handles.Models.param(1).names(3),handles.Models.param(1).names(1)};
else
    param_prior = [models_prior.param2, models_prior.param3, models_prior.param4, models_prior.thick];
    param_pre = {handles.Models.param(1).names(2),handles.Models.param(1).names(3),handles.Models.param(1).names(4),handles.Models.param(1).names(1)};
end

param_names = {};
for i = 1 : nb_param,
    for j = 1 : nb_layer,
        layer = ['(Layer ' num2str(j) ')'];
        param_names{(i-1)*nb_layer+j} = ([param_pre{i} layer]);
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
addpath([pwd '/GlobalFunctions']);
addpath([pwd '/General']);
%% 1) Forward modelling:
models = handles.Solution.model;
function_used = handles.ForwardModel;
X_true = handles.Data(:,1);
nb_param = sum(not(isnan(handles.Models.param(1).min)));
Y = zeros(1000,length(X_true));
if nb_param <2,
    param = [models.thick];
elseif nb_param < 3,
    param = [models.thick, models.param2];
elseif nb_param < 4,
    param = [models.thick, models.param2, models.param3];
else
    param = [models.thick, models.param2, models.param3, models.param4];
end
for j = 1 : 1000,      
    Y(j,:) = function_used(X_true,param(j,:));
end
handles.Solution.model.results = Y;
guidata(hObject, handles);
data_true = handles.Data(:,2)';
data_post = handles.Solution.model.results;

%% 3) Computing RMS
RMS = rms(data_post-repmat(data_true,size(data_post,1),1),2);

[RMS,I] = sort(RMS,'descend');

RMS_fact = (RMS)./(max(RMS)+0.005);

HpostCoeff = [models.thick];
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
if nb_param > 1,
    subplot(10,3,[1:3:7*3-1]);% V_p
    hold on;
    for i = 1 : 10 : size(HpostCoeff,1),
        if ~isnan(RMS(i)),
            stairs([models.param2(I(i),1) models.param2(I(i),:)],[0 cumsum(models.thick(I(i),1:end)) 2000],'Color',map(RMS_fact_map(i),:),'LineWidth',(1-RMS_fact(i))*3);
        end
    end
    set (gca,'Ydir','reverse');
    xlabel(handles.Models.param(1).names(2));
    ylabel('Depth [m]');
    hold off;
    grid on;
end
if nb_param > 2,
    subplot(10,3,[2:3:7*3]);% V_s
    hold on;
    for i = 1 : 10 : size(HpostCoeff,1),
        if ~isnan(RMS(i)),
            stairs([models.param3(I(i),1) models.param3(I(i),:)],[0 cumsum(models.thick(I(i),1:end)) 2000],'Color',map(RMS_fact_map(i),:),'LineWidth',(1-RMS_fact(i))*3);
        end
    end
    set (gca,'Ydir','reverse');
    xlabel(handles.Models.param(1).names(3));
    ylabel('Depth [m]');
    hold off;
    grid on;
end
if nb_param > 3,
    subplot(10,3,[3:3:7*3+1]);% Rho
    hold on;
    for i = 1 : 10 : size(HpostCoeff,1),
        if ~isnan(RMS(i)),
            stairs([models.param4(I(i),1) models.param4(I(i),:)],[0 cumsum(models.thick(I(i),1:end)) 2000],'Color',map(RMS_fact_map(i),:),'LineWidth',(1-RMS_fact(i))*3);
        end
    end
    set (gca,'Ydir','reverse');
    xlabel(handles.Models.param(1).names(4));
    ylabel('Depth [m]');
    hold off;
    grid on;
end

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
xlabel(['RMS error on ' handles.Data_Yname]);


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MBEL.Data = handles.Data;
MBEL.Data_Xname = handles.Data_Xname;
MBEL.Data_Yname = handles.Data_Yname;
MBEL.Models = handles.Models;
MBEL.Solution = handles.Solution;
save([handles.Files.path '\Results_PFA.mpfa'],'MBEL','-v7.3');
[file,path] = uiputfile('*.mbel');
resp = movefile([handles.Files.path '\Results_PFA.mbel'],[path file]);
if resp,
    fprintf('File save at : \n \t %s%s \n',path,file);
else
    fprintf('Error while saving the file! Try again . . . \n');
end


% --- Executes on selection change in listbox3.
function listbox3_Callback(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(hObject,'String'));% returns listbox3 contents as cell array
num = contents{get(hObject,'Value')};% returns selected item from listbox3
num = num(end);
num = str2double(num);

tmp = {handles.Models.param(num).names' handles.Models.param(num).min' handles.Models.param(num).max'};
for i = 1 : 4,
    for j = 1 : 3,
        if j ~= 1,
            data{i,j} = tmp{j}(i);
        else
            data{i,j} = tmp{j}{i};
        end
    end
end
set(handles.uitable_prior,'Data',data);



% --- Executes during object creation, after setting all properties.
function listbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
handles.ForwardModel = get(hObject,'String');
handles.ForwardModel = str2func(handles.ForwardModel);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in uitable_prior.
function uitable_prior_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable_prior (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(handles.listbox3,'String'));% returns listbox3 contents as cell array
num = contents{get(handles.listbox3,'Value')};% returns selected item from listbox3
num = num(end);
num = str2double(num);
if eventdata.Indices(2) == 1,
    for i = 1 : 6,
        handles.Models.param(i).names(eventdata.Indices(1)) = cellstr(eventdata.NewData);
    end
elseif eventdata.Indices(2) == 2,
    handles.Models.param(num).min(eventdata.Indices(1)) = eventdata.NewData;
else
    handles.Models.param(num).max(eventdata.Indices(1)) = eventdata.NewData;
end
guidata(hObject, handles);



function editX_Callback(hObject, eventdata, handles)
% hObject    handle to editX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editX as text
%        str2double(get(hObject,'String')) returns contents of editX as a double
handles.Data_Xname = get(hObject,'String');
guidata(hObject,handles);
axes(handles.axes2)
xlabel(handles.Data_Xname);




% --- Executes during object creation, after setting all properties.
function editX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editY_Callback(hObject, eventdata, handles)
% hObject    handle to editY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editY as text
%        str2double(get(hObject,'String')) returns contents of editY as a double
handles.Data_Yname = get(hObject,'String');
guidata(hObject,handles);
axes(handles.axes2)
ylabel(handles.Data_Yname);

% --- Executes during object creation, after setting all properties.
function editY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
