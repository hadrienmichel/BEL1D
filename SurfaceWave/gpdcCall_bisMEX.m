function [models, ignored] = gpdcCall_bisMEX(models, f_compute)
% GPDCCALL is used to perform the forward modelling of the dispersions
% curves for the models described in the models structure.
%
% In order to run, the folder containing the 
%
% The functions takes as input
%       - models: a structure that contains the models informations
%           >> models.nbLayers -> the number of layers in the models
%           >> models.N -> the number of models to simulate
%           >> models.model -> a structure with the models constititions
%                   models.model.thick: a N-by-(nbLayers-1) matrix 
%                   models.model.Vp: a N-by-nbLayers matrix
%                   models.model.Vs: a N-by-nbLayers matrix
%                   models.model.rho: a N-by-nbLayers matrix
%
% The function outputs the solution of the forward model
%       - solution: a structure containing the results
%
tic;
models.model.results = zeros(models.N,length(f_compute));
H = waitbar(0,'Please wait, forward modelling in progress!');
ignored = false(1,models.N);
max_error = 0;

%% New new version (mex file (MuLTI))
if ispc
    addpath([pwd '\SurfaceWave\MuLTI\MuLTI-master\gpdc mex file WINDOWS\']);% Windows systems
elseif isunix
    warning('This version need to be compiled by the user!');
    addpath([pwd '\SurfaceWave\MuLTI\MuLTI-master\gpdc mex file LINUX\']);% UNIX systems: NOT TESTED!!!
else
    error('No version of gpdc in MEX interface is suported for mac users!');
end

if size(f_compute,1) > 1,
    f_compute = f_compute';
end
for i = 1 : models.N,
    waitbar(i/models.N,H,'Please wait . . .');
    try
        out = gpdc([models.model.thick(i,:) 0],models.model.Vp(i,:),models.model.Vs(i,:),models.model.rho(i,:),'fV',f_compute);
        out = out(:,2); % Only intersted in the fundamental mode (Rayleight)
        if ~isempty(find(isnan(out),1)),
            ignored(i) = true;
        else
            models.model.results(i,:) = (out').^(-1);
        end
    catch
        ignored(i) = true;
    end
end

fprintf('%d models were not computed (on %d)\n',sum(ignored),models.N);
toc
close(H);
rmpath([pwd '\SurfaceWave\MuLTI\MuLTI-master\gpdc mex file WINDOWS\']);

end