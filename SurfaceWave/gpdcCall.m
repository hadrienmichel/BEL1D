function [models, ignored] = gpdcCall(models, f_compute)
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

%% New new version (java cmd)
addpath([pwd '\SurfaceWave\matlab-jsystem-master\src']);

%% Old version:
% vect_compute = 1 : 12 : models.N;
% if vect_compute(end) ~= models.N,
%     vect_compute(end+1) = models.N;
% end
% tmp = 2;
% for j = vect_compute(1:end-1)
%     waitbar(j/models.N,H,sprintf('Model %d on %d . . .',j,models.N));
%     % 1) Create the *.model' file (12 at a time)
%     fileID = fopen('tmp.model','w');
%     for k = j : vect_compute(tmp),
%         fprintf(fileID,'%d\n',models.nbLayers);
%         for i = 1 : models.nbLayers-1
%             fprintf(fileID,'%f\t%f\t%f\t%f\n',models.model.thick(k,i), models.model.Vp(k,i), models.model.Vs(k,i), models.model.rho(k,i));
%         end
%         fprintf(fileID,'%f\t%f\t%f\t%f\n',0, models.model.Vp(k,end), models.model.Vs(k,end), models.model.rho(k,end));
%     end
%     fclose(fileID);
%     n_samples = 1000;
%     n_workers = 6;
%     % 2) Running the forward modelling through external command
%     command = ['gpdc -n ' num2str(n_samples) ' -s frequency -min ' num2str(min(f_compute)) ' -max ' num2str(max(f_compute)) ' -jobs ' num2str(n_workers) ' < tmp.model > tmp.disp'];
%     [status, cmdout] = system(command);
%     if status ~=0 && ~isempty(cmdout) && contains(cmdout,'** Warning ** : dispersion curve impossible to obtain for mode 0'),
%         % Detect the number of computed disp:
%         
%         warning('The model number %d to %d were not computed!',j,vect_compute(tmp));
%         models.model.results(j:vect_compute(tmp),:) = nan;
%         ignored(j:vect_compute(tmp)) = true;
%     else
%         % 3) Reading the generated file
%         for i = 1 : length(j : vect_compute(tmp)),
%             R1 = 3*i + (i-1)*n_samples;
%             R2 = R1 + 999;
%             C1 = 0;
%             C2 = 1;
%             data = dlmread('tmp.disp',' ',[R1 C1 R2 C2]);
%             for k = 1 : length(f_compute)
%                 models.model.results(j+i-1,k) = data(find(abs(data(:,1)-f_compute(k))==min(abs(data(:,1)-f_compute(k))),1,'first'),2)^(-1);
%                 if min(data(:,1)-f_compute(k)) > max_error,
%                     max_error = min(abs(data(:,1)-f_compute(k)));
%                 end
%             end
%         end
%     end
%     tmp = tmp + 1;
%     delete('tmp.model');
% end

%% New version - january 2019
remainings = 1;% To models.N;
computed = false;
while ~computed
    waitbar(remainings/models.N,H,sprintf('Model %d on %d . . .',remainings,models.N));
    fileID = fopen('tmp.model','w');
    for i = remainings : models.N
        fprintf(fileID,'%d\n',models.nbLayers);
        for j = 1 : models.nbLayers-1
            fprintf(fileID,'%f\t%f\t%f\t%f\n',models.model.thick(i,j), models.model.Vp(i,j), models.model.Vs(i,j), models.model.rho(i,j));
        end
        fprintf(fileID,'%f\t%f\t%f\t%f\n',0, models.model.Vp(i,end), models.model.Vs(i,end), models.model.rho(i,end));
    end
    fclose(fileID);
    % The tmp.model is completed
    n_samples = 1000;
    n_workers = feature('numcores')-1;
    % 2) Running the forward modelling through external command
    % command = ['gpdc -n ' num2str(n_samples) ' -s frequency -min ' num2str(min(f_compute)) ' -max ' num2str(max(f_compute)) ' -jobs ' num2str(n_workers) ' < tmp.model > tmp.disp'];
    command = ['gpdc -n ' num2str(n_samples) ' -s frequency -min ' num2str(min(f_compute)) ' -max ' num2str(max(f_compute)) ' -jobs ' num2str(n_workers) ' < tmp.model'];
    % [status, cmdout] = system(command);
    [status, cmdout, stderr] = jsystem(command);
    cmdout = splitlines(cmdout);
    n = length(cmdout);
    nbOK = floor(n/(n_samples+3));
    f_computed = linspace(min(f_compute),max(f_compute),n_samples);
    if nbOK~=(models.N-remainings)+1%status~=0 && contains(stderr,'** Warning ** : dispersion curve impossible to obtain for mode 0')        
%         fid = fopen('tmp.disp','r');
%         n = 0;
%         tline = fgetl(fid);
%         while ischar(tline)
%           tline = fgetl(fid);
%           n = n+1;
%         end
%         nbOK = ceil(n/1003);
%         fclose(fid);
        tmp = 1;
        if nbOK~=0
            for i = remainings : remainings+nbOK-1
                R1 = 4 + (tmp-1)*(n_samples+3);
                R2 = R1 + 999;
%                 R1 = 3 + (tmp-1)*(n_samples+3);
%                 R2 = R1 + 999;
%                 C1 = 0;
%                 C2 = 1;
%                 current = [R1 C1 R2 C2];
%                 disp(current);
%                 data = dlmread('tmp.disp',' ',[R1 C1 R2 C2]);
                data = cmdout(R1:R2);% Data is a cell array with each cell corresponding to a char with '%f %f' formatting.
                for k = 1 : length(f_compute)
                    index = find(abs(f_computed-f_compute(k))==min(abs(f_computed-f_compute(k))),1,'first');
                    dataTmp = str2num(data{index});
                    models.model.results(i,k) = dataTmp(2)^(-1);
                    if abs(dataTmp(1)-f_compute(k)) > max_error
                        max_error = (abs(dataTmp(1)-f_compute(k)));
                    end
                end
                tmp = tmp + 1;
            end
            models.model.results(remainings+nbOK,:) = nan;
            ignored(remainings+nbOK) = true;
            remainings = remainings + nbOK+1;% We skip the uncomputable model
        else
            models.model.results(remainings,:) = nan;
            ignored(remainings) = true;
            remainings = remainings + 1;
        end
    else
        tmp = 1;
        for i = remainings : models.N
            R1 = 4 + (tmp-1)*(n_samples+3);
            R2 = R1 + 999;
            data = cmdout(R1:R2);% Data is a cell array with each cell corresponding to a char with '%f %f' formatting.
            for k = 1 : length(f_compute)
                index = find(abs(f_computed-f_compute(k))==min(abs(f_computed-f_compute(k))),1,'first');
                dataTmp = str2num(data{index});
                models.model.results(i,k) = dataTmp(2)^(-1);
                if abs(dataTmp(1)-f_compute(k)) > max_error
                    max_error = (abs(dataTmp(1)-f_compute(k)));
                end
            end
            tmp = tmp + 1;
%             R1 = 3*tmp + (tmp-1)*n_samples;
%             R2 = R1 + 999;
%             C1 = 0;
%             C2 = 1;
%             data = dlmread('tmp.disp',' ',[R1 C1 R2 C2]);
%             for k = 1 : length(f_compute)
%                 models.model.results(i,k) = data(find(abs(data(:,1)-f_compute(k))==min(abs(data(:,1)-f_compute(k))),1,'first'),2)^(-1);
%                 if min(data(:,1)-f_compute(k)) > max_error
%                     max_error = min(abs(data(:,1)-f_compute(k)));
%                 end
%             end
%             tmp = tmp + 1;
        end
        computed = true;
    end
    delete('tmp.model');%,'tmp.disp');
end

fprintf('%d models were not computed (on %d)\nMaximum error on the frequency = %f\n',sum(ignored),models.N,max_error);
toc
close(H);
rmpath([pwd '\SurfaceWave\matlab-jsystem-master\src']);

end