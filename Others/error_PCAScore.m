function [isNoisy, CovNoise, elapsedTime] = error_PCAScore(results, noise, path_file)
% Analysis of the effect of noise on thescores of PCA for a given noise level
%
%   Hadrien MICHEL (hadrien.michel[at]student.uliege.be): february 2018
%
%   1) Compute PCA scores for noise-free FID
%   2) Adding noise to the FID (assuming gaussian noise)
%   3) Compute PCA scoes for noisy FID
%   4) Analysisng changes in PCA scores (has the noise some influence on
%   the process?) using graphs
%
%   The function requieres as input:
%       - results: a n-by-m-by-p matrix containing the FID computed for all
%           models (n models). n is the number of models, m is the nmber of
%           pulse moment x 2 (fisrt half: real part, second half: imainary 
%           part) and p is the number of time steps for the acquisition
%       - noise: the level of noise applied to the signals (in V)
%       - path_file: the path to a file directory which has large space
%           available. I STRONGLY recommend to use a SSD for this, due to
%           numerous reading/writing operations that could damage a HDD.
%           USB keys are OK (prefer USB 3.0 or better)
%
%   The function will write in the referenced directory (in no directory is
%   referenced, in current directory).
%
% _____________


if nargin < 3
    path_file = '';
end

% Computing PCA scores for noise-free FID
%% Changes for pendulum
% Prior_D = [];
% for i = 1 : size(results,2)
%     Prior_D = [Prior_D squeeze(results(:,i,:))];
% end

Prior_D = results;
%% End changes

clearvars models;

% Computing memory needs for the whole dataset
memoryNeeds = 0;
memoryNeeds = memoryNeeds + 2 * size(Prior_D,1)*999*8; % for the score computation (noise-free + noisy)
memoryNeeds = memoryNeeds + size(Prior_D,1)*999*999*8; % for the covariance matrix of the differences
memoryNeeds = memoryNeeds + size(Prior_D,1)*999*8; % for the score differences
memoryNeeds = memoryNeeds + 1000 * 8; % for safety (all other small variables)

[uV, ~] = memory;

memoryAvailable = uV.MemAvailableAllArrays;
memoryOneArray = uV.MaxPossibleArrayBytes;

limit_exp_mem = (memoryNeeds*2)/memoryAvailable; % In order to keep at least half the memory free to be able to work
limit_exp_array = (size(Prior_D,1)*999*999*8)/memoryOneArray;

%nb_exp = ceil(max([limit_exp_mem, limit_exp_array])/2)*2;
nb_exp = 1;

warning ('off','all');
[~,score_D,~,~,~,~]=pca(Prior_D);
size_explored = size(score_D,1)/nb_exp;
%score_D_mod = zeros(size_explored,size(score_D,2));
%% Changed for pendulum
%Cov_df = zeros(size_explored,size(score_D,2)-1,size(score_D,2)-1);
%diff_score = zeros(size_explored,size(score_D,2)-1);
Cov_df = zeros(size_explored,size(score_D,2),size(score_D,2));
diff_score = zeros(size_explored,size(score_D,2));
%% End changes

tic
w_bar = waitbar(0,'Please wait... This may take some time');
s = clock;

% for j = 1 : nb_exp,
%     name_file = ['File_PCA_error_noise_' num2str(noise*1e9) 'nV_' num2str(j) '.mat'];
    
    min_explored = 1;%(j-1)*size_explored+1;
    max_explored = size_explored;%min_explored + size_explored-1;
    
    %Computing the scores for noisy data
    for i = min_explored : max_explored
        if mod(i,2)==0
            est_time_tot = etime(clock,s)/i*(size(Prior_D,1))-etime(clock,s);
            seconds = round(mod(est_time_tot,60));
            minutes = round(mod(est_time_tot,60*60)/60);
            hour = round(mod(est_time_tot,60*60*60)/(60*60));
            waitbar(i/size(Prior_D,1),w_bar,sprintf('Please wait... Model nb. %d \n Estimated remaining time : %02u:%02u:%02u', i, hour, minutes, seconds));
        end
        data_real = Prior_D(i,:)+randn(1,size(Prior_D,2))*noise;
        data_synthe = Prior_D;
        data_synthe(i,:) = [];
        %data_(i-min_explored+1,:)=data_real(i-min_explored+1,:)+randn(1,size(Prior_D,2))*noise;
        %[coeff,score_D_tmp,~,~,~,~]=pca(data_synthe);
        [coeff,~,~,~,~,~]=pca(data_synthe);
        score_D_tmp = (data_real - mean(data_synthe,1))*coeff;
        %score_D_mod(i,:) = score_D_tmp(i,:);
        %% Changed for Pendulum
        %Cov_df(i,:,:) = cov([score_D(i,1:end-1);score_D_tmp]);
        Cov_df(i,:,:) = cov([score_D(i,:);score_D_tmp]);
        %diff_score(i,:) = score_D(i,1:end-1)-score_D_tmp;
        diff_score(i,:) = score_D(i,:)-score_D_tmp;
    end
    %save([path_file name_file],'Cov_df','diff_score','min_explored','max_explored','-v7.3');
    
% end

delete(w_bar);

clearvars -except score_D score_D_mod Cov_df diff_score noise path_file nb_exp;

% Rebuilding the big matrix :
%name_big_data = ['File_tot.mat'];
%m = matfile([path_file name_big_data],'Writable',true);
% m.Cov_df = [];
% m.diff_score = [];
% m.score_D = score_D;
% clearvars score_D
% for j = 1 : nb_exp,
%     fprintf('Reading/writing from data file nb. %d. \n',j);
%     name_file = ['File_PCA_error_noise_' num2str(noise*1e9) 'nV_' num2str(j) '.mat'];
%     m_small = matfile([path_file name_file]);
%     m.Cov_df = [m.Cov_df; m_small.Cov_df];
%     m.diff_score = [m.diff_score; m_small.diff_score];
% end

% clearvars -except name_big_data path_file;
% m = matfile([path_file name_big_data]);

Cdf_mean = squeeze(mean(Cov_df,1));
Cov_dif = cov(diff_score);
dif_mean = mean(diff_score,1);
%% Change Pendulum
%dif_rel = diff_score./score_D(1:size(diff_score,1),1:end-1)*100;
dif_rel = diff_score./score_D(1:size(diff_score,1),:)*100;
%% End changes
dif_rel_mean = mean(dif_rel,1);

clearvars -except Cdf_mean Cov_dif dif_mean dif_rel_mean

figure
plot(1:size(dif_mean,2),dif_mean,'ko');
xlabel('Component','Fontsize',16);ylabel('Mean PCA score difference','Fontsize',16);
set(gca,'FontSize',16)
axis([0 100.5 min(dif_mean) max(dif_mean)]);

figure
plot(1:size(dif_mean,2),dif_rel_mean,'ko');
xlabel('Component','Fontsize',16);ylabel('Mean PCA score relative error','Fontsize',16);
set(gca,'FontSize',16)
axis([0 100.5 min(dif_rel_mean) max(dif_rel_mean)]);

figure
image(Cdf_mean,'CDataMapping','scaled')
% map = colormap(jet(20));
% tmp = ((0:size(map,1)-1)./(size(map,1)-1)) .* (max(max(Cdf_mean)-min(min(Cdf_mean))));
% tmp = tmp+min(min(Cdf_mean));
% idx = find(abs(tmp) == min(abs(tmp)),1,'last');
% if isempty(idx) ~= 0,
%     map(1,:) = [1 1 1];
% else
%     map(idx,:) = [1 1 1];
% end
map = colormap(jet);
map(1,:) = 1;
colormap(map);
c = colorbar('eastoutside');
axis([0.5 50.5 0.5 50.5])
xlabel('Component j','Fontsize',16);ylabel('Component i','Fontsize',16);
title('Mean Covariance matrix C_{df,ij}','FontSize',16)
set(gca,'FontSize',16)

figure
image(Cov_dif,'CDataMapping','scaled')
colormap(map);
h=colorbar('eastoutside');
axis([0.5 50.5 0.5 50.5])
xlabel('Component j','Fontsize',16);ylabel('Component i','Fontsize',16);
title('Covariance matrix of score difference C_{df,ij}','FontSize',16)
set(gca,'FontSize',16)

warning ('on','all');

isNoisy = 0;
CovNoise = Cdf_mean;

elapsedTime = toc;

end
