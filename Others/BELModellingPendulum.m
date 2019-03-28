close all
clear all
clc
path_full = pwd;
pathBEL = extractBefore(path_full,'\Others');
addpath([pathBEL '\GlobalFunctions']);
rng(0);

%% Producing the data :

l_true = 3;
h_true = 7.5;
M_true = 40;
time = (10:1:120)';
Y_true = ForwardPendulum(time,[l_true h_true M_true]);
% Adding noise to the data:
    Y_true = Y_true + normrnd(0,0.1,size(Y_true,1),size(Y_true,2));
% End noise
Data = [time,Y_true];

figure;
plot(Data(:,1),Data(:,2),'k','linewidth',2');
xlabel('Time [sec]');
ylabel('Y [m]');

%% Creating and sampling the prior model space:

l = [1 9];
h = [1 9];
M = [0 50];

type = 2; % Uniformely distributed
N = 10000; % Number of models
parameters = [[l; 0 0], [h; M], [0 0; 0 0], [0 0; 0 0]];
nb_layer = 2;
nb_param = 2;

models = ModelGenerator(type, N, parameters, nb_layer);

Models.model.l = models.thick;
Models.model.h = models.param2(:,1);
Models.model.M = models.param2(:,2);

Models.nbLayers = 2;

clear models

%% Modelling the forward response
tic

w = waitbar(0,{'Computing the forward model . . .','Please wait'});
param = [Models.model.l Models.model.h Models.model.M];
time = Data(:,1);

for j = 1 : N,
    if (mod(j,50)==0),
        waitbar(j/N,w);
    end        
    Y(j,:) = ForwardPendulum(time,param(j,:));
end
Models.model.results = Y;
clear Y;
close(w);

%% Dimension reduction:

d_real_obs = Data(:,2);
Prior_h = param;
Prior_d = Models.model.results;

PCA_h = Prior_h - repmat(mean(Prior_h,1),size(Prior_h,1),1);
coeff_h = eye(size(Prior_h,2));
dimh = size(PCA_h,2);

level_var = 0.1;
warning('off','all');
[coeff_d,score,~,~,explained_F,~]=pca(Prior_d);
warning('on','all');
nb_PC = max([find(explained_F>level_var,1,'last'),dimh]);% To be sure that there is at least more or equal d than h dim
PCA_d = score(:,1:nb_PC);
dimd = nb_PC;

dobs_f = (d_real_obs'-(mean(Prior_d)))*coeff_d(:,1:nb_PC);

[A,B,~,Dc,Hc] = canoncorr(PCA_d,PCA_h);
dobs_c = (dobs_f-mean(PCA_d))*A;

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
figure;
for i = 1 : (3*nb_layer-1),
   c{i} = ['M^c_' num2str(i)];
end   
bar((abs(B)'./repmat(sum(abs(B)',1),size(B,1),1)),'stacked');
set(gca,'xticklabel',c);
set(gca,'yticklabel',[]);
legend('l','h','M');

%% Kernel density estimator

disc = 10000;% Number of points in the pdf
values = -10 : 20/(disc-1) : 10;
for i = 1 : size(Dc,2)
    fi(:,i) = KernelDensity([Dc(:,i) Hc(:,i)],dobs_c(i),values,[0.1 0.1]);
end

cdf = fi.*0;
for i = 2 : disc,
    cdf(i,:) = trapz(values(1:i),fi(1:i,:));
end

%% Sampling and back-transformation

nb_sample = 10000;
CCAi = zeros(nb_sample,size(fi,2));
HpostCoeff = zeros(nb_sample, size(fi,2));
nb_trial = 0;

Pmin = [min(l) min(h) min(M)];
Pmax = [max(l) max(h) max(M)];
i = 0;
while i < nb_sample,
    i = i + 1;
    for j = 1 : size(cdf,2),
        ni = rand(1);
        idx = find(cdf(:,j)<=ni,1,'last');
        CCAi(i,j) = values(idx);
    end
    HpostCoeff(i,:) = CCAi(i,:)*pinv(B)+repmat(mean(PCA_h,1)',1,1)';% Return to PCA space
    HpostCoeff(i,:) = HpostCoeff(i,:)*pinv(coeff_h)+repmat(mean(Prior_h,1)',1,1)';% Return to real space
    if not(isempty(find(HpostCoeff(i,:)<Pmin,1))) || not(isempty(find(HpostCoeff(i,:)>Pmax,1))) || not(isempty(find(HpostCoeff(i,1)+HpostCoeff(i,2)<10,1))),
        i = i - 1;
        nb_trial = nb_trial + 1;
        if nb_trial >= 100000,
            warning('Unable to sample correctly! \n \t The data are outside the physical boundaries of the domain!');
            break
        end
    end

end

models = struct('l',HpostCoeff(:,1),'h',HpostCoeff(:,2),'M',HpostCoeff(:,3));

Solution.model = models;

toc

%% Displaying models:

models = Solution.model;
models_prior = Models.model;
nb_layer = Models.nbLayers;

param = [models.l, models.h models.M];
param_true = [l_true h_true M_true];
param_prior = [models_prior.l, models_prior.h models_prior.M];
param_names = {'l [m]','h [m]','M [kg]'};

nb_param = 3;
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
            histogram(param_prior(:,j),10,'DisplayStyle','stairs','Normalization','pdf','LineWidth',2.5);
            histogram(param(:,j),'DisplayStyle','stairs','Normalization','pdf','LineWidth',2.5);
            yl = ylim;
            plot([param_true(j) param_true(j)],yl,':k','LineWidth',2.5);
            if i == nb_param,
                xlabel(param_names{j});
                w = legend('Prior','Posterior','Real');
                w.Position = position;
                w.FontSize = 12;
            elseif j == 1,
                ylabel(param_names{j});
            end
        end
    end
end

%% Displaying data:


w = waitbar(0,{'Computing the forward model for the solutions. . .','Please wait'});
param = [Solution.model.l Solution.model.h Solution.model.M];
time = Data(:,1);

for j = 1 : nb_sample,
    if (mod(j,50)==0),
        waitbar(j/nb_sample,w);
    end        
    Y(j,:) = ForwardPendulum(time,param(j,:));
end
Solution.model.results = Y;
clear Y;
close(w);

%%
figure;
hold on;
Y_prior = Models.model.results;
Y_post = Solution.model.results;

% Adding prior data:
pas = 100;
for i = 1 : pas : size(Y_prior,1),
    plot(time,Y_prior(i,:),'color',[0.5 0.5 0.5],'linewidth',0.2);
end

% Adding posterior data:
for i = 1 : pas : size(Y_post,1),
    plot(time,Y_post(i,:),'b','linewidth',0.2);
end

% Adding true data:
plot(Data(:,1),Data(:,2),'r','linewidth',1);

xlabel('Time [sec]');
ylabel('Y [m]');

RMS = rms(Y_post-repmat(Data(:,2)',size(Y_post,1),1),2);

%% RMS plot with the position of departure:
X_prior = tan(acos((10-Models.model.h)./Models.model.l)).*(10-Models.model.h);
X_soluce = tan(acos((10-Solution.model.h)./Solution.model.l)).*(10-Solution.model.h);
X_true = tan(acos((10-h_true)./l_true)).*(10-h_true);

XY = [X_soluce Solution.model.h];
XY_prior = [X_prior Models.model.h];

[RMS,I] = sort(RMS,'descend');

RMS_fact = (RMS)./(max(RMS)+0.005);
XY = XY(I,:);
map = colormap(jet);
% Adding feature for the same number of models per share of color
nb_map = length(map);
Q_used = (0 : nb_map-1)./(nb_map-1);
thresh = max(RMS);
RMS_Interval = quantile(RMS(RMS<thresh),Q_used);
RMS_fact_map = ones(1,length(RMS));
Q_act = 1;
i = length(RMS);
while RMS(i-1)<thresh,
    RMS_fact_map(i) = Q_act;
    if RMS(i-1) >= RMS_Interval(Q_act),
        Q_act = Q_act + 1;
    end
    i = i - 1;
end
RMS_fact_map(1:i-1) = nb_map;

figure('Units','centimeters','OuterPosition',[1.5 1.5 25 15]);
subplot(10,2,1:2:15);
hold on;
scatter(XY_prior(:,1),XY_prior(:,2),10,[0.5 0.5 0.5],'filled');
sz = ones(length(XY),1);
c = ones(length(XY),3);
for i = 1 : length(XY),
    sz(i) = (1-RMS_fact(i))*20;
    c(i,:) = map(RMS_fact_map(i),:);
end
scatter(XY(:,1),XY(:,2),sz,c,'filled');
s = scatter(X_true,h_true,50);
axis equal;
axis tight;
s.LineWidth = 2;
s.MarkerEdgeColor = 'w';
s.MarkerFaceColor = 'none';
xlim([0 10]);
ylim([0 10]);
xlabel('X [m]');
ylabel('Y [m]');
title('RMS error for the initial position of the bob');
%set(gca,'fontsize', 12);

subplot(10,2,19);
tmp = min(RMS) : (max(RMS)-min(RMS))/999 : max(RMS);
for i = 1 : length(tmp),
    tmp(i) = find(RMS_Interval<=tmp(i),1,'last');
end
imagesc(RMS_Interval,0:1,[(tmp)' (tmp)']');
% End addition
% imagesc([min(RMS) : (max(RMS)-min(RMS))/(size(map,1)-1) : max(RMS)],0:1,[(1:size(map,1))' (1:size(map,1))']');
colormap(map);
set(gca,'YTick',[]);
xlabel('RMS error on the data [m]');
%set(gca,'fontsize', 12);

%% Probability plot for the position of departure:
step = 0.1;
x = 0:step:10;
y = 0:step:10;
proba = zeros(length(x)-1,length(y)-1);
mean_rms = proba;
proba_prior = proba;

area = step^2;

for i = 1:length(x)-1
    for j = 1:length(y)-1
        index1 = XY>repmat([x(i) y(j)],size(XY,1),1);
        index2 = XY<repmat([x(i+1) y(j+1)],size(XY,1),1);
        index = index1&index2;
        index = sum(index,2);
        index = index==2;
        proba(i,j) = sum(index)/size(XY,1);
        if proba(i,j)>0
            mean_rms(i,j) = mean(RMS(index));
        else
            mean_rms(i,j) = nan;
        end
        index1 = XY_prior>repmat([x(i) y(j)],size(XY_prior,1),1);
        index2 = XY_prior<repmat([x(i+1) y(j+1)],size(XY_prior,1),1);
        index = index1&index2;
        index = sum(index,2);
        index = index==2;
        proba_prior(i,j) = sum(index)/size(XY_prior,1);
    end
end

proba(proba==0) = nan;
proba_prior(proba_prior==0) = nan;
proba = proba./area;
proba_prior = proba_prior./area;

position = subplot(10,2,20);
position.Visible = 'off';
h2 = subplot(10,2,2:2:16);
imagesc((x(1:end-1)+x(2:end))./2,(y(1:end-1)+y(2:end))./2,proba','alphadata',~isnan(proba'));
myColormap = colormap('jet');
colormap(myColormap);
set(h2,'Ydir','normal')
axis equal
axis tight
col = colorbar;
col.Location = 'southoutside';
col.Position = position.Position;%'southoutside';
col.Label.String = 'Estimated probability [/]';
hold on;
s = scatter(X_true,h_true,50);
s.LineWidth = 2;
s.MarkerEdgeColor = 'w';
s.MarkerFaceColor = 'none';
title('Probability plot for the initial position of the bob');
xlabel('X [m]');
ylabel('Y [m]');

figure;
h2 = axes;
imagesc((x(1:end-1)+x(2:end))./2,(y(1:end-1)+y(2:end))./2,mean_rms','alphadata',~isnan(mean_rms'));
myColormap = flipud(colormap('jet'));
colormap(myColormap);
set(h2,'Ydir','normal')
axis equal
axis tight
colorbar;
hold on;
X_true = tan(acos((10-h_true)./l_true)).*(10-h_true);
s = scatter(X_true,h_true,50);
s.LineWidth = 2;
s.MarkerEdgeColor = 'k';
s.MarkerFaceColor = 'none';
title('RMS error for the data from different initial positions');
xlabel('X [m]');
ylabel('Y [m]');

% %% Analyzing the results:
% 
% clear Y;
% 
% w = waitbar(0,{'Computing the full forward modeling . . .','Please wait'});
% time = 0.01 : 0.01 : max(time);
% param = [Models.model.l Models.model.h Models.model.M];
% for j = 1 : N,
%     if (mod(j,50)==0),
%         waitbar(j/N,w);
%     end        
%     Y(j,:) = ForwardPendulum(time,param(j,:));
% end
% Models.model.full = Y;
% Models.model.time = time;
% 
% clear Y;
% close(w);
% 
% w = waitbar(0,{'Computing the full forward modeling . . .','Please wait'});
% param = [Solution.model.l Solution.model.h Solution.model.M];
% for j = 1 : nb_sample,
%     if (mod(j,50)==0),
%         waitbar(j/nb_sample,w);
%     end        
%     Y(j,:) = ForwardPendulum(time,param(j,:));
% end
% Solution.model.full = Y;
% Solution.model.time = time;
% 
% clear Y;
% close(w);
% 
% Y_true_full = ForwardPendulum(time,[l_true h_true M_true]);
% 
% %% Displaying the full data (unknown to the experiment):
% 
% figure;
% hold on;
% Y_prior = Models.model.full;
% Y_post = Solution.model.full;
% 
% % Adding prior data:
% pas = 100;
% for i = 1 : pas : size(Y_prior,1),
%     plot(time,Y_prior(i,:),'color',[0.5 0.5 0.5],'linewidth',0.2);
% end
% 
% % Adding posterior data:
% pas  = 20;
% for i = 1 : pas : size(Y_post,1),
%     plot(time,Y_post(i,:),'b','linewidth',0.2);
% end
% 
% % Adding true data:
% plot(time,Y_true_full,'r','linewidth',1);
% plot(Data(:,1),Data(:,2),'.k');
% 
% xlabel('Time [sec]');
% ylabel('Y [m]');

rmpath([pathBEL '\GlobalFunctions']);