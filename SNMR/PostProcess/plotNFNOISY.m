% Pre load BEL results and disp
% Load NF dataset
[file, path] = uigetfile('*.mrsbel','Select the noise-free dataset');
load([path file],'-mat');
models_NF = SNMR.Solution.model;
% Load Noisy dataset
[file, path] = uigetfile('*.mrsbel','Select the noisy dataset');
load([path file],'-mat');
models_NOISY = SNMR.Solution.model;
models_prior = SNMR.Models.model;


nb_layer = 2;
param_NF = [models_NF.water models_NF.T2 models_NF.thick];
param_NOISY = [models_NOISY.water models_NOISY.T2 models_NOISY.thick];
param_prior = [models_prior.water models_prior.T2, models_prior.thick];
nb_param = 3;
param_pre = {'W_{','T^*_{2,','e_{'};
units = {'[m^3/m^3]','[ms]','[m]'};
true_values = [5/100 25/100 100 200 5];
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
            scatter(param_NOISY(:,j),param_NOISY(:,i),'.');
            scatter(param_NF(:,j),param_NF(:,i),'.');
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
            histogram(param_NOISY(:,j),'DisplayStyle','stairs','Normalization','pdf');
            histogram(param_NF(:,j),'DisplayStyle','stairs','Normalization','pdf');
            yl = ylim;
            plot([true_values(i) true_values(i)],yl,':r');
            if i == nb_param,
                xlabel(param_names{j});
                w = legend('Prior','Posterior Noisy','Posterior Noise-free');
                w.Position = position;
                w.FontSize = 12;
            elseif j == 1,
                ylabel(param_names{j});
            end
        end
    end
end

%% Adding the graph of distributions:

load('ResultsTx50Rx50_2L_NF_1000m.mrsSNMR','-mat');

data_true = SNMR.Data.d_real_obs;
data_post = SNMR.Solution.model.results;
models = SNMR.Solution.model;
tmp = [];
for i = 1 : size(data_post,2),
    tmp = [tmp squeeze(data_post(:,i,:))];
end
data_post = tmp;
clear tmp;

% Computing RMS
RMS = rms(data_post.*1e9-repmat(data_true.*1e9,size(data_post,1),1),2);

[RMS,I] = sort(RMS,'descend');

RMS_fact = (RMS)./(max(RMS)+0.005);

HpostCoeff = [models.thick, models.water, (models.T2)];
HpostCoeff = HpostCoeff(I,:);
% Presenting the model in 1D profil

map = colormap(jet);
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

figure;
subplot(10,2,[1:2:13]);% Water content
hold on;
for i = 1 : 1 : size(HpostCoeff,1),
    stairs([models.water(I(i),1) models.water(I(i),:)],[models.depth(I(i),:) 100],'Color',map(RMS_fact_map(i),:),'LineWidth',(1-RMS_fact(i))*3);
end
stairs([0.05 0.05 0.25],[0 5 100],':w','LineWidth',2); 
stairs([0.035 0.035 0.1],[0 7.5 100],':k','LineWidth',2);
stairs([0.1 0.1 0.3],[0 2.5 100],':k','LineWidth',2);
set (gca,'Ydir','reverse');
xlabel('Water content (\phi) [/]');
xlim([0 0.4]);
ylim([0 20]);
ylabel('Depth [m]');
grid on;
hold off;
set(gca,'FontSize',12);

subplot(10,2,[2:2:14]);% Relaxation time
hold on;
for i = 1 : 1 : size(HpostCoeff,1),
    stairs([models.T2(I(i),1) models.T2(I(i),:)],[models.depth(I(i),:) 100],'Color',map(RMS_fact_map(i),:),'LineWidth',(1-RMS_fact(i))*3);
end
stairs([100 100 200],[0 5 100],':w','LineWidth',2); 
stairs([5 5 5],[0 7.5 100],':k','LineWidth',2);
stairs([350 350 350],[0 2.5 100],':k','LineWidth',2);
set (gca,'Ydir','reverse');
xlabel('Relaxation time (T_2^*) [ms]');
xlim([0 400]);
ylim([0 20]);
ylabel('Depth [m]');
grid on;
hold off;
set(gca,'FontSize',12);

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
% hold on;
% yl = ylim;
% plot([35 35],yl,'k','linewidth',3);
xlabel('RMS error [nV]');
set(gca,'FontSize',10);