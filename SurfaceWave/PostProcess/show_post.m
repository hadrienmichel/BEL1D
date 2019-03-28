%% 1) Load data
[file, path] = uigetfile('*.maswbel');
load([path file],'-mat');

%% 2) Show models
models = MASW.Solution.model;
models_prior = MASW.Models.model;
nb_layer = MASW.Models.nbLayers;
figure('Units','normalized','OuterPosition',[0.1 0.1 0.5 0.6]);
subplot(nb_layer,8, 8*nb_layer-1:8*nb_layer);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
set (gca,'Visible','off');
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
set(findall(gcf,'-property','FontSize'),'FontSize',16);

%% Correlations
R = corrcoef([models.Vp models.Vs models.rho models.thick])
R_prior =  corrcoef([models_prior.Vp models_prior.Vs models_prior.rho models_prior.thick])

% %%
% figure;
% histogram(sum(models_prior.thick,2),'Normalization','pdf');
% hold on
% histogram(sum(models.thick,2),'Normalization','pdf');
% ylimits = ylim;
% plot([118 118],ylimits,'k','LineWidth',4);
% legend('Prior','Posterior','Borehole');
% xlabel('Depth to the bedrock [m]');
% ylabel('Probability estimation [/]');

%%
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
            set(gca,'xtick',[])
            set(gca,'ytick',[])
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
            set(gca,'xtick',[])
            set(gca,'ytick',[])
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