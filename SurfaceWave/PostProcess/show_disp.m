%% 1) Load data
[file, path] = uigetfile('*.maswbel');
load([path file],'-mat');
data_post = MASW.Solution.model.results;
data_true = MASW.Data.Curve.ct;
models = MASW.Solution.model;
%% 2) Computing RMS
RMS = rms(data_post-repmat(data_true',size(data_post,1),1),2);

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
for i = 1 : 1 : size(HpostCoeff,1),
    stairs([models.Vp(I(i),1) models.Vp(I(i),:)],[0 cumsum(models.thick(I(i),1:end)) 1000],'Color',map(RMS_fact_map(i),:),'LineWidth',(1-RMS_fact(i))*3);
end
% % Adding the true model if existing:
% if isfield(handles.Data,'true');
%     stairs([handles.Data.true.Vp(1) handles.Data.true.Vp],[0 cumsum(handles.Data.true.thick) 1000],':w','LineWidth',2);
% end
set (gca,'Ydir','reverse');
xlabel('V_p [m/s]');
xlim([0 5000]);
ylim([0 500]);
ylabel('Depth [m]');
hold off;
grid on;
subplot(10,3,[2:3:7*3]);% V_s
hold on;
for i = 1 : 1 : size(HpostCoeff,1),
    stairs([models.Vs(I(i),1) models.Vs(I(i),:)],[0 cumsum(models.thick(I(i),1:end)) 1000],'Color',map(RMS_fact_map(i),:),'LineWidth',(1-RMS_fact(i))*3);
end
% % Adding the true model if existing:
% if isfield(handles.Data,'true');
%     stairs([handles.Data.true.Vs(1) handles.Data.true.Vs],[0 cumsum(handles.Data.true.thick) 1000],':w','LineWidth',2);
% end
set (gca,'Ydir','reverse');
xlabel('V_s [m/s]');
xlim([0 2500]);
ylim([0 500]);
ylabel('Depth [m]');
hold off;
grid on;
subplot(10,3,[3:3:7*3+1]);% Rho
hold on;
for i = 1 : 1 : size(HpostCoeff,1),
    stairs([models.rho(I(i),1) models.rho(I(i),:)],[0 cumsum(models.thick(I(i),1:end)) 1000],'Color',map(RMS_fact_map(i),:),'LineWidth',(1-RMS_fact(i))*3);
end
% % Adding the true model if existing:
% if isfield(handles.Data,'true'),
%     stairs([handles.Data.true.rho(1) handles.Data.true.rho],[0 cumsum(handles.Data.true.thick) 1000],':w','LineWidth',2);
% end
set (gca,'Ydir','reverse');
xlabel('\rho [kg/m^3]');
xlim([0 3000]);
ylim([0 500]);
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

%% Show some dispersion curves:

c_true = MASW.Data.Curve.ct;
f_true = MASW.Data.Curve.f;
c_prior = MASW.Models.model.results;
c_post = MASW.Solution.model.results;

figure;
hold on;
% Adding prior data:
pas = 1;
for i = 1 : pas : length(c_prior),
    if max(c_prior(i,:))<1500 && min(c_prior(i,:))>0,
        plot(f_true,c_prior(i,:),'color',[0.5 0.5 0.5],'linewidth',0.2);
    end
end

% Adding posterior data for the 10 most probable models (lowest RMS):
for i = 1 : length(c_post),
%     if max(c_post(I(i),:))<1500 && min(c_post(I(i),:))>0,
        plot(f_true,c_post(I(i),:),'Color',map(RMS_fact_map(i),:),'linewidth',(1-RMS_fact(i))*3);
%     end
end

plot(f_true,c_true,'w','linewidth',2);
xlabel('Frequency [Hz]');
ylabel('Phase Velocity [m/s]');
