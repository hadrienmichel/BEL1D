% PLOTQTVSBEL is a function that displays the results from the QT inversion
% on top of the produced models from the BEL modeling framework. 

close all;
clear all;
clc;

%% 1) Get BELmod results:

    [file,path] = uigetfile('*.mrsbel','Select the BEL modeling result file');
    SNMR = load([path file],'-mat');
    SNMR = SNMR.SNMR;

%% 2) Get QTInv results:

    [file,path] = uigetfile('*.mrsi','Select the QT inversion result file');
    idata = load([path file],'-mat');
    idata = idata.idata;
    irun = 1;
    T2 = idata.inv1Dqt.smoothMono.solution(irun).T2.*1000; % T2 in ms
    W = idata.inv1Dqt.smoothMono.solution(irun).w ; % W in /

%% 3) Display the results:

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

    pas = 1;
    
    figure(1);
    subplot(10,2,[1:2:13]);% Water content
    hold on;
    for i = 1 : pas : size(HpostCoeff,1),
        stairs([models.water(I(i),1) models.water(I(i),:)],[models.depth(I(i),:) 100],'Color',map(RMS_fact_map(i),:),'LineWidth',(1-RMS_fact(i))*3);
    end
    stairs([W(1); W], [0 idata.inv1Dqt.smoothMono.z],'Color',[0.5 0.5 0.5],'LineStyle',':','LineWidth',2);
    set (gca,'Ydir','reverse');
    xlabel('Water content (\phi) [/]');
    xlim([0 1]);
    ylim([0 10]);
    ylabel('Depth [m]');
    grid on;
    hold off;
    set(gca,'FontSize',12);

    subplot(10,2,[2:2:14]);% Relaxation time
    hold on;
    for i = 1 : pas : size(HpostCoeff,1),
        stairs([models.T2(I(i),1) models.T2(I(i),:)],[models.depth(I(i),:) 100],'Color',map(RMS_fact_map(i),:),'LineWidth',(1-RMS_fact(i))*3);
    end
    stairs([T2(1); T2], [0 idata.inv1Dqt.smoothMono.z],'Color',[0.5 0.5 0.5],'LineStyle',':','LineWidth',2);
    set (gca,'Ydir','reverse');
    xlabel('Relaxation time (T_2^*) [ms]');
    xlim([0 400]);
    ylim([0 10]);
    ylabel('Depth [m]');
    grid on;
    hold off;
    set(gca,'FontSize',12);

    subplot(10,2,19:20);
    tmp = min(RMS) : (max(RMS)-min(RMS))/999 : max(RMS);
    for i = 1 : length(tmp),
        tmp(i) = find(RMS_Interval<=tmp(i),1,'last');
    end
    imagesc(RMS_Interval,0:1,[(tmp)' (tmp)']');
    colormap(map);
    set(gca,'YTick',[]);
    answer = inputdlg('Enter noise level [nV]:','Noise input',1,{'0'});
    noise = str2double(answer{1});
    yl = ylim;
    if noise ~= 0,
        hold on;
        plot([noise noise],yl,'k','LineWidth',2);
    end
    xlabel('RMS error [nV]');
    set(gca,'FontSize',10);
    
%% 4) Maximum likelihood models and QT inversion (only for 2-layers models)

if SNMR.Models.nbLayers == 2,
    answer = inputdlg('Enter the approximate time of computation [min]:','Time input',1,{'1'});
    time_max = str2double(answer{1})*60;
    time_10_5 = 5;
    time_10_0 = time_10_5/(10^5);
    nb_compute = time_max/time_10_0;
    n = floor(nb_compute^(1/5));
    fprintf('The number of discretization is %d \n',n);
    nb_param = 5;
    proba = zeros(n,n,n,n,n);
    % For RMS plot
    RMS = rms(data_post.*1e9-repmat(data_true.*1e9,size(data_post,1),1),2);

    min_val = (min(HpostCoeff))-0.1.*(min(HpostCoeff));
    min_val(min_val<0) = 0;
    max_val = (max(HpostCoeff))+0.1.*(max(HpostCoeff));
    
    Models_proba = zeros(n^nb_param,nb_param);
    Models_proba_exist = false(n^nb_param,1);
    Models_proba_value = zeros(n^nb_param,1);
    
    for i = 1 : nb_param,
        tmp = min_val(i) : (max_val(i)-min_val(i))/n : max_val(i);
        bounds_min(:,i) = tmp(1:end-1)';
        bounds_max(:,i) = tmp(2:end)';
    end
    mids = (bounds_min+bounds_max)./2;
    tic;
    tmp = 0;
    w = waitbar(0,{'Please wait . . .','Computing the probability'});
    for h = 1 : n,
        for i = 1 : n,
            for j = 1 : n,
                for k = 1 : n,
                    for l = 1 : n,
                        tmp = tmp + 1;
                        if mod(tmp,100)==0,
                            waitbar(tmp/n^nb_param);
                        end
                        vect_min = [bounds_min(h,1) bounds_min(i,2) bounds_min(j,3) bounds_min(k,4) bounds_min(l,5)];
                        vect_max = [bounds_max(h,1) bounds_max(i,2) bounds_max(j,3) bounds_max(k,4) bounds_max(l,5)];
                        index1 = HpostCoeff>repmat(vect_min,size(HpostCoeff,1),1);
                        index2 = HpostCoeff<repmat(vect_max,size(HpostCoeff,1),1);
                        index = index1&index2;
                        index = sum(index,2);
                        index = index==nb_param;
                        proba(h,i,j,k,l) = sum(index)/size(HpostCoeff,1);
                        Models_proba(tmp,:) = [mids(h,1) mids(i,2) mids(j,3) mids(k,4) mids(l,5)];
                        if proba(h,i,j,k,l) > 0,
                            Models_proba_exist(tmp) = true;
                            Models_proba_value(tmp) = sum(index);
                        else
                            Models_proba_value(tmp) = nan;
                        end
                    end
                end
            end
        end
    end
    proba(proba==0)=nan;
    elapsedTime = toc;
    delete(w);
    
%% 5) Display the models with probability

    [Models_proba_value_NONAN,I] = sort(Models_proba_value(~isnan(Models_proba_value)),'ascend');
    Models_proba = Models_proba(~isnan(Models_proba_value),:);
    map = colormap(jet(256));
    nb_map = length(map);
    Q_used = (0 : nb_map-1)./(nb_map-1);
    RMS_Interval = quantile(Models_proba_value_NONAN,Q_used);
    RMS_fact_map = ones(1,length(Models_proba_value_NONAN));
    Q_act = 1;
    for i = length(Models_proba): -1 : 2,
        RMS_fact_map(i) = Q_act;
        if Models_proba_value_NONAN(i-1) >= RMS_Interval(Q_act),
            Q_act = Q_act + 1;
        end
    end
    RMS_fact_map(1) = nb_map;
    RMS_fact = (Models_proba_value_NONAN./(max(Models_proba_value_NONAN)+0.005));

    figure(2);
    subplot(10,2,[1:2:13]);% Water content
    hold on;
    for i = 1 : 1 : size(Models_proba_value_NONAN,1),
        stairs([Models_proba(I(i),2) Models_proba(I(i),[2 3])],[0 Models_proba(I(i),1) 100],'Color',map(RMS_fact_map(i),:),'LineWidth',(1-RMS_fact(i))*5/3);
    end
    stairs([W(1); W], [0 idata.inv1Dqt.smoothMono.z],'Color',[0.5 0.5 0.5],'LineStyle',':','LineWidth',2);
    set (gca,'Ydir','reverse');
    xlabel('Water content (\phi) [/]');
    xlim([0 1]);
    ylim([0 10]);
    ylabel('Depth [m]');
    grid on;
    hold off;
    set(gca,'FontSize',12);

    subplot(10,2,[2:2:14]);% Relaxation time
    hold on;
    for i = 1 : 1 : size(Models_proba_value_NONAN,1),
        stairs([Models_proba(I(i),4) Models_proba(I(i),[4 5])],[0 Models_proba(I(i),1) 100],'Color',map(RMS_fact_map(i),:),'LineWidth',(1-RMS_fact(i))*5/3);
    end
    stairs([T2(1); T2], [0 idata.inv1Dqt.smoothMono.z],'Color',[0.5 0.5 0.5],'LineStyle',':','LineWidth',2);
    set (gca,'Ydir','reverse');
    xlabel('Relaxation time (T_2^*) [ms]');
    xlim([0 400]);
    ylim([0 10]);
    ylabel('Depth [m]');
    grid on;
    hold off;
    set(gca,'FontSize',12);
    
    %%
    subplot(10,2,19:20);
    tmp = max(Models_proba_value_NONAN) : -(max(Models_proba_value_NONAN)-min(Models_proba_value_NONAN))/4999 : min(Models_proba_value_NONAN);
    for i = 1 : length(tmp),
        tmp(i) = find(RMS_Interval<=tmp(i),1,'last');
    end
    imagesc(RMS_Interval,0:1,[(tmp)' (tmp)']');
    colormap(map);
    set(gca,'YTick',[]);
    xlabel('Probability estimation [/]');
    set(gca,'FontSize',10);

end