close all
clear all
clc
multiple = true;
answer = questdlg('Would you like to observe the impact of the number of models?','Yes','No');
if strcmp(answer,'No')
    multiple = false;
end

if multiple
    finished = 0;
    tmp = 1;

    while finished ~= 1
        [file, path] = uigetfile('*.mrsbel');

        load([path file],'-mat');

        % true
        timing = SNMR.Data.proclog1.Q(1).rx(1).sig(2).t;
        data_true = SNMR.Data.d_real_obs;
        q = [SNMR.Data.proclog1.Q(:).q];
        nb_q = length(q);
        q_used = input(['Which pulse moment to use (nb = ' num2str(nb_q) ') : ']);
        v_used = length(timing)*(q_used-1)+1:length(timing)*q_used;
        data_t = data_true(v_used);

        % prior
        data_prior = squeeze(SNMR.Models.model.results(:,q_used,:));

        %post
        data_post = squeeze(SNMR.Solution.model.results(:,q_used,:));

        figure;
        hold on;
        %% Adding prior data:
        for i = 1 : 5 : size(data_prior,1)
            plot(timing,data_prior(i,:).*1e9,'color',[0.5 0.5 0.5],'linewidth',0.2);
        end

        %% Adding posterior data:
        for i = 1 : 5 : size(data_post,1)
            plot(timing,data_post(i,:).*1e9,'b','linewidth',0.2);
        end

        %% Adding true data:
        plot(timing,data_t.*1e9,'r','linewidth',1);

        xlabel('Time [s]');
        ylabel('Amplitude [nV]');

        data_post = SNMR.Solution.model.results;
        tmp1 = [];

        for i = 1 : size(data_post,2)
            tmp1 = [tmp1 squeeze(data_post(:,i,:))];
        end
        data_post = tmp1;
        clear tmp1;

        n_model(tmp) = size(data_prior,1);
        RMS = rms(data_post.*1e9-repmat(data_true.*1e9,size(data_post,1),1),2);
        min_RMS(tmp) = min(RMS);
        P5_RMS(tmp) = quantile(RMS,0.05);
        P95_RMS(tmp) = quantile(RMS,0.95);
        max_RMS(tmp) = max(RMS);
        param = [SNMR.Solution.model.thick SNMR.Solution.model.water SNMR.Solution.model.T2];
        std_param(tmp,:) = std(param);
        mean_param(tmp,:) = mean(param);

        tmp = tmp + 1;
        clear data_post data_prior

        answer = questdlg('Would you like to add a solution?','Yes','No');
        if strcmp(answer,'No')
            finished = 1;
        end

    end
    %%
    intervals = [2.5 7.5; 3.5/100 10/100; 10/100 30/100; 5 350; 5 350];
    std_default = sqrt(((intervals(:,2)-intervals(:,1)).^2)/12);

    figure
    plot(n_model,min_RMS,'linewidth',3);
    hold on;
    plot(n_model,P5_RMS,'linewidth',3);
    plot(n_model,P95_RMS,'linewidth',3);
    plot(n_model,max_RMS,'linewidth',3);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('Number of models in the prior population');
    ylabel({'RMS error in the ','posterior population [nV]'});
    legend('Minimum','Percentile 5','Percentile 95','Maximum');
    set(gca,'FontSize',14)
    grid on

    figure; 
    plot(n_model,std_param./repmat(std_default',size(std_param,1),1),'linewidth',3);
    set(gca,'xscale','log')
    xlabel('Number of models in the prior population');
    ylabel('Normalized standard deviation');
    legend('e_1','W_{1}','W_{2}','T_{2,1}^*','T_{2,2}^*')
    grid on
    set(gca,'FontSize',14)

    figure; 
    plot(n_model,mean_param./repmat([5 5/100 25/100 100 200],size(mean_param,1),1),'linewidth',3);
    set(gca,'xscale','log')
    xlabel('Number of models in the prior population');
    ylabel('Normalized mean value');
    legend('e_1','W_{1}','W_{2}','T_{2,1}^*','T_{2,2}^*')
    grid on
    set(gca,'FontSize',14)
else
    [file, path] = uigetfile('*.mrsbel');

    load([path file],'-mat');

    % true
    timing = SNMR.Data.proclog1.Q(1).rx(1).sig(2).t;
    data_true = SNMR.Data.d_real_obs;
    q = [SNMR.Data.proclog1.Q(:).q];
    nb_q = length(q);
    %q_used = input(['Which pulse moment to use (nb = ' num2str(nb_q) ') : ']);
    figure('units','normalized','outerposition',[0 0 1 1]);
    nb_line = 4;
    nb_disp = 4;
    nb_col = ceil(nb_disp/nb_line);
    tmp = 1;
    for j = 1 : floor(nb_q/nb_disp) : nb_q,
        subplot(nb_line,nb_col,tmp);
        q_used = j;
        v_used = length(timing)*(q_used-1)+1:length(timing)*q_used;
        data_t = data_true(v_used);

        % prior
        data_prior = squeeze(SNMR.Models.model.results(:,q_used,:));

        %post
        data_post = squeeze(SNMR.Solution.model.results(:,q_used,:));

        hold on;
        %% Adding prior data:
        for i = 1 : 5 : size(data_prior,1)
            plot(timing,data_prior(i,:).*1e9,'color',[0.5 0.5 0.5],'linewidth',0.2);
        end

        %% Adding posterior data:
        for i = 1 : 5 : size(data_post,1)
            plot(timing,data_post(i,:).*1e9,'b','linewidth',0.2);
        end

        %% Adding true data:
        plot(timing,data_t.*1e9,'r','linewidth',1);
        
        title(sprintf('Pulse moment nb. %d',j));
        xlabel('Time [s]');
        ylabel('Amplitude [nV]');
        set(gca,'FontSize',10);
        
        tmp = tmp +1;
    end
end
