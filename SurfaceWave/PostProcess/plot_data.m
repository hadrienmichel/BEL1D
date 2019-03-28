[file, path] = uigetfile('*.maswbel');

load([path file],'-mat');

% lambda  = MASW.Data.Curve.lambda;
c_true  = MASW.Data.Curve.ct;
f_true  = MASW.Data.Curve.f;
% lambda = lambda';

c_prior = MASW.Models.model.results;
% if size(lambda,1) ~= 1,
%     lambda = lambda';
% end
% f_prior = c_prior./repmat(lambda,size(c_prior,1),1);

c_post  = MASW.Solution.model.results;
% f_post = c_post./repmat(lambda,size(c_post,1),1);

% figure;
% hold on;
% %% Adding prior data:
% pas = 1;
% for i = 1 : pas : length(c_prior),
%     if max(c_prior(i,:))<1500 && min(c_prior(i,:))>0,
%         plot(lambda,c_prior(i,:),'color',[0.5 0.5 0.5],'linewidth',0.2);
%     end
% end
% 
% %% Adding posterior data:
% for i = 1 : pas : length(c_post),
%     if max(c_post(i,:))<1500 && min(c_post(i,:))>0,
%         plot(lambda,c_post(i,:),'b','linewidth',0.2);
%     end
% end
% 
% %% Adding true data:
% plot(lambda,c_true,'r','linewidth',1);
% 
% xlabel('Wavelength [m]');
% ylabel('Phase Velocity [m/s]');

% figure;
% hold on;
% %% Adding prior data:
% for i = 1 : pas : length(c_prior),
%     if max(c_prior(i,:))<1500 && min(c_prior(i,:))>0,
%         plot(f_prior(i,:),c_prior(i,:),'color',[0.5 0.5 0.5],'linewidth',0.2);
%     end
% end
% 
% %% Adding posterior data:
% for i = 1 : pas : length(c_post),
%     if max(c_post(i,:))<1500 && min(c_post(i,:))>0,
%         plot(f_post(i,:),c_post(i,:),'b','linewidth',0.2);
%     end
% end
% 
% %% Adding true data:
% plot(f_true,c_true,'r','linewidth',1);
% 
% xlabel('Frequency [Hz]');
% ylabel('Phase Velocity [m/s]');

figure;
hold on;
%% Adding prior data:
pas = 1;
for i = 1 : pas : length(c_prior),
    if max(c_prior(i,:))<1500 && min(c_prior(i,:))>0,
        plot(f_true,c_prior(i,:),'color',[0.5 0.5 0.5],'linewidth',0.2);
    end
end

%% Adding posterior data:
for i = 1 : pas : length(c_post),
    if max(c_post(i,:))<1500 && min(c_post(i,:))>0,
        plot(f_true,c_post(i,:),'b','linewidth',0.2);
    end
end

%% Adding true data:
plot(f_true,c_true,'r','linewidth',1);

xlabel('Frequency [Hz]');
ylabel('Phase Velocity [m/s]');

% close all
% clear all
% clc
% finished = 0;
% tmp = 1;
% 
% while finished ~= 1,
%     [file, path] = uigetfile('*.maswpfa');
% 
%     load([path file],'-mat');
% 
%     X  = MASW.Data.Curve.lambda;
%     Y  = MASW.Data.Curve.ct;
% 
%     Y_prior = MASW.Models.model.results;
%     Y_post = MASW.Solution.model.results;
% %     figure;
% %     hold on;
% %     %% Adding prior data:
% %     pas = 5;
% %     for i = 1 : pas : size(Y_prior,1),
% %         plot(X,Y_prior(i,:),'color',[0.5 0.5 0.5],'linewidth',0.2);
% %     end
% % 
% %     %% Adding posterior data:
% %     for i = 1 : pas : length(Y_post),
% %         plot(X,Y_post(i,:),'b','linewidth',0.2);
% %     end
% % 
% %     %% Adding true data:
% %     plot(X,Y,'r','linewidth',1);
% % 
% %     xlabel('Wavenumber [/]');
% %     ylabel('Phase velocity [m/s]');
%     
%     n_model(tmp) = size(Y_prior,1);
%     RMS = rms(Y_post-repmat(Y',size(Y_post,1),1),2);
%     min_RMS(tmp) = min(RMS);
%     P5_RMS(tmp) = quantile(RMS,0.05);
%     P95_RMS(tmp) = quantile(RMS,0.95);
%     max_RMS(tmp) = max(RMS);
%     param = [MASW.Solution.model.thick MASW.Solution.model.Vp MASW.Solution.model.Vs MASW.Solution.model.rho];
%     std_param(tmp,:) = std(param);
%     mean_param(tmp,:) = mean(param);
%     
%     clear('X','Y','Y_post','Y_prior','MASW');
%     tmp = tmp + 1;
%     
%     answer = questdlg('Would you like to add a solution?','Yes','No');
%     if strcmp(answer,'No')
%         finished = 1;
%     end
% end
% 
% %% Graphs of evolution:
% 
% intervals = [2.5 7.5; 300 1200; 1200 2500; 100 500; 400 1200; 1500 2000; 1900 2400];
% std_default = sqrt(((intervals(:,2)-intervals(:,1)).^2)/12);
% 
% figure
% plot(n_model,min_RMS,'linewidth',3);
% hold on;
% plot(n_model,P5_RMS,'linewidth',3);
% plot(n_model,P95_RMS,'linewidth',3);
% plot(n_model,max_RMS,'linewidth',3);
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% xlabel('Number of models in the prior population');
% ylabel('Minimum RMS in the posterior population [m/s]');
% set(gca,'FontSize',14)
% grid on
% 
% figure; 
% plot(n_model,std_param./repmat(std_default',size(std_param,1),1),'linewidth',3);
% set(gca,'xscale','log')
% xlabel('Number of models in the prior population');
% ylabel('Normalized standard deviation');
% legend('e_1','V_{P,1}','V_{P,2}','V_{S,1}','V_{S,2}','\rho_1','\rho2')
% grid on
% set(gca,'FontSize',14)
% 
% figure; 
% plot(n_model,mean_param./repmat([5 600 2000 300 1000 1600 2000],size(mean_param,1),1),'linewidth',3);
% set(gca,'xscale','log')
% xlabel('Number of models in the prior population');
% ylabel('Normalized mean value');
% legend('e_1','V_{P,1}','V_{P,2}','V_{S,1}','V_{S,2}','\rho_1','\rho2')
% grid on
% set(gca,'FontSize',14)


	

