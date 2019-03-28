close all;
clear all;
clc;
path_full = pwd;
pathBEL = extractBefore(path_full,'\Others');
addpath([pathBEL '\GlobalFunctions']);

%% Initialization:
rng(0);
X = rand(10000,1).*10;
Y = X.*5+normrnd(0,2,10000,1);

X_true = 5;
error = 1;
Y_comp = -10:70/9999:60;
time = zeros(3,1);%[t_true, t_errorB, t_errorN]

figure;
plot(X,Y,'.');
hold on
yl = ylim;
plot([X_true X_true],yl);
xlabel('X');
ylabel('Y');
legend('Dataset','True value of X');
set(gca,'FontSize',14);

%% Kernel density estimation on X_true:
tic;
KD_True = KernelDensity([X Y],X_true,Y_comp,[0.1 0.1]);
time(1) = toc;

%% Kernel density estimation on X_true+error (bandwidth):
tic;
KD_ErrorBand = KernelDensity([X Y],X_true,Y_comp,[error 0.1]);
time(2) = toc;

%% Kernel density estimation on X_true+error (distribution):
X_compute = normrnd(X_true,error,100,1);
tmp = zeros(100,length(Y_comp));
tic;
for i = 1 : 100
    tmp(i,:) = KernelDensity([X Y],X_compute(i),Y_comp,[0.1 0.1]);
end
KD_ErrorNorm = mean(tmp);
time(3) = toc;

%% Figure
figure;
plot(Y_comp,KD_True);
hold on
plot(Y_comp,KD_ErrorBand);
plot(Y_comp,KD_ErrorNorm);
legend_text = {sprintf('Noise-free (t = %f)',time(1)),sprintf('Noisy: bandwidth (t = %f)',time(2)),sprintf('Noisy: distribution (t = %f)',time(3))};
legend(legend_text,'FontSize',12);
xlabel('Y');
ylabel('Probability estimate');
set(gca,'FontSize',14);

rmpath([pathBEL '\GlobalFunctions']);
