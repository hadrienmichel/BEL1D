%% Hadrien MICHEL : february 2018
%
% Exploring the data-set may be complicated due to the large number of
% values. This script enables a fast and easy way to observe directly the
% raw data with graphs. 
%
% __________

% Loading the file containing the data

[file, path] = uigetfile('*.mrsd');
load([path file],'-mat');
%load('mergedTx50_2L_NF_NOISY2.mrsd','-mat');

% Constituting the graphs for the different receptors

pm_vect = [proclog.Q(:).q];
for j = 1 : 1,%length(proclog.rxinfo),
    %fig = figure;
    %subplot(2,1,1);
    hold on;
    for i = 1 : length(pm_vect),
        plot3(proclog.Q(i).rx(j).sig(2).t,ones(1,length(proclog.Q(i).rx(j).sig(2).t))*pm_vect(i),real(proclog.Q(i).rx(j).sig(2).V)*1e9,'r','Linewidth',1.5);
    end
    xlabel('Time [s]');
    ylabel('Pulse Moment [As]');
    zlabel('FID (real part) [nV]');
    %title(['Rx = Channel ' num2str(proclog.rxinfo(j).channel)]);
    grid on;
    view(45,45);
%     subplot(2,1,2);
%     hold on;
%     for i = 1 : length(pm_vect),
%         plot3(proclog.Q(i).rx(j).sig(2).t,ones(1,length(proclog.Q(i).rx(j).sig(2).t))*pm_vect(i),imag(proclog.Q(i).rx(j).sig(2).V)*1e9);
%     end
%     xlabel('Time [s]');
%     ylabel('Pulse Moment [As]');
%     zlabel('FID (imaginary part) [nV]');
%     %title(['Rx = Channel ' num2str(proclog.rxinfo(j).channel)]);
%     grid on;
%     view(45,45);
%     subtitle(['Rx = Channel ' num2str(proclog.rxinfo(j).channel)]);
%     subplot(2,2,3);
%     hold on;
%     for i = 1 : length(kdata.measure.pm_vec),
%         plot3(proclog.Q(i).rx(j).sig(2).t,ones(1,length(proclog.Q(i).rx(j).sig(2).t))*kdata.measure.pm_vec(i),real(proclog.Q(i).rx(j).sig(2).E));
%     end
%     xlabel('Time [s]');
%     ylabel('Pulse Moment [As]');
%     zlabel('E (real part) [?]');
%     %title(['Rx = Channel ' num2str(proclog.rxinfo(j).channel)]);
%     grid on;
%     subplot(2,2,4);
%     hold on;
%     for i = 1 : length(kdata.measure.pm_vec),
%         plot3(proclog.Q(i).rx(j).sig(2).t,ones(1,length(proclog.Q(i).rx(j).sig(2).t))*kdata.measure.pm_vec(i),imag(proclog.Q(i).rx(j).sig(2).E));
%     end
%     xlabel('Time [s]');
%     ylabel('Pulse Moment [As]');
%     zlabel('E (imaginary part) [?]');
%     grid on;
    
end