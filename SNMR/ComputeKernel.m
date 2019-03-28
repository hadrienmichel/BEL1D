function [kdataall] = ComputeKernel(proclog, channelrx, resEarth, model)
% COMPUTEKERNEL is a function that computes the kernels for a given
% constitution of the earth's subsurface.
%
% It is a modified version of the MRSMatlab  MRSKernel GUI (Müller-Petke et
% al., 2016).
%
% ------------------------------- ATTENTION!-------------------------------
% In order to run, the function needs acces to the MRSMatlab
% functions!!!
% -------------------------------ATTENTION!-------------------------------
%
% The function takes as arguments:
%   - proclog: a loaded MRS dataset in the MRSMatlab formatting
%   - channelrx: the channel of the computed receiver
%   - resEarth: a flag for resistive earth (1 = inf resistivity, 0 = custom
%               resistivity, given in the model structure)
%   - model (optional): a structure containing the informations about the
%                       subsurface that are relevent to the computation of
%                       the kernels. If used, this structure MUST conform
%                       to the struct elements:
%
%                                            dh: 1
%                                          dphi: 0.0175
%                                          zmax: 75
%                                       z_space: 1
%                                            dz: 1
%                                       zlogmin: 0.1000
%                                       zlogmax: 1
%                                         hmaxx: 200
%                                         hmaxy: 200
%                                            nz: 96
%                                     sinh_zmin: 0.0200
%                                      LL_dzmin: 0.0200
%                                      LL_dzmax: 0.2000
%                                       LL_dlog: 2
%                                             z: [1×N double]
%                                            Dz: [1×N double]
%
%                       where N is the number of discretisations in the z
%                       direction.
%
% The fucntion outputs the kdataall structure, containing the kernel in the
% MRSMatlab formatting.

if length(channelrx) == 1 && nargin < 4,
    %% Kernel computation
    kdata.gammaH             = +0.267518*1e9;
    kdata.K                  = [];
    kdata.mes_conf           = 1;
    
    % Loop characteristics
    kdata.loop.shape         = proclog.rxinfo(channelrx).looptype; % circ:1 eight:2
    kdata.loop.eightoritn    = 0; % inline:1 normal:2
    kdata.loop.eightsep      = 0; % for 8 only
    % For inloop parameter, insert [Dtx, Drx]
    % Accessing to tx loopsiz : proclog.txinfo(1).loopsize and proclog.txinfo(1).loopturns
    kdata.loop.size          = [proclog.rxinfo(channelrx).loopsize]; % diameter for circular, side for square
    kdata.loop.turns         = [proclog.txinfo(1).loopturns proclog.rxinfo(channelrx).loopturns]; % number of turns [transmitter receiver]
    kdata.loop.rmax          = 300;
    kdata.loop.dr            = 24;
    
    % Earth model
    kdata.earth.inkl         = 60;% inclinaison of ambiant magnetici field
    kdata.earth.decl         = 0;
    erdt = 0;%input('Please, refer the experimental ambiant magnetic field [nT] (if default values, type 0): ');
    if erdt ~= 0,
        kdata.earth.erdt     = erdt*1e-9;
    else
        kdata.earth.erdt         = 48000*1e-9;% B_0 : ambiant magnetic field (in nT)
    end
    kdata.earth.f            = -kdata.gammaH*kdata.earth.erdt/(2*pi);% Larmor frequency computed (rad/s)
    kdata.earth.w_rf         = -kdata.gammaH*kdata.earth.erdt;% Larmor frequency computed Hz
    kdata.earth.res          = resEarth;%1 = inf resistivity, 0 = custom resistivity
    kdata.earth.nl           = 1;
    % If wanna insert resistivity model --> HERE
    % res = input('Do you want to add a resistivity model (0=no, 1=yes)?');
    % if res ~= 0,
    %   kdata.earth.zm = input('Refer depth to interfaces (if no, press ENTER): ');
    %   kdata.earth.sm = input('Refer resistivity of layers (at least one): ');
    %   kdata.earth.sm = kdata.earth.sm.^(-1);
    % else
    kdata.earth.zm           = []; % Depth
    kdata.earth.sm		    = 0.01; % Conductivity (0.01 => 100 Ohm.m)
    % end
    kdata.earth.fm		    = 1.00; % Water content
    kdata.earth.temp         = 281; % aquifer temperature
    kdata.earth.type         = 1; % 1: single pulse; 2: double pulse kernel; 3: single-pulse df kernel
    
    % Model characteristics
    kdata.model.dh           =        1;   %
    kdata.model.dphi         = 1*pi/180;   % for pol cordinates   [deg]
    kdata.model.zmax         = 1.5*proclog.txinfo(1).loopsize;   % for both             [m]
    kdata.model.z_space      =        1;   % sinh:1 loglin:2
    kdata.model.dz           =        1;   % layer thickness, thickness of first layer for log
    kdata.model.zlogmin      =      0.1;
    kdata.model.zlogmax      =        1;
    kdata.model.hmaxx        =      200;
    kdata.model.hmaxy        =      200;
    
    % Measure characteristics
    kdata.measure.dim        =        1;   % Dimension of Kernel
    kdata.measure.pm_vec     =  [proclog.Q(:).q];% Magnitude of induces magnetic pulses!
    kdata.measure.pm_vec_2ndpulse =  kdata.measure.pm_vec; % second pulse pulse moment
    kdata.measure.taud         = 1e9; % inter pulse delay for double pulse
    kdata.measure.taup1        = proclog.Q(1).timing.tau_p1; % duration pulse 1
    % Attention, modified for synthetical data (taup1 instead of taup2)
    kdata.measure.taup2        = proclog.Q(1).timing.tau_p1; % duration pulse 2
    kdata.measure.df           = 0; % frequency offset [Hz]
    kdata.measure.pulsesequence= 'FID';
    
    if kdata.measure.df == 0
        flag_df = 0;
    else
        flag_df = 1;
    end
    
    kdata.model.nz        =  4*length(kdata.measure.pm_vec);% Number of layers in the kernel model = 4*nb_pulses
    kdata.model.sinh_zmin =  kdata.loop.size(1)/500;% Thickness of first layer if sinh
    kdata.model.LL_dzmin  =  kdata.loop.size(1)/500;% Thickness of first laer if loglin
    kdata.model.LL_dzmax  =  kdata.loop.size(1)/50;% Thickness of last layer if loglin
    kdata.model.LL_dlog   =  kdata.loop.size(1)/5;% Depth of log
    kdata.model           =  MakeZvec(kdata.model);% Constructing the depth model
    if proclog.txinfo.channel ~= proclog.rxinfo(channelrx).channel,
        %If wanna compute for an inloop :
        kdata.loop.shape  = 5;
        kdata.loop.size   = [proclog.txinfo(1).loopsize proclog.rxinfo(channelrx).loopsize]; % diameter for circular, side for square
    end
    %             % Checking that the model is computable :
    %             [Computable,depth_max] = Check_computability_limit(kdata.loop.size(1),kdata.loop.size(2),1.5);
    %             if Computable == false,
    %                 index_limit = find(kdata.model.z < depth_max);
    %                 save.model = kdata.model;
    %                 kdata.model.nz = length(index_limit);
    %                 kdata.model.z = save.model.z(index_limit);
    %                 kdata.model.Dz = save.model.Dz(index_limit);
    %             end
    %         else
    %             Computable = true;
    %         end
    % [kdata.K, J, kdata.B1] = MakeKernel(kdata.loop, kdata.model, kdata.measure,kdata.earth);
    
    % Trying to avoid the special cases :
    Computable = true;
    IS_computed = false;
    save.model = kdata.model;
    % nz_previous = save.model.nz;
    nz_current = save.model.nz;
    while not(IS_computed),
        kdata.model.nz = nz_current;
        kdata.model.z = save.model.z(1:kdata.model.nz);
        kdata.model.Dz = save.model.Dz(1:kdata.model.nz);
        try
            [kdata.K, J, kdata.B1] = MakeKernel(kdata.loop, kdata.model, kdata.measure,kdata.earth);
        catch
            Computable = false;
            % Not computable --> too deep --> reducing max.depth
            nz_current = nz_current - 1;
            continue;
        end
        IS_computed = true;
    end
    
    
    % End - Trying to avoid the special cases
    if Computable == false,
        K_tmp = zeros(length(kdata.measure.pm_vec),4*length(kdata.measure.pm_vec));
        K_tmp(:,1:kdata.model.nz) = kdata.K;
        kdata.K = K_tmp;
        kdata.model = save.model;
    end
    
    %% Show kernel
    figure;
    weight = repmat(kdata.model.Dz,size(kdata.K,1),1);
    pcolor(kdata.measure.pm_vec, kdata.model.z, (abs(kdata.K)./weight).');
    % Weight gives more importance to the shallow layers
    axis ij
    shading flat
    box on
    grid on
    set(gca, 'layer', 'top')
    title(['Kernel for Tx/Rx' num2str(channelrx)])
    xlabel('Pulse moment q [As]')
    ylabel('Depth [m]')
    hold on;
    
    % Find max
    max_K = ones(1,length(kdata.measure.pm_vec));
    for i = 1 : length(kdata.measure.pm_vec),
        max_K(i) = find(abs(kdata.K(i,:))./weight(i,:) == max(abs(kdata.K(i,3:end))./weight(i,3:end)));
    end
    plot(kdata.measure.pm_vec, kdata.model.z(max_K),'r');
    %% Return kernel
    kdataall = kdata;
else
    channelsave = channelrx;
    passed = 0;
    for channelrx = channelsave,
        if passed ~= 1,
            %% Kernel computation
            kdata.gammaH             = +0.267518*1e9;
            kdata.K                  = [];
            kdata.mes_conf           = 1;
        end
        % Loop characteristics
        kdata.loop.shape         = proclog.rxinfo(channelrx).looptype; % circ:1 eight:2
        kdata.loop.eightoritn    = 0; % inline:1 normal:2
        kdata.loop.eightsep      = 0; % for 8 only
        % For inloop parameter, insert [Dtx, Drx]
        % Accessing to tx loopsiz : proclog.txinfo(1).loopsize and proclog.txinfo(1).loopturns
        kdata.loop.size          = [proclog.rxinfo(channelrx).loopsize]; % diameter for circular, side for square
        kdata.loop.turns         = [proclog.txinfo(1).loopturns proclog.rxinfo(channelrx).loopturns]; % number of turns [transmitter receiver]
        kdata.loop.rmax          = 300;
        kdata.loop.dr            = 24;
        if passed ~= 1,
            % Earth model
            kdata.earth.inkl         = 60;% inclinaison of ambiant magnetici field
            kdata.earth.decl         = 0;
            erdt = 0;%input('Please, refer the experimental ambiant magnetic field [nT] (if default values, type 0): ');
            if erdt ~= 0,
                kdata.earth.erdt     = erdt*1e-9;
            else
                kdata.earth.erdt         = 48000*1e-9;% B_0 : ambiant magnetic field (in nT)
            end
            kdata.earth.f            = -kdata.gammaH*kdata.earth.erdt/(2*pi);% Larmor frequency computed (rad/s)
            kdata.earth.w_rf         = -kdata.gammaH*kdata.earth.erdt;% Larmor frequency computed Hz
            kdata.earth.nl           = 1;
            % If wanna insert resistivity model --> HERE
            res = input('Do you want to add a resistivity model (0=no, 1=yes)?');
            kdata.earth.res          = not(res);%1 = inf resistivity, 0 = custom resistivity
            if res ~= 0,
                kdata.earth.zm = input('Refer depth to interfaces (if no, press ENTER): ');
                kdata.earth.sm = input('Refer resistivity of layers (at least one): ');
            else
                kdata.earth.zm = [];
                kdata.earth.sm = 0.01;
            end
            % end
            kdata.earth.fm		    = 1.00; % Water content
            kdata.earth.temp         = 281; % aquifer temperature
            kdata.earth.type         = 1; % 1: single pulse; 2: double pulse kernel; 3: single-pulse df kernel
        end
        
        % Measure characteristics
        kdata.measure.dim        =        1;   % Dimension of Kernel
        kdata.measure.pm_vec     =  [proclog.Q(:).q];% Magnitude of induces magnetic pulses!
        kdata.measure.pm_vec_2ndpulse =  kdata.measure.pm_vec; % second pulse pulse moment
        kdata.measure.taud         = 1e9; % inter pulse delay for double pulse
        kdata.measure.taup1        = proclog.Q(1).timing.tau_p1; % duration pulse 1
        % Attention, modified for synthetical data (taup1 instead of taup2)
        kdata.measure.taup2        = proclog.Q(1).timing.tau_p1; % duration pulse 2
        kdata.measure.df           = 0; % frequency offset [Hz]
        kdata.measure.pulsesequence= 'FID';
        
        if kdata.measure.df == 0
            flag_df = 0;
        else
            flag_df = 1;
        end
        
        kdata.model = model;
        
        if proclog.txinfo.channel ~= proclog.rxinfo(channelrx).channel,
            %If wanna compute for an inloop :
            kdata.loop.shape  = 5;
            kdata.loop.size   = [proclog.txinfo(1).loopsize proclog.rxinfo(channelrx).loopsize]; % diameter for circular, side for square
        end
        % [kdata.K, J, kdata.B1] = MakeKernel(kdata.loop, kdata.model, kdata.measure,kdata.earth);
        
        % Trying to avoid the special cases :
        Computable = true;
        IS_computed = false;
        save.model = kdata.model;
        % nz_previous = save.model.nz;
        nz_current = save.model.nz;
        while not(IS_computed),
            kdata.model.nz = nz_current;
            kdata.model.z = save.model.z(1:kdata.model.nz);
            kdata.model.Dz = save.model.Dz(1:kdata.model.nz);
            try
                [kdata.K, J, kdata.B1] = MakeKernel(kdata.loop, kdata.model, kdata.measure,kdata.earth);
            catch
                Computable = false;
                % Not computable --> too deep --> reducing max.depth
                nz_current = nz_current - 1;
                continue;
            end
            IS_computed = true;
        end
        
        
        % End - Trying to avoid the special cases
        if Computable == false,
            K_tmp = zeros(length(kdata.measure.pm_vec),save.model.nz);
            K_tmp(:,1:kdata.model.nz) = kdata.K;
            kdata.K = K_tmp;
            kdata.model = save.model;
        end
        %% Show kernel
        figure;
        weight = repmat(kdata.model.Dz,size(kdata.K,1),1);
        pcolor(kdata.measure.pm_vec, kdata.model.z, (abs(kdata.K)./weight).');
        % Weight gives more importance to the shallow layers
        axis ij
        shading flat
        box on
        grid on
        set(gca, 'layer', 'top')
        title(['Kernel for Tx/Rx' num2str(channelrx)])
        xlabel('Pulse moment q [As]')
        ylabel('Depth [m]')
        hold on;
        
        % Find max
        max_K = ones(1,length(kdata.measure.pm_vec));
        for i = 1 : length(kdata.measure.pm_vec),
            max_K(i) = find(abs(kdata.K(i,:))./weight(i,:) == max(abs(kdata.K(i,3:end))./weight(i,3:end)));
        end
        plot(kdata.measure.pm_vec, kdata.model.z(max_K),'r');
        drawnow
        %% Return kernel
        eval(['kdataall.kdata' num2str(channelrx) '= kdata;']);
        passed = 1;
    end
end
end
