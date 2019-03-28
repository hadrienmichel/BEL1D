function mdata = ForwardModelling_FID(mdata,kdata)
% FORAWRDMODELLING_FID is a fucntion that computes the SNMR response for a
% given soil constritution and experiment configuration. 
%
% It is a modified version of the MRSMatlab  MRSModelling GUI (Müller-Petke
% et al., 2016).
%
% The function takes as arguments:
%   - mdata: a structure containing the model data (Z discretization, water
%            contents, etc.) in the MRSmatlab formatting
%   - kdata: a structure containing the kernel fucntion in the MRSmatlab
%            formatting
%
% It outputs an updated fversion of the mdata structure containing the SNMR
% response in mdata.dat.fid1.

mdata.dat.fid1=[];
mdata.dat.fid1_rot=[];
mdata.dat.fid2=[];
mdata.dat.fid2_rot=[];

fi   = zeros(length(kdata.model.z),1);
T2si = zeros(length(kdata.model.z),1);
T1i  = zeros(length(kdata.model.z),1);
for ilayer = 1:mdata.mod.Nlayer
    inthislayer       = find(mdata.mod.zlayer(ilayer) <= kdata.model.z & kdata.model.z < mdata.mod.zlayer(ilayer+1));
    fi(inthislayer)   = mdata.mod.f(ilayer);
    T2si(inthislayer) = mdata.mod.T2s(ilayer);
    T1i(inthislayer)  = mdata.mod.T1(ilayer);
end
fi(end)   = mdata.mod.f(end);
T2si(end) = mdata.mod.T2s(end);
T1i(end)  = mdata.mod.T1(end);

switch kdata.measure.pulsesequence
    case 'FID'
        ft = zeros(length(kdata.model.z), length(mdata.mod.tfid1));
        for m = 1:length(kdata.model.z)
            ft(m,:) = fi(m) * exp(-mdata.mod.tfid1./T2si(m));
        end
        mdata.dat.v0     = kdata.K * fi;
        mdata.dat.fid1   = kdata.K * ft;
        mdata.dat.fid1   = mdata.dat.fid1 + complex(mdata.mod.noise*randn(size(mdata.dat.fid1)),...
                           mdata.mod.noise*randn(size(mdata.dat.fid1)));
    case 'T1'
        % fid1
        % restrict length to tau
        mdata.mod.tfid1(mdata.mod.tfid1 > mdata.mod.ctau)=[];
        ft = zeros(length(kdata.model.z), length(mdata.mod.tfid1));
        for m = 1:length(kdata.model.z)
            ft(m,:) = fi(m) * exp(-mdata.mod.tfid1./T2si(m));
        end
        mdata.dat.v0     = kdata.K * fi;
        mdata.dat.fid1   = kdata.K * ft;
        mdata.dat.fid1   = mdata.dat.fid1 + complex(mdata.mod.noise*randn(size(mdata.dat.fid1)),...
                           mdata.mod.noise*randn(size(mdata.dat.fid1)));

        % fid2
        ft = zeros(length(kdata.model.z), length(mdata.mod.tfid2));
        for m = 1:length(kdata.model.z)
            ft(m,:) = fi(m) * exp(-mdata.mod.tfid2./T2si(m));
        end
        mdata.dat.v0fid2   = kdata.KT1 * fi;
        mdata.dat.fid2   = kdata.KT1 * ft;
        mdata.dat.fid2   = mdata.dat.fid2 + complex(mdata.mod.noise*randn(size(mdata.dat.fid2)),...
                           mdata.mod.noise*randn(size(mdata.dat.fid2)));
    case 'T2'
        ft        = zeros(length(kdata.model.z), length(mdata.mod.tfid1));
        echotimes = mdata.mod.tau;
        nE        = numel(echotimes);
        for m = 1:length(kdata.model.z)
            for iE=1:nE
                t_echo  = mdata.mod.tfid1 - echotimes(iE) - mdata.mod.tfid1(1);
                ft(m,:) = max(ft(m,:),fi(m) * exp(-echotimes(iE)/T1i(m)).*exp(-(t_echo.^2)/T2si(m)^2));
            end
        end
        mdata.dat.v0     = kdata.K * fi;
        mdata.dat.fid1   = kdata.K * ft;
        mdata.dat.fid1   = mdata.dat.fid1 + complex(mdata.mod.noise*randn(size(mdata.dat.fid1)),...
                           mdata.mod.noise*randn(size(mdata.dat.fid1)));
end

end