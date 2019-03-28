function [Y] = forwardVSP(X,param)
% FORWARDVSP is a function that performs the forward modelling for a
% given experiment (in this case, Vertical Seismic Profile [VSP]). 
% The fucntion takes as arguments:
%       - The X coordonates for which the model must be calculated (X)
%       - The parameters of the model in a vector form [1 x ((n x m)-1)]
%           o n is the number of layers
%           o m is the number of parameters
%           The vector is expresed as follow:
% >> param = 
%      
%        [e_1, e_2, ..., e_n-1, param2_1, ..., param2_n, param3_1, ...,
%        param3_n, ..., param4_1, ..., param4_n]
%
% In this case: 
%
% X = dobs : depth of observation point
% Y = t : time of first arrival corresponding to the observed depth.

% param = [e_1, e_2, V_1, V_2, V_3];

N = length(X);
M = N;
dd = X(2)-X(1);
G=tril(dd*ones(N,M));% Forward opperator

nb_layer = ceil(length(param)/2);
V = ones(size(X,1),size(X,2));
% Constructing layered V profile:
for i = 1 : nb_layer,
    if i == 1,
        V((X-X(1))<=param(1)) = param(nb_layer);
    elseif i == nb_layer,
        V((X-X(1))>param(nb_layer-1)) = param(end);
    else
        V((X-X(1))<=param(i) & (X-X(1))>param(i-1)) = param(nb_layer+(i-1));
    end
end
    

S = V.^(-1);
% Forward model
Y = G*S;

end