function [Y] = forwardModel(X,param)
% FORWARDMODEL is a function that performs the forward modelling for a
% given experiment. 
% The fucntion takes as arguments:
%       - The X coordonates for which the model must be calculated (X)
%       - The parameters of the model in a vector form [1 x ((n x m)-1)]
%           o n is the number of layers
%           o m is the number of parameters
%           The vector is expresed as follow:
% >> param = 
%      
%        [e_1, e_2, ..., e_n, param2_1, ..., param2_n, param3_1, ...,
%        param3_n, ..., param4_1, ..., param4_n]
%

% In this example, the forward model is just a non-linear equation

Y = ForwardPendulum(length(X),param(1),param(2),param(3),1.5);
% True model: param = [6 8 20];

%Y = param(1) + X.*param(2) + X.^param(3);

end