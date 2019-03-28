function [models] = ModelGenerator(type, N, parameters, nb_layer)
% MODELGENERATOR is a fucntion that generates models for the SNMR
% experiment problem.
%
% The function takes as arguments:
%   - type: 1 = Uniformly distributed
%           2 = Uniformly distributed with a Latin-Hypercube sampler
%           3 = Normal distribution
%           4 = Normal distribution with a Latin-Hypercube sampler
%
% The normaly distributed variables are computed with the parameters:
%   - mean = (min+max)/2
%   - std = (max-min)/4
% 
% The function outputs a structure model containing the sampled models.

models = struct('thick',zeros(N,nb_layer-1),'depth',zeros(N,nb_layer),'water',zeros(N,nb_layer),'T2',zeros(N,nb_layer),'res',ones(N,1));
sigma = 100;%input('Enter the constant resistivity for the model [in Ohm.m] : '); % in Ohm.m

thick = parameters(1:end-1,1:2);
water = parameters(:,3:4);
T2 = parameters(:,5:6);

if type == 1,% Uniformly distributed
    % Computing parameters
    models.thick = ones(N,nb_layer-1)*diag(thick(:,1)) + rand(N,nb_layer-1)*(diag(thick(:,2)-thick(:,1)));
    models.water = ones(N,nb_layer) * diag(water(:,1)) + rand(N,nb_layer)*(diag(water(:,2)-water(:,1)));
    models.T2 = ones(N,nb_layer)*diag(T2(:,1)) + rand(N,nb_layer)*(diag(T2(:,2)-T2(:,1)));
    models.res(:) = ones(N,1)*sigma;
elseif type == 2,% Latin Hypercube Sampler uniformly distributed
    % Computing parameters
    temporary = lhsdesign(N,nb_layer*3); % Latin hypercube sampling
    models.thick = ones(N,nb_layer-1)*diag(thick(:,1)) + temporary(:,1:nb_layer-1)*(diag(thick(:,2)-thick(:,1)));
    models.water = ones(N,nb_layer) * diag(water(:,1)) + temporary(:,nb_layer:(2*nb_layer)-1)*(diag(water(:,2)-water(:,1)));
    models.T2 = ones(N,nb_layer)*diag(T2(:,1)) + temporary(:,2*nb_layer:end-1)*(diag(T2(:,2)-T2(:,1)));
    models.res(:) = ones(N,1)*sigma;
elseif type == 3,% Normal distribution (col1 = mu, col2 = sigma)
    tmp1 = (thick(:,1)+thick(:,2))./2;
    tmp2 = abs(thick(:,2)-thick(:,1))./4;
    thick = [tmp1 tmp2];
    tmp1 = (water(:,1)+water(:,2))./2;
    tmp2 = abs(water(:,2)-water(:,1))./4;
    water = [tmp1 tmp2];
    tmp1 = (T2(:,1)+T2(:,2))./2;
    tmp2 = abs(T2(:,2)-T2(:,1))./4;
    T2 = [tmp1 tmp2];
    % Computing parameters
    models.thick = normrnd(repmat(thick(:,1)',N,1),repmat(thick(:,2)',N,1));
    models.water = normrnd(repmat(water(:,1)',N,1),repmat(water(:,2)',N,1));
    models.T2(:,:) = normrnd(repmat(T2(:,1)',N,1),repmat(T2(:,2)',N,1));
    models.thick(models.thick<0) = 0;
    models.water(models.water<0) = 0;
    models.T2(models.T2<0) = 0;
    models.res(:) = ones(N,1)*sigma;
elseif type == 4,% Latin hypecube sampler Normal distribution (col1 = mu, col2 = sigma)
    tmp1 = (thick(:,1)+thick(:,2))./2;
    tmp2 = abs(thick(:,2)-thick(:,1))./4;
    thick = [tmp1 tmp2];
    tmp1 = (water(:,1)+water(:,2))./2;
    tmp2 = abs(water(:,2)-water(:,1))./4;
    water = [tmp1 tmp2];
    tmp1 = (T2(:,1)+T2(:,2))./2;
    tmp2 = abs(T2(:,2)-T2(:,1))./4;
    T2 = [tmp1 tmp2];
    % Computing parameters
    temporary = lhsnorm([thick(:,1)', water(:,1)', T2(:,1)'],diag([thick(:,2)', water(:,2)', T2(:,2)']),N);
    models.thick = temporary(:,1:nb_layer-1);
    models.water = temporary(:,nb_layer:(nb_layer*2)-1);
    models.T2 = temporary(:,nb_layer*2:end);
    models.thick(models.thick<0) = 0;
    models.water(models.water<0) = 0;
    models.T2(models.T2<0) = 0;
    models.res(:) = ones(N,1)*sigma;
else
    error('You didn''t choose a satistical law');
end

% Constructing depth
for i = 2 : nb_layer,
    models.depth(:,i) = models.depth(:,i-1)+ models.thick(:,i-1);
end

end