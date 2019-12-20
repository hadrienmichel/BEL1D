function [Dist,out] = KernelDensity(dataset, X_true, Y_comp, error)
% KERNELDENSITY is a function that computes the probability density
% function of a value in a one dimensional space from a two dimensional
% dataset (X,Y), constrained to a given value of X_true (and possibly an
% increased impact of the surrounding through the error parameter) at the Y
% coordonates given in Y_comp.
%
% The function computes the sum of the contribiutions of all the points in
% the dataset, with a different load according to the direction (the load
% along Y is fixed but the load along X can vary according to the error).
%
% The arguments are:
%   - dataset: a N-by-2 matrix containing the dataset in the form [X, Y]
%   - X_true: a scalar value presenting the true value of X constraining
%             the obtained distribution
%   - Y_comp: the Y-coordonates of the points on which the distribution is
%             estimated.
%   - error: either
%               o a scalar value that attributes an error (standard 
%                 deviation) to the value of X_true
%               o a 1-by-2 vector containing [std_X, std_Y]
%
% The default values of the bandwidth are 0.1 in both directions.
% The standard deviations in both directions are adapted with the reduced
% dataset described by the bandwidth.
%
% To call this function, you can type:
%
%   KernelDensity(dataset,X_true): the function outputs the kernel density
%                                  estimator for the probability density 
%                                  function of Y, constrained to X_true.
%   KernelDensity(dataset,X_true,error): the function outputs the
%                                        kernel density estimator of the
%                                        probability density function of
%                                        Y, constrained to X_true,
%                                        accounting for the error
%                                        associated to the value X_true
%                                        through an increased bandwidth in
%                                        the X direction.

bandwidth_X = 0.1;% The standard deviation along the X direction
bandwidth_Y = 0.1;% The standard deviation along the Y direction
if nargin == 4
    if max(size(error)) == 1
        bandwidth_X = error;
    elseif max(size(error)) == 2
        bandwidth_X = error(1);
        bandwidth_Y = error(2);
    else
        error('The error input must be either a scalar or a 1-by-2 vector!');
    end
end
out = false;

% Discarding the points in the dataset with negligeable influence
dataset_modif = dataset(dataset(:,1)<X_true+3*bandwidth_X,:);
dataset_modif = dataset_modif(dataset_modif(:,1)>X_true-3*bandwidth_X,:);
disp('Testing bandwidth:');
test = 0;
while size(dataset_modif,1) < 0.01*size(dataset,1) && test < 15,
    test = test + 1;
    bandwidth_X = bandwidth_X*2;
    bandwidth_Y = bandwidth_Y*2;
    disp('Bandwidth multiplied by 2 in X and Y direction!');
    dataset_modif = dataset(dataset(:,1)<X_true+3*bandwidth_X,:);
    dataset_modif = dataset_modif(dataset_modif(:,1)>X_true-3*bandwidth_X,:);
end
if size(dataset_modif,1) < 10,
    disp('The true value is outside the prior! Returning to initial bandwidthds...');
    bandwidth_X = 0.1;% The standard deviation along the X direction
    bandwidth_Y = 0.1;% The standard deviation along the Y direction
    if nargin == 4
        if max(size(error)) == 1
            bandwidth_X = error;
        elseif max(size(error)) == 2
            bandwidth_X = error(1);
            bandwidth_Y = error(2);
        else
            error('The error input must be either a scalar or a 1-by-2 vector!');
        end
    end
end
disp('Bandwith is OK');

values = Y_comp;
Dist = zeros(size(values,1),size(values,2));
if isempty(dataset_modif)
    out = true;
    min_dist = min(dataset(:,2))-3*bandwidth_Y;
    max_dist = max(dataset(:,2))+3*bandwidth_Y;
    index_max = values<max_dist;
    index_min = values>min_dist;
    index_Y = index_min&index_max;
    index_Y = find(index_Y~=0);
    % Computing the full distribution
%     tic;
    z = @(x,y,X,Y) (x-X)^2/bandwidth_X^2 + (y-Y)^2/bandwidth_Y^2 - (2*(x-X)*(y-Y)/bandwidth_X*bandwidth_Y);
    pdf_compute = @(x,y,X,Y) (1/(2*pi()*bandwidth_X*bandwidth_Y)*exp(-z(x,y,X,Y)/2));
    values_X = (min(dataset(:,1))-4*bandwidth_X : 2*bandwidth_X : max(dataset(:,1))+4*bandwidth_X);
    tmp = zeros(length(Y_comp),length(values_X));
    for j = 1 : length(values_X),
        dataset_modif = dataset(dataset(:,1)<values_X(j)+3*bandwidth_X,:);
        dataset_modif = dataset_modif(dataset_modif(:,1)>values_X(j)-3*bandwidth_X,:);
        if ~isempty(dataset_modif)
            for i = 1 : length(index_Y),
                for k = 1 : size(dataset_modif,1),
                    tmp(index_Y(i),j) = tmp(index_Y(i),j) + pdf_compute(values_X(j),values(index_Y(i)),dataset_modif(k,1),dataset_modif(k,2));
                end
            end
        end
    end
    Dist = sum(tmp,2);
%     toc
    % End: computing the full distribution
    % Dist(index_Y) = 1;
    Dist = Dist./trapz(Y_comp,Dist);
    Dist(isnan(Dist)) = 0;% When no information is avalaible on a band: divided by 0 but equal to 0.
    
    fprintf('The true value is outside of the dataset or the number of points \n in the dataset is too low to correctly recover a consistent distribution!\n')
    warning(sprintf('The true value is outside of the dataset or the number of points \n in the dataset is too low to correctly recover a consistent distribution!\n'));
    return
else
    min_compute = min(dataset_modif(:,2))-4*bandwidth_Y;
    max_compute = max(dataset_modif(:,2))+4*bandwidth_Y;
    index_min = values>min_compute;
    index_max = values<max_compute;
    index = index_min&index_max;
    if length(find(index))<10,
        out = true;
        min_dist = min(dataset(:,2))-3*bandwidth_Y;
        max_dist = max(dataset(:,2))+3*bandwidth_Y;
        index_max = values<max_dist;
        index_min = values>min_dist;
        index_Y = index_min&index_max;
        index_Y = find(index_Y~=0);
        % Computing the full distribution
    %     tic;
        z = @(x,y,X,Y) (x-X)^2/bandwidth_X^2 + (y-Y)^2/bandwidth_Y^2 - (2*(x-X)*(y-Y)/bandwidth_X*bandwidth_Y);
        pdf_compute = @(x,y,X,Y) (1/(2*pi()*bandwidth_X*bandwidth_Y)*exp(-z(x,y,X,Y)/2));
        values_X = (min(dataset(:,1))-4*bandwidth_X : 2*bandwidth_X : max(dataset(:,1))+4*bandwidth_X);
        tmp = zeros(length(Y_comp),length(values_X));
        for j = 1 : length(values_X),
            dataset_modif = dataset(dataset(:,1)<values_X(j)+3*bandwidth_X,:);
            dataset_modif = dataset_modif(dataset_modif(:,1)>values_X(j)-3*bandwidth_X,:);
            if ~isempty(dataset_modif)
                for i = 1 : length(index_Y),
                    for k = 1 : size(dataset_modif,1),
                        tmp(index_Y(i),j) = tmp(index_Y(i),j) + pdf_compute(values_X(j),values(index_Y(i)),dataset_modif(k,1),dataset_modif(k,2));
                    end
                end
            end
        end
        Dist = sum(tmp,2);
    %     toc
        % End: computing the full distribution
        % Dist(index_Y) = 1;
        Dist = Dist./trapz(Y_comp,Dist);
        Dist(isnan(Dist)) = 0;% When no information is avalaible on a band: divided by 0 but equal to 0.

        fprintf('The true value is outside of the dataset or the number of points \n in the dataset is too low to correctly recover a consistent distribution!\n')
        warning(sprintf('The true value is outside of the dataset or the number of points \n in the dataset is too low to correctly recover a consistent distribution!\n'));
        return
    end
    dataset = dataset_modif;
    clear('dataset_modif');
end

% Computing the kernel density estimate:
z = @(x,y,X,Y) (x-X)^2/bandwidth_X^2 + (y-Y)^2/bandwidth_Y^2 - (2*(x-X)*(y-Y)/bandwidth_X*bandwidth_Y);
pdf_compute = @(x,y,X,Y) (1/(2*pi()*bandwidth_X*bandwidth_Y)*exp(-z(x,y,X,Y)/2));
index = find(index~=0);
for j = 1 : length(index)
    for i = 1 : size(dataset,1)
        Dist(index(j)) = Dist(index(j))' + pdf_compute(X_true,values(index(j)),dataset(i,1),dataset(i,2));
    end
end

Dist = Dist./trapz(Y_comp,Dist);

end