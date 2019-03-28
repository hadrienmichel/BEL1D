function [models] = ModelGenerator(type, N, parameters, nb_layer)

    thick = parameters(1:end-1,1:2);
    param2 = parameters(:,3:4);
    param3 = parameters(:,5:6);
    param4 = parameters(:,7:8);
    OK = false;
    passed = false;
    while ~OK,
        models = struct('thick',zeros(N,nb_layer),'param2',zeros(N,nb_layer),'param3',zeros(N,nb_layer),'param4',zeros(N,nb_layer));
        if type == 1,% Uniformly distributed
            % Computing parameters
            models.thick = ones(N,nb_layer-1)*diag(thick(:,1)) + rand(N,nb_layer-1)*(diag(thick(:,2)-thick(:,1)));
            models.param2 = ones(N,nb_layer) * diag(param2(:,1)) + rand(N,nb_layer)*(diag(param2(:,2)-param2(:,1)));
            models.param3 = ones(N,nb_layer)*diag(param3(:,1)) + rand(N,nb_layer)*(diag(param3(:,2)-param3(:,1)));
            models.param4 = ones(N,nb_layer)*diag(param4(:,1)) + rand(N,nb_layer)*(diag(param4(:,2)-param4(:,1)));
        elseif type == 2,% Latin Hypercube Sampler uniformly distributed
            % Computing parameters
            temporary = lhsdesign(N,nb_layer*4-1); % Latin hypercube sampling
            models.thick = ones(N,nb_layer-1)*diag(thick(:,1)) + temporary(:,1:nb_layer-1)*(diag(thick(:,2)-thick(:,1)));
            models.param2 = ones(N,nb_layer) * diag(param2(:,1)) + temporary(:,nb_layer:(2*nb_layer)-1)*(diag(param2(:,2)-param2(:,1)));
            models.param3 = ones(N,nb_layer)*diag(param3(:,1)) + temporary(:,2*nb_layer:(3*nb_layer)-1)*(diag(param3(:,2)-param3(:,1)));
            models.param4 = ones(N,nb_layer)*diag(param4(:,1)) + temporary(:,3*nb_layer:end)*(diag(param4(:,2)-param4(:,1)));
        elseif type == 3,% Normal distribution (col1 = mu, col2 = sigma)
            tmp1 = (thick(:,1)+thick(:,2))./2;
            tmp2 = abs(thick(:,2)-thick(:,1))./4;
            thick = [tmp1 tmp2];
            tmp1 = (param2(:,1)+param2(:,2))./2;
            tmp2 = abs(param2(:,2)-param2(:,1))./4;
            param2 = [tmp1 tmp2];
            tmp1 = (param3(:,1)+param3(:,2))./2;
            tmp2 = abs(param3(:,2)-param3(:,1))./4;
            param3 = [tmp1 tmp2];
            tmp1 = (param4(:,1)+param4(:,2))./2;
            tmp2 = abs(param4(:,2)-param4(:,1))./4;
            param4 = [tmp1 tmp2];
            % Computing parameters
            models.thick = normrnd(repmat(thick(:,1)',N,1),repmat(thick(:,2)',N,1));
            models.param2 = normrnd(repmat(param2(:,1)',N,1),repmat(param2(:,2)',N,1));
            models.param3(:,:) = normrnd(repmat(param3(:,1)',N,1),repmat(param3(:,2)',N,1));
            models.param4(:,:) = normrnd(repmat(param4(:,1)',N,1),repmat(param4(:,2)',N,1));
            models.thick(models.thick<0) = 0;
            models.param2(models.param2<0) = 0;
            models.param3(models.param3<0) = 0;
            models.param4(models.param4<0) = 0;
        elseif type == 4,% Latin hypecube sampler Normal distribution (col1 = mu, col2 = sigma)
            tmp1 = (thick(:,1)+thick(:,2))./2;
            tmp2 = abs(thick(:,2)-thick(:,1))./4;
            thick = [tmp1 tmp2];
            tmp1 = (param2(:,1)+param2(:,2))./2;
            tmp2 = abs(param2(:,2)-param2(:,1))./4;
            param2 = [tmp1 tmp2];
            tmp1 = (param3(:,1)+param3(:,2))./2;
            tmp2 = abs(param3(:,2)-param3(:,1))./4;
            param3 = [tmp1 tmp2];
            tmp1 = (param4(:,1)+param4(:,2))./2;
            tmp2 = abs(param4(:,2)-param4(:,1))./4;
            param4 = [tmp1 tmp2];
            % Computing parameters
            temporary = lhsnorm([thick(:,1)', param2(:,1)', param3(:,1)', param4(:,1)'],diag([thick(:,2)', param2(:,2)', param3(:,2)'], param4(:,2)),N);
            models.thick = temporary(:,1:nb_layer-1);
            models.param2 = temporary(:,nb_layer:(nb_layer*2)-1);
            models.param3 = temporary(:,nb_layer*2:(nb_layers*3)-1);
            models.param4 = temporary(:,nb_layer*3:end);
            models.thick(models.thick<0) = 0;
            models.param2(models.param2<0) = 0;
            models.param3(models.param3<0) = 0;
            models.param4(models.param4<0) = 0;
        else
            error('You didn''t choose a satistical law');
        end
        %% Add here your constraints on the sampled models (/!\ may not work with latin-hypercube sampler), else put OK = true;
%         if passed
%             tmp = models.thick+models.param2(:,1);
%             N = N - sum(tmp>=10);
%             models_save.thick = [models_save.thick; models.thick(tmp>=10)];
%             models_save.param2 = [models_save.param2; models.param2(tmp>=10,:)];
%         else
%             tmp = models.thick+models.param2(:,1);
%             N = N - sum(tmp>=10);
%             models_save = models;
%             models_save.thick = models_save.thick(tmp>=10);
%             models_save.param2 = models_save.param2(tmp>=10,:);
%         end
%         passed = true;
%         if N == 0
%             OK = true;
%         end
        %% No conditions
        models_save = models;
        OK = true;
    end
    models = models_save;

end