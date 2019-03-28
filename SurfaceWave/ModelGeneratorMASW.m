function [models] = ModelGeneratorMASW(type, N, parameters, nb_layer, increasing)

    models = struct('thick',zeros(N,nb_layer-1),'Vp',zeros(N,nb_layer),'Vs',zeros(N,nb_layer),'rho',zeros(N,nb_layer));
    
    thick = parameters(1:end-1,1:2); % Thickness
    Vp = parameters(:,3:4); % P-Wave velocity
    Nu = parameters(:,5:6);% Poisson ratio (if NaNs, not used).
    Vs = parameters(:,7:8);% S-Wave velocity
    rho = parameters(:,9:10);% Density

    if ~increasing,
        %% For purely random values
        if type == 1,% Uniformly distributed
            % Computing parameters
            models.thick = ones(N,nb_layer-1)*diag(thick(:,1)) + rand(N,nb_layer-1)*(diag(thick(:,2)-thick(:,1)));
            models.Vp = ones(N,nb_layer) * diag(Vp(:,1)) + rand(N,nb_layer)*(diag(Vp(:,2)-Vp(:,1)));
            if ~isnan(Nu(1,1)),
                Nu_tmp = ones(N,nb_layer) * diag(Nu(:,1)) + rand(N,nb_layer)*(diag(Nu(:,2)-Nu(:,1)));
                models.Vs = models.Vp.*(((ones(N,nb_layer)-Nu_tmp)./((ones(N,nb_layer).*0.5)-Nu_tmp)).^(-0.5));
            else
                models.Vs = ones(N,nb_layer)*diag(Vs(:,1)) + rand(N,nb_layer)*(diag(Vs(:,2)-Vs(:,1)));
            end
            models.rho = ones(N,nb_layer)*diag(rho(:,1)) + rand(N,nb_layer)*(diag(rho(:,2)-rho(:,1)));
        elseif type == 2,% Latin Hypercube Sampler uniformly distributed
            % Computing parameters
            temporary = lhsdesign(N,nb_layer*4-1); % Latin hypercube sampling
            models.thick = ones(N,nb_layer-1)*diag(thick(:,1)) + temporary(:,1:nb_layer-1)*(diag(thick(:,2)-thick(:,1)));
            models.Vp = ones(N,nb_layer) * diag(Vp(:,1)) + temporary(:,nb_layer:(2*nb_layer)-1)*(diag(Vp(:,2)-Vp(:,1)));
            if ~isnan(Nu(1,1)),
                Nu_tmp = ones(N,nb_layer)*diag(Nu(:,1)) + temporary(:,2*nb_layer:(3*nb_layer)-1)*(diag(Nu(:,2)-Nu(:,1)));
                models.Vs = models.Vp.*(((ones(N,nb_layer)-Nu_tmp)/((ones(N,nb_layer).*0.5)-Nu_tmp)).^(-0.5));
            else
                models.Vs = ones(N,nb_layer)*diag(Vs(:,1)) + temporary(:,2*nb_layer:(3*nb_layer)-1)*(diag(Vs(:,2)-Vs(:,1)));
            end
            models.Vs = ones(N,nb_layer)*diag(Vs(:,1)) + temporary(:,2*nb_layer:(3*nb_layer)-1)*(diag(Vs(:,2)-Vs(:,1)));
            models.rho = ones(N,nb_layer)*diag(rho(:,1)) + temporary(:,3*nb_layer:end)*(diag(rho(:,2)-rho(:,1)));
        elseif type == 3,% Normal distribution (col1 = mu, col2 = sigma)
            tmp1 = (thick(:,1)+thick(:,2))./2;
            tmp2 = abs(thick(:,2)-thick(:,1))./4;
            thick = [tmp1 tmp2];
            tmp1 = (Vp(:,1)+Vp(:,2))./2;
            tmp2 = abs(Vp(:,2)-Vp(:,1))./4;
            Vp = [tmp1 tmp2];
            tmp1 = (Nu(:,1)+Nu(:,2))./2;
            tmp2 = abs(Nu(:,2)-Nu(:,1))./4;
            Nu = [tmp1 tmp2];
            tmp1 = (Vs(:,1)+Vs(:,2))./2;
            tmp2 = abs(Vs(:,2)-Vs(:,1))./4;
            Vs = [tmp1 tmp2];
            tmp1 = (rho(:,1)+rho(:,2))./2;
            tmp2 = abs(rho(:,2)-rho(:,1))./4;
            rho = [tmp1 tmp2];
            % Computing parameters
            models.thick = normrnd(repmat(thick(:,1)',N,1),repmat(thick(:,2)',N,1));
            models.Vp = normrnd(repmat(Vp(:,1)',N,1),repmat(Vp(:,2)',N,1));
            if ~isnan(Nu(1,1)),
                Nu_tmp = normrnd(repmat(Nu(:,1)',N,1),repmat(Nu(:,2)',N,1));
                models.Vs = models.Vp.*(((ones(N,nb_layer)-Nu_tmp)/((ones(N,nb_layer).*0.5)-Nu_tmp)).^(-0.5));
            else
                models.Vs(:,:) = normrnd(repmat(Vs(:,1)',N,1),repmat(Vs(:,2)',N,1));
            end
            models.rho(:,:) = normrnd(repmat(rho(:,1)',N,1),repmat(rho(:,2)',N,1));
            models.thick(models.thick<0) = 0;
            models.Vp(models.Vp<0) = 0;
            models.Vs(models.Vs<0) = 0;
            models.rho(models.rho<0) = 0;
        elseif type == 4,% Latin hypecube sampler Normal distribution (col1 = mu, col2 = sigma)
            tmp1 = (thick(:,1)+thick(:,2))./2;
            tmp2 = abs(thick(:,2)-thick(:,1))./4;
            thick = [tmp1 tmp2];
            tmp1 = (Vp(:,1)+Vp(:,2))./2;
            tmp2 = abs(Vp(:,2)-Vp(:,1))./4;
            Vp = [tmp1 tmp2];
            tmp1 = (Nu(:,1)+Nu(:,2))./2;
            tmp2 = abs(Nu(:,2)-Nu(:,1))./4;
            Nu = [tmp1 tmp2];
            tmp1 = (Vs(:,1)+Vs(:,2))./2;
            tmp2 = abs(Vs(:,2)-Vs(:,1))./4;
            Vs = [tmp1 tmp2];
            tmp1 = (rho(:,1)+rho(:,2))./2;
            tmp2 = abs(rho(:,2)-rho(:,1))./4;
            rho = [tmp1 tmp2];
            % Computing parameters
            temporary = lhsnorm([thick(:,1)', Vp(:,1)', Vs(:,1)', rho(:,1)'],diag([thick(:,2)', Vp(:,2)', Vs(:,2)'], rho(:,2)),N);
            models.thick = temporary(:,1:nb_layer-1);
            models.Vp = temporary(:,nb_layer:(nb_layer*2)-1);
            if ~isnan(Nu(1,1)),
                Nu_tmp = temporary(:,nb_layer*2:(nb_layers*3)-1);
                models.Vs = models.Vp.*(((ones(N,nb_layer)-Nu_tmp)/((ones(N,nb_layer).*0.5)-Nu_tmp)).^(-0.5));
            else
                models.Vs = temporary(:,nb_layer*2:(nb_layers*3)-1);
            end
            models.rho = temporary(:,nb_layer*3:end);
            models.thick(models.thick<0) = 0;
            models.Vp(models.Vp<0) = 0;
            models.Vs(models.Vs<0) = 0;
            models.rho(models.rho<0) = 0;
        else
            error('You didn''t choose a satistical law');
        end
    else
        %% For increaing values (no Latin-Hypercube possible!!!)
        if type == 1,% Uniformly distributed
            % Computing parameters
            NOTOK = true;
            index_change = true(N,1);
            while NOTOK,
                models.thick(index_change,:) = ones(sum(index_change),nb_layer-1)*diag(thick(:,1)) + rand(sum(index_change),nb_layer-1)*(diag(thick(:,2)-thick(:,1)));
                models.Vp(index_change,:) = ones(sum(index_change),nb_layer) * diag(Vp(:,1)) + rand(sum(index_change),nb_layer)*(diag(Vp(:,2)-Vp(:,1)));
                if ~isnan(Nu(1,1)),
                    Nu_tmp = ones(sum(index_change),nb_layer) * diag(Nu(:,1)) + rand(sum(index_change),nb_layer)*(diag(Nu(:,2)-Nu(:,1)));
                    models.Vs(index_change,:) = models.Vp(index_change,:).*(((ones(sum(index_change),nb_layer)-Nu_tmp)./((ones(sum(index_change),nb_layer).*0.5)-Nu_tmp)).^(-0.5));
                else
                    models.Vs(index_change,:) = ones(sum(index_change),nb_layer)*diag(Vs(:,1)) + rand(sum(index_change),nb_layer)*(diag(Vs(:,2)-Vs(:,1)));
                end
                models.rho(index_change,:) = ones(sum(index_change),nb_layer)*diag(rho(:,1)) + rand(sum(index_change),nb_layer)*(diag(rho(:,2)-rho(:,1)));
                index_change = false(N,1);
                for i = 1 : nb_layer,
                    if i < nb_layer,
                        index_change = (index_change | models.Vp(:,i)>models.Vp(:,i+1));
                        index_change = (index_change | models.Vs(:,i)>models.Vs(:,i+1));
                        index_change = (index_change | models.rho(:,i)>models.rho(:,i+1));
                    end
                    index_change = (index_change | models.Vs(:,i)<Vs(i,1) | models.Vs(:,i)>Vs(i,2));
                end
                if sum(index_change) == 0,
                    NOTOK = false;
                end
            end
        elseif type == 2,% Latin Hypercube Sampler uniformly distributed
            % Computing parameters
            warning(sprintf('You are attempting to use a Latin-Hypercube sampler on a rejection sampler.\nThe sampler has been changed to the non latin hypercube equivalent!\n'));
            NOTOK = true;
            index_change = true(N,1);
            while NOTOK,
                models.thick(index_change,:) = ones(sum(index_change),nb_layer-1)*diag(thick(:,1)) + rand(sum(index_change),nb_layer-1)*(diag(thick(:,2)-thick(:,1)));
                models.Vp(index_change,:) = ones(sum(index_change),nb_layer) * diag(Vp(:,1)) + rand(sum(index_change),nb_layer)*(diag(Vp(:,2)-Vp(:,1)));
                if ~isnan(Nu(1,1)),
                    Nu_tmp = ones(sum(index_change),nb_layer) * diag(Nu(:,1)) + rand(sum(index_change),nb_layer)*(diag(Nu(:,2)-Nu(:,1)));
                    models.Vs(index_change,:) = models.Vp(index_change,:).*(((ones(sum(index_change),nb_layer)-Nu_tmp)./((ones(sum(index_change),nb_layer).*0.5)-Nu_tmp)).^(-0.5));
                else
                    models.Vs(index_change,:) = ones(sum(index_change),nb_layer)*diag(Vs(:,1)) + rand(sum(index_change),nb_layer)*(diag(Vs(:,2)-Vs(:,1)));
                end
                models.rho(index_change,:) = ones(sum(index_change),nb_layer)*diag(rho(:,1)) + rand(sum(index_change),nb_layer)*(diag(rho(:,2)-rho(:,1)));
                index_change = false(N,1);
                for i = 1 : nb_layer-1,
                    if i < nb_layer,
                        index_change = (index_change | models.Vp(:,i)>models.Vp(:,i+1));
                        index_change = (index_change | models.Vs(:,i)>models.Vs(:,i+1));
                        index_change = (index_change | models.rho(:,i)>models.rho(:,i+1));
                    end
                    index_change = (index_change | models.Vs(:,i)<Vs(i,1) | models.Vs(:,i)>Vs(i,2));
                end
                if sum(index_change) == 0,
                    NOTOK = false;
                end
            end
        elseif type == 3,% Normal distribution (col1 = mu, col2 = sigma)
            tmp1 = (thick(:,1)+thick(:,2))./2;
            tmp2 = abs(thick(:,2)-thick(:,1))./4;
            thick = [tmp1 tmp2];
            tmp1 = (Vp(:,1)+Vp(:,2))./2;
            tmp2 = abs(Vp(:,2)-Vp(:,1))./4;
            Vp = [tmp1 tmp2];
            tmp1 = (Nu(:,1)+Nu(:,2))./2;
            tmp2 = abs(Nu(:,2)-Nu(:,1))./4;
            Nu = [tmp1 tmp2];
            tmp1 = (Vs(:,1)+Vs(:,2))./2;
            tmp2 = abs(Vs(:,2)-Vs(:,1))./4;
            Vs = [tmp1 tmp2];
            tmp1 = (rho(:,1)+rho(:,2))./2;
            tmp2 = abs(rho(:,2)-rho(:,1))./4;
            rho = [tmp1 tmp2];
            % Computing parameters
            NOTOK = true;
            index_change = true(N,1);
            while NOTOK,
                models.thick(index_change,:) = normrnd(repmat(thick(:,1)',sum(index_change),1),repmat(thick(:,2)',sum(index_change),1));
                models.Vp(index_change,:) = normrnd(repmat(Vp(:,1)',sum(index_change),1),repmat(Vp(:,2)',sum(index_change),1));
                if ~isnan(Nu(1,1)),
                    Nu_tmp = normrnd(repmat(Nu(:,1)',sum(index_change),1),repmat(Nu(:,2)',sum(index_change),1));
                    models.Vs(index_change,:) = models.Vp(index_change,:).*(((ones(sum(index_change),nb_layer)-Nu_tmp)./((ones(sum(index_change),nb_layer).*0.5)-Nu_tmp)).^(-0.5));
                else
                    models.Vs(index_change,:) = normrnd(repmat(Vs(:,1)',sum(index_change),1),repmat(Vs(:,2)',sum(index_change),1));
                end
                models.rho(index_change,:) = normrnd(repmat(rho(:,1)',sum(index_change),1),repmat(rho(:,2)',sum(index_change),1));
                models.thick(models.thick<0) = 0;
                models.Vp(models.Vp<0) = 0;
                models.Vs(models.Vs<0) = 0;
                models.rho(models.rho<0) = 0;
                index_change = false(N,1);
                for i = 1 : nb_layer-1,
                    if i < nb_layer,
                        index_change = (index_change | models.Vp(:,i)>models.Vp(:,i+1));
                        index_change = (index_change | models.Vs(:,i)>models.Vs(:,i+1));
                        index_change = (index_change | models.rho(:,i)>models.rho(:,i+1));
                    end
                    index_change = (index_change | models.Vs(:,i)<Vs(i,1) | models.Vs(:,i)>Vs(i,2));
                end
                if sum(index_change) == 0,
                    NOTOK = false;
                end
            end
        elseif type == 4,% Latin hypecube sampler Normal distribution (col1 = mu, col2 = sigma)
            warning(sprintf('You are attempting to use a Latin-Hypercube sampler on a rejection sampler.\nThe sampler has been changed to the non latin hypercube equivalent!\n'));
            tmp1 = (thick(:,1)+thick(:,2))./2;
            tmp2 = abs(thick(:,2)-thick(:,1))./4;
            thick = [tmp1 tmp2];
            tmp1 = (Vp(:,1)+Vp(:,2))./2;
            tmp2 = abs(Vp(:,2)-Vp(:,1))./4;
            Vp = [tmp1 tmp2];
            tmp1 = (Nu(:,1)+Nu(:,2))./2;
            tmp2 = abs(Nu(:,2)-Nu(:,1))./4;
            Nu = [tmp1 tmp2];
            tmp1 = (Vs(:,1)+Vs(:,2))./2;
            tmp2 = abs(Vs(:,2)-Vs(:,1))./4;
            Vs = [tmp1 tmp2];
            tmp1 = (rho(:,1)+rho(:,2))./2;
            tmp2 = abs(rho(:,2)-rho(:,1))./4;
            rho = [tmp1 tmp2];
            % Computing parameters
            NOTOK = true;
            index_change = true(N,1);
            while NOTOK,
                models.thick(index_change,:) = normrnd(repmat(thick(:,1)',sum(index_change),1),repmat(thick(:,2)',sum(index_change),1));
                models.Vp(index_change,:) = normrnd(repmat(Vp(:,1)',sum(index_change),1),repmat(Vp(:,2)',sum(index_change),1));
                if ~isnan(Nu(1,1)),
                    Nu_tmp = normrnd(repmat(Nu(:,1)',sum(index_change),1),repmat(Nu(:,2)',sum(index_change),1));
                    models.Vs(index_change,:) = models.Vp(index_change,:).*(((ones(sum(index_change),nb_layer)-Nu_tmp)./((ones(sum(index_change),nb_layer).*0.5)-Nu_tmp)).^(-0.5));
                else
                    models.Vs(index_change,:) = normrnd(repmat(Vs(:,1)',sum(index_change),1),repmat(Vs(:,2)',sum(index_change),1));
                end
                models.rho(index_change,:) = normrnd(repmat(rho(:,1)',sum(index_change),1),repmat(rho(:,2)',sum(index_change),1));
                models.thick(models.thick<0) = 0;
                models.Vp(models.Vp<0) = 0;
                models.Vs(models.Vs<0) = 0;
                models.rho(models.rho<0) = 0;
                index_change = false(N,1);
                for i = 1 : nb_layer-1,
                    if i < nb_layer,
                        index_change = (index_change | models.Vp(:,i)>models.Vp(:,i+1));
                        index_change = (index_change | models.Vs(:,i)>models.Vs(:,i+1));
                        index_change = (index_change | models.rho(:,i)>models.rho(:,i+1));
                    end
                    index_change = (index_change | models.Vs(:,i)<Vs(i,1) | models.Vs(:,i)>Vs(i,2));
                end
                remains = sum(index_change)
                if sum(index_change) == 0,
                    NOTOK = false;
                end
            end
        else
            error('You didn''t choose a satistical law');
        end
    end

end