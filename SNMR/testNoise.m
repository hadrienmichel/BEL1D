
function [isNoisy, Std_CCA] = testNoise(Prior_d, PCA_d_base, A, default_bandwidth, noise_level)

% Getting the number of PCA dim kept:
dimD = size(PCA_d_base,2);
dimH = size(A,2);

% Testing the impact of noise (on 50 random samples):
nbTest = 50;
Cov_df = zeros(nbTest,dimD,dimD);
index = ceil(rand(1,nbTest).*size(Prior_d,1));
warning('off','all');
h = waitbar(0,'Please wait . . .');
tic;
for i = 1 : nbTest,
    data_real = Prior_d(index(i),:)+randn(1,size(Prior_d,2))*noise_level;
    data_synthe = Prior_d;
    data_synthe(index(i),:) = [];
    [coeff,~,~,~,~,~]=pca(data_synthe,'NumComponents',dimD);
    score_D_tmp = (data_real - mean(data_synthe,1))*coeff;
    score_D_tmp = score_D_tmp(1:dimD);
    Cov_df(i,:,:) = cov([PCA_d_base(index(i),:);score_D_tmp]);
    time = toc;
    time = time/i;% Time per pass
    remain = time*(nbTest-i);
    waitbar(i/nbTest,h,['Remaining time: ' num2str(remain) ' sec']);
end
warning('on','all');

C_f = squeeze(mean(Cov_df,1));

figure
image(C_f,'CDataMapping','scaled')
map = colormap(jet);
%map(1,:) = 1;
colormap(map);
colorbar('eastoutside');
axis([0.5 dimD+0.5 0.5 dimD+0.5])
xlabel('Component j','Fontsize',16);ylabel('Component i','Fontsize',16);
title('Mean Covariance matrix C^f','FontSize',16)
set(gca,'FontSize',16)

C_c = A'*C_f*A;

figure
image(C_c,'CDataMapping','scaled')
map = colormap(jet);
%map(1,:) = 1;
colormap(map);
colorbar('eastoutside');
axis([0.5 dimH+0.5 0.5 dimH+0.5])
xlabel('Component j','Fontsize',16);ylabel('Component i','Fontsize',16);
title('Mean Covariance matrix C^c','FontSize',16)
set(gca,'FontSize',16)

std_default = ones(dimH,1).*default_bandwidth;
if max(diag(C_c))>(default_bandwidth^2),
    isNoisy = true;
    Std_CCA = diag(C_c).^(1/2);% We are neglecting the components outside the diagonal
    Std_CCA = max(Std_CCA, std_default);
else
    isNoisy = false;
    Std_CCA = std_default;
end

delete(h);

end