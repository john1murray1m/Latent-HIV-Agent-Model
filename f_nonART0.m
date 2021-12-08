%% generate the defective integrants at the beginning
function [L,tt,ncopy]=f_nonART0(n1,flag_Weib,params)

mu=params(1);
lm=params(2);
lv=params(3);
alpha=params(4);
p_alpha=params(5);
n_alpha=params(6);
subsource=params(7);

if flag_Weib
    pdweib=makedist('Weibull','a',lm,'b',lv);
    L0=random(pdweib,n1,1); % Weibull dist
else
    lmean=log(lm^2/sqrt(lv+lm^2));
    lsd=sqrt(log(lv/lm^2+1));
    L0=lognrnd(lmean,lsd,n1,1); % lognormal dist
end
iij=L0>0;
L=L0(iij);

% their time of calculation
tt=zeros(size(L0));

% their copy number
ncopy=ones(size(L0));
