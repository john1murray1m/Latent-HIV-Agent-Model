%% generate the intact integrants in the  beginning
function [L,tt,ncopy]=f_nonART0_intact(n1,flag_Weib,params_intact)

lm=params_intact(2);
lv=params_intact(3);

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
