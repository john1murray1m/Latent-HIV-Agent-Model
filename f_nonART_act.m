%% generate the integrants in the absence of ART, for each time period
function [L,tt,ncopy]=f_nonART_act(n1,timestep,flag_Weib,params,L0,tt0,ncopy0,ta,tb)

mu=params(1);
lm=params(2);
lv=params(3);
alpha=params(4);
p_alpha=params(5);
n_alpha=params(6);
subsource=params(7);
mact=params(8);

s=n1/subsource; % number of new integrants per year

if flag_Weib
    pdweib=makedist('Weibull','a',lm,'b',lv);
else
    lmean=log(lm^2/sqrt(lv+lm^2));
    lsd=sqrt(log(lv/lm^2+1));
end

L=L0;
ncopy=ncopy0;
tt=tt0;
tdiff=tb-ta; % number of years in this group

    for itime=1:round(tdiff/timestep) % number of timesteps. 
        %% death. 
        % Cells die off at rate mu and proportional to copy number
        ncopy1=binornd(ncopy,mu*timestep);
        ncopy=ncopy-ncopy1;
        iij=ncopy>0; % these are still alive

        % update
        L=L(iij);
        ncopy=ncopy(iij);
        tt=tt(iij);

        %% proliferation
        ncopy1=binornd(ncopy,L*timestep);
        ncopy=ncopy+ncopy1;

        %% activation
        lsj=length(ncopy);
        bG=ones(lsj,1);
        
        ncopy1=binornd(bG,alpha*timestep);
        % choose those that have activated and then the fraction who will
        % proceed to expand
        iij=ncopy1==0; % these are not activated
        
        nfred=ncopy(iij); % not activated
        Lfred=L(iij);
        tfred=tt(iij); 

        % those that are activated
        lsj=sum(~iij);
        if lsj>0 % there are some that are activated

            ngeorge=ncopy(~iij); %  activated
            Lgeorge=L(~iij);
            tgeorge=tt(~iij); 
            
            bG=ones(lsj,1);
            nact=binornd(bG,1-p_alpha);
            iij1=nact==1; % those that survive activation
            
            % these do not get expanded but some might survive if
            % ncopy>2^mact
            nnact=sum(~iij1); 
            nharry=ngeorge(~iij1)-2^mact;
            Lharry=Lgeorge(~iij1);
            tharry=tgeorge(~iij1);
            
           %remove those that are <=zero
           iij3=nharry>0;
           nharry=nharry(iij3);
           Lharry=Lharry(iij3);
           tharry=tharry(iij3);
            
            niij1=sum(iij1);
            if niij1>0 % some expand
                nruthy=max(ngeorge(iij1)-2^mact,zeros(niij1,1));
                Lruthy=Lgeorge(iij1);
                truthy=tgeorge(iij1);

                    ndiv=poissrnd(n_alpha,niij1,1);
                    % Assume they expand to numbers independent
                    % of starting size and growth rate L;
                    % the remainder are removed

                    nruthy=nruthy+2.^ndiv;
                iij4=nruthy>0;
                nruthy=nruthy(iij4);
                Lruthy=Lruthy(iij4);
                truthy=truthy(iij4);
                 
                ncopy=[nfred; nharry; nruthy];
                L=[Lfred; Lharry; Lruthy];
                tt=[tfred; tharry; truthy];
                

             else % none expand
                ncopy=[nfred; nharry];
                L=[Lfred; Lharry];
                tt=[tfred; tharry];

            end
        else % none expand
                ncopy=nfred;
                L=Lfred;
                tt=tfred; 

        end

        % now add in new integrants
        if flag_Weib
            Lnew=random(pdweib,round(s*timestep),1); % Weibull dist
        else
            Lnew=lognrnd(lmean,lsd,round(s*timestep),1); % lognormal dist
        end

        iij=Lnew>0;
        Lnew=Lnew(iij);
        ttnew=itime*ones(size(Lnew))*timestep+ta;
        ncopynew=ones(size(Lnew));

        % update the total
        L=[L; Lnew];
        tt=[tt; ttnew];
        ncopy=[ncopy; ncopynew];

    end
