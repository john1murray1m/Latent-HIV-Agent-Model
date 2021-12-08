%% import the LHS parameter set and run/analyse - include a flag (flag_array) which performs a single or array job
% set for array processing on a cluster



flag_array=1; % 1 if performing an array job where the index is set by the environment variable. The parameter 
              % values are contained in the files n.dat where n=myindex possible values. 
              % For the array jobs only output summary data to a csv file.
              % Otherwise 0 where output more detailed info to a .mat file for a single index
                % For the single index must specify the index and
                % parameters or call a particular parameter set.


%% timing - most of these data are set to emulate Patient ! of Maldarelli et al. Science 2014.
% CD4 count at start of ART was 22 cells so assume this is approximately 10
% years from infection. Time is in years.
 t0=10;
t_nonART=[0 0.2 t0/2 t0]; % this has been changed with the multiple calcs to keep L0

% update each timestep
timestep=1/10;


% run for t0 years without ART and then [t1, t2, t4] years with ART where
% samples stored at these time points to emulate the sample times of Patient 1.
% Also save the 7.2 year time point for pat 3 of Maldarelli and use the 0.2 time as surrogate for
% the actual preART sample at time 0 on ART.
t1=0.2;
t2=4.8;
t3=7.2;
t4=11.4;

t_ART=[t1 t2 t3 t4];
times=[t_nonART t0+t_ART];

%% use the following as an estimate of CD4 count prior to ART
bcd4=[387.6 -176.49];
cd4_nonArt=(bcd4(1)+bcd4(2)*log10(t_nonART(2:end)));

% replace the last  nonArt value with 22 copies to match Pat1 from
% Maldarelliand then the CD4 counts on ART
cd4_nonArt(end)=22;
cd4_Art=[48 185 190 193];

%% the sample sizes for patients 1 and 3
m_ART=[287 70 1314 1314];
m_ART3=[47 268 268 268];

%% whether to generate the proliferative values lambda using a Weibull or lognormal gmdistribution
% The simulations were perfored with a lognormal dist
flag_Weib=0; % so lognormal


%% The clones do not interact so the simulations of their proliferation/activation/death can be performed
% independently on subsets of clones. Split the entire clone size into subsets of size n1 
 n1=1e3; 
 
% how much starting value of agents is reduced to speed up calculation. 
% Some calculations for reservoir size dynamics are OK with a smaller total size in which case
% reduce the total size by fold_redn, eg set it to 100. For a full calculation this is 1
fold_redn=1; 

% the starting reservoir size (modified by fold_redn)
n0=3e8/fold_redn; % to start with. Must be multiple of n1

% Set the intact reservoir to be 100 fold lower than total
subn_intact=100; % the divisor for conversion of risk of intact versus defective
n0_intact=n0/subn_intact;

%% convert body measures to values per ml by assuming 5.5% of CD4 cells are in peripheral blood
% and there are 5 litres of blood in an individual
scale_ml=5.5/100/5000;

% random start for random numbers
rng shuffle

%% input the particular parameter set
 
 if flag_array
     myindex=getenv("PBS_ARRAY_INDEX");
    param0=readmatrix([num2str(myindex),'.dat']);
    
   %  specify a save file name for the summary data
   file_out=1; % whether to output the summary data to a file
   file_name='summary_output_';
  else
   myindex=1000;
   % param0=readmatrix([num2str(myindex),'.dat']); % if inputting from a file
   param0=[0.256805682,1.55169e-05,0.007735612,0.992550438,0.994661246,10,49,7]; 
  
   % also specify a save file name 
   file_name='single_output_';
 end


%% intact calculations
nki=n0_intact/n1; % this is the number of parallel calculations
nk=n0/n1; % this is the number of clones in each parallel calculation at the start

L_nonART_intact=cell(length(t_nonART),nki); % store the proliferative values for each clone calculated at the preART times
tt_nonART_intact=cell(length(t_nonART),nki); % store the times these clones arose
ncopy_nonART_intact=cell(length(t_nonART),nki); % store the number of copies for each clone - starts at 1

% now add in ART. Same as above. The only difference in the calculations will be no additional integrants
L_ART_intact=cell(length(t_ART),nki);
tt_ART_intact=cell(length(t_ART),nki);
ncopy_ART_intact=cell(length(t_ART),nki);

if flag_array               % Store the number of cells turning from latent to productive 
                            % from each of proliferation that fails to stay latent (I_lam_ART), or from the clone being
                            % activated (I_act_ART).Only do this during ART times since otherwise it is overshadowed from 
                            % ongoing productive infection. If calculating summary data (flag_array=1) then just store at 
                            % the end of ART and the 4 ART times. Otherwise store at each timestep.
  I_lam_ART=NaN(5,nki);
  I_act_ART=NaN(5,nki);
else
  I_lam_ART=NaN(round(t_ART(end)/timestep),nki);
  I_act_ART=NaN(round(t_ART(end)/timestep),nki);
end

%% Combine the parallel intact calculations
L_nonART_intact1=cell(length(t_nonART),1);
tt_nonART_intact1=cell(length(t_nonART),1);
ncopy_nonART_intact1=cell(length(t_nonART),1);

L_ART_intact1=cell(length(t_ART),1);
tt_ART_intact1=cell(length(t_ART),1);
ncopy_ART_intact1=cell(length(t_ART),1);
I_lam_ART1=NaN(size(I_lam_ART,1),1); % number of prod infected cells arising from proliferation
I_act_ART1=NaN(size(I_act_ART,1),1); % number of prod infected cells arising from activation

%% specify the parameter values
     
    mu=param0(1);
    lv=param0(2);

    % the lognormal values
    lm=1.*mu; % mean
    lmean=log(lm^2/sqrt(lv+lm^2));
    lsd=sqrt(log(lv/lm^2+1));

    alpha=param0(3);
    p_alpha=param0(4);
    p_alpha_intact=p_alpha;

    p_lam=param0(5);

    n_alpha=param0(6);  
    n_alpha_intact=n_alpha;

    subsource=param0(7);

    mact=n_alpha+param0(8);

    %% parameters for defective clones
    params(1)=mu;
    params(2)=lm;
    params(3)=lv;
    params(4)=alpha;
    params(5)=p_alpha;
    params(6)=n_alpha;
    params(7)=subsource;
    params(8)=mact;

    %% parameters for intact clones 
    params_intact(1)=mu;
    params_intact(2)=lm;
    params_intact(3)=lv;
    params_intact(4)=alpha;
    params_intact(5)=p_alpha_intact;
    params_intact(6)=n_alpha_intact;
    params_intact(7)=subsource;
    params_intact(8)=subn_intact;
    params_intact(9)=p_lam;
    params_intact(10)=mact;


    %% now run the intact calculation
    parfor ni=1:nki
        [L0,tt0,ncopy0]=f_nonART0_intact(n1,flag_Weib,params_intact); % first infection

        % infection over nonART subintervals
        [L1,tt1,ncopy1,I_lam1,I_act1]=f_nonART_act_intact(n1,timestep,flag_Weib,params_intact,L0,tt0,ncopy0,t_nonART(1),t_nonART(2));
        [L2,tt2,ncopy2,I_lam2,I_act2]=f_nonART_act_intact(n1,timestep,flag_Weib,params_intact,L1,tt1,ncopy1,t_nonART(2),t_nonART(3));
        [L3,tt3,ncopy3,I_lam3,I_act3]=f_nonART_act_intact(n1,timestep,flag_Weib,params_intact,L2,tt2,ncopy2,t_nonART(3),t_nonART(4));

    % infection over ART subintervals
        [LA1,ttA1,ncopyA1,I_lamA1,I_actA1]=f_ART_act_intact(timestep,params_intact,L3,tt3,ncopy3,t_nonART(end)+0,t_nonART(end)+t_ART(1));
        [LA2,ttA2,ncopyA2,I_lamA2,I_actA2]=f_ART_act_intact(timestep,params_intact,LA1,ttA1,ncopyA1,t_nonART(end)+t_ART(1),t_nonART(end)+t_ART(2));
        [LA3,ttA3,ncopyA3,I_lamA3,I_actA3]=f_ART_act_intact(timestep,params_intact,LA2,ttA2,ncopyA2,t_nonART(end)+t_ART(2),t_nonART(end)+t_ART(3));
        [LA4,ttA4,ncopyA4,I_lamA4,I_actA4]=f_ART_act_intact(timestep,params_intact,LA3,ttA3,ncopyA3,t_nonART(end)+t_ART(3),t_nonART(end)+t_ART(4));


        Lfred={L0,L1, L2, L3}';
        L_nonART_intact(:,ni)=Lfred;
        ttfred={tt0, tt1, tt2, tt3}';
        tt_nonART_intact(:,ni)=ttfred;
        nfred={ncopy0, ncopy1, ncopy2, ncopy3}';
        ncopy_nonART_intact(:,ni)=nfred;


        Lfred={LA1, LA2, LA3, LA4}';
        L_ART_intact(:,ni)=Lfred;
        ttfred={ttA1, ttA2, ttA3, ttA4}';
        tt_ART_intact(:,ni)=ttfred;
        nfred={ncopyA1, ncopyA2, ncopyA3, ncopyA4}';
        ncopy_ART_intact(:,ni)=nfred;
        if flag_array
            I_lam_ART(:,ni)=[I_lamA1([1 end]); I_lamA2(end); I_lamA3(end); I_lamA4(end)];
            I_act_ART(:,ni)=[I_actA1([1 end]); I_actA2(end); I_actA3(end); I_actA4(end)];
        else
            I_lam_ART(:,ni)=[I_lamA1; I_lamA2; I_lamA3; I_lamA4];
            I_act_ART(:,ni)=[I_actA1; I_actA2; I_actA3; I_actA4];
        end
    end
    
    % now combine into one 
    for j=1:length(t_nonART)
        L_nonART_intact1(j)={cell2mat(L_nonART_intact(j,:)')};
        ncopy_nonART_intact1(j)={cell2mat(ncopy_nonART_intact(j,:)')};
        tt_nonART_intact1(j)={cell2mat(tt_nonART_intact(j,:)')};
    end
    for j=1:length(t_ART)
        L_ART_intact1(j)={cell2mat(L_ART_intact(j,:)')};
        ncopy_ART_intact1(j)={cell2mat(ncopy_ART_intact(j,:)')};
        tt_ART_intact1(j)={cell2mat(tt_ART_intact(j,:)')};
    end

    I_lam_ART1=squeeze(nansum(I_lam_ART,2));
    I_act_ART1=squeeze(nansum(I_act_ART,2));

 %    disp('intact finished')
    
    %% now run the defective calculations
    parfor ni=1:nk
        [L0,tt0,ncopy0]=f_nonART0(n1,flag_Weib,params); % first infection

        % infection over nonART subintervals
        [L1,tt1,ncopy1]=f_nonART_act(n1,timestep,flag_Weib,params,L0,tt0,ncopy0,t_nonART(1),t_nonART(2));
        [L2,tt2,ncopy2]=f_nonART_act(n1,timestep,flag_Weib,params,L1,tt1,ncopy1,t_nonART(2),t_nonART(3));
        [L3,tt3,ncopy3]=f_nonART_act(n1,timestep,flag_Weib,params,L2,tt2,ncopy2,t_nonART(3),t_nonART(4));

    % infection over ART subintervals
        [LA1,ttA1,ncopyA1]=f_ART_act(timestep,params,L3,tt3,ncopy3,t_nonART(end)+0,t_nonART(end)+t_ART(1));
        [LA2,ttA2,ncopyA2]=f_ART_act(timestep,params,LA1,ttA1,ncopyA1,t_nonART(end)+t_ART(1),t_nonART(end)+t_ART(2));
        [LA3,ttA3,ncopyA3]=f_ART_act(timestep,params,LA2,ttA2,ncopyA2,t_nonART(end)+t_ART(2),t_nonART(end)+t_ART(3));
        [LA4,ttA4,ncopyA4]=f_ART_act(timestep,params,LA3,ttA3,ncopyA3,t_nonART(end)+t_ART(3),t_nonART(end)+t_ART(4));


        Lfred={L0,L1, L2, L3}';
        L_nonART(:,ni)=Lfred;
        ttfred={tt0, tt1, tt2, tt3}';
        tt_nonART(:,ni)=ttfred;
        nfred={ncopy0, ncopy1, ncopy2, ncopy3}';
        ncopy_nonART(:,ni)=nfred;

        Lfred={LA1, LA2, LA3, LA4}';
        L_ART(:,ni)=Lfred;
        ttfred={ttA1, ttA2, ttA3, ttA4}';
        tt_ART(:,ni)=ttfred;
        nfred={ncopyA1, ncopyA2, ncopyA3, ncopyA4}';
        ncopy_ART(:,ni)=nfred;

    end
    for j=1:length(t_nonART)
        L_nonART1(j)={cell2mat(L_nonART(j,:)')};
        ncopy_nonART1(j)={cell2mat(ncopy_nonART(j,:)')};
         tt_nonART1(j)={cell2mat(tt_nonART(j,:)')};
    end

    for j=1:(length(t_ART))
        L_ART1(j)={cell2mat(L_ART(j,:)')};
        ncopy_ART1(j)={cell2mat(ncopy_ART(j,:)')};
         tt_ART1(j)={cell2mat(tt_ART(j,:)')};
    end


if ~flag_array % the file will be about 5GB
    save([file_name,num2str(myindex)],'myindex','param0','fold_redn','L_nonART1','ncopy_nonART1','tt_nonART1',...
        'L_ART1','ncopy_ART1','tt_ART1','L_nonART_intact1','ncopy_nonART_intact1','tt_nonART_intact1',...
        'L_ART_intact1','ncopy_ART_intact1','tt_ART_intact1','I_lam_ART1','I_act_ART1','-v7.3')
else
    %% determine the clone profiles during ART - just for total
    clones_sampled=cell(length(t_ART),2); % collection of all clones sampled at each time
                % the first is for pat1 and second for pat3
    tbl_size=cell(length(t_ART),2); % table of clones sampled at each time
    out_sims=cell(7,1); % reservoir sizes, half-lives, cross-sampled clones for pat 1 then pat 3
                    % 5 is clone num, size info, 6 is single_perc info, 7
                    % is sampled single_perc info


for k=1:2
    if k==1 % sampling for patient 1, otherwise patient 3
        mmA=m_ART;
    else
        mmA=m_ART3;
    end
        
    for j=1:(length(t_ART))
        ndefect=length(ncopy_ART1{j}); % the number of defective clones
        nndefect=sum(ncopy_ART1{j}); % the number of defective integrants
        nintact=length(ncopy_ART_intact1{j}); % the number of intact clones
        nnintact=sum(ncopy_ART_intact1{j}); % the number of intact integrants

        nn=nndefect+nnintact;
            nsample=randsample(nn,mmA(j));
            nsample=sort(nsample);

            % determine the number of samples for each clone
            ncopy_sum=[0; cumsum([ncopy_ART1{j}; ncopy_ART_intact1{j}])];
             clone_num=NaN(mmA(j),1);
            lambda_num=NaN(mmA(j),1); % identify clones at different times through their
                                        % lambda value (which is a surrogate
                                        % for integration point; 
            fred=[L_ART1{j}; L_ART_intact1{j}];
            for i=1:mmA(j)
                 clone_num(i)=sum(nsample(i)>ncopy_sum);
                lambda_num(i)=fred(clone_num(i));
            end

            % total
            tbl0=tabulate(lambda_num);
            if ~isempty(tbl0)
                iij=tbl0(:,2)>0;
                tbl0=tbl0(iij,:);
                clones_sampled(j,k)={tbl0(:,1:2)}; % this is the clone identifier and 
                tbl=tabulate(tbl0(:,2));
                iij=tbl(:,2)>0;
                if ~isempty(iij)
                    tbl_size(j,k)={tbl(iij,1:2)};
                else
                    tbl_size(j,k)={[]};
                end
                    
            end
    end
end

    Res_size_defect=[cellfun(@sum, ncopy_nonART1) cellfun(@sum, ncopy_ART1)]*fold_redn;
    Res_size_intact=[cellfun(@sum, ncopy_nonART_intact1)' cellfun(@sum, ncopy_ART_intact1)']*fold_redn;
    Res_size=Res_size_defect+Res_size_intact;
    Res_size_ml=Res_size*scale_ml;
    Res_size_defect_ml=Res_size_defect*scale_ml;
    Res_size_intact_ml=Res_size_intact*scale_ml;
    Res_size_cd4=Res_size_ml/1e3./[1000 cd4_nonArt cd4_Art]*1e6;
    Res_size_defect_cd4=Res_size_defect_ml/1e3./[1000 cd4_nonArt cd4_Art]*1e6;
    Res_size_intact_cd4=Res_size_intact_ml/1e3./[1000 cd4_nonArt cd4_Art]*1e6;
    % estimate doubling time from the last 2 nonART times
    kgrowth=log(Res_size(4)/Res_size(3))/(times(4)-times(3));
    tgrowth=log(2)/kgrowth;

    kdecay=-log(Res_size_cd4(end-1)/Res_size_cd4(end-2))/(times(end-1)-times(end-2));
    thalf_cd4=log(2)/kdecay;
    kdecay=-log(Res_size_defect_cd4(end-1)/Res_size_defect_cd4(end-2))/(times(end-1)-times(end-2));
    thalf_defect_cd4=log(2)/kdecay;
    kdecay=-log(Res_size_intact_cd4(end-1)/Res_size_intact_cd4(end-2))/(times(end-1)-times(end-2));
    thalf_intact_cd4=log(2)/kdecay;

    %% now incorporate clone sizes - this will only be a rough estimate for fold_redn>1
    Clone_num_defect=[cellfun(@length, L_nonART1) cellfun(@length, L_ART1)]*fold_redn;
    Max_clone_defect=[cellfun(@max, ncopy_nonART1) cellfun(@max, ncopy_ART1)];
    Med_clone_defect=[cellfun(@median, ncopy_nonART1) cellfun(@median, ncopy_ART1)];

    Clone_num_intact=[cellfun(@length, L_nonART_intact1)' cellfun(@length, L_ART_intact1)']*fold_redn;
    Max_clone_intact=[cellfun(@max, ncopy_nonART_intact1)' cellfun(@max, ncopy_ART_intact1)'];
    Med_clone_intact=[cellfun(@median, ncopy_nonART_intact1)' cellfun(@median, ncopy_ART_intact1)'];

    Clone_num=Clone_num_defect+Clone_num_intact;
    Max_clone=nanmax([(Max_clone_defect); (Max_clone_intact)],[],1);

    out_sims(5)={[Clone_num_defect; Clone_num_intact; Max_clone_defect; Max_clone_intact; ...
                    Med_clone_defect; Med_clone_intact]};
                
                
%% The proportion of singletons in reservoir and in sample. 

% singleton clones as a fraction of all clones
Single_num_defect=[cellfun(@(x) sum(x==1), ncopy_nonART1) cellfun(@(x) sum(x==1), ncopy_ART1)]*fold_redn;
Single_num_intact=[cellfun(@(x) sum(x==1), ncopy_nonART_intact1)' cellfun(@(x) sum(x==1), ncopy_ART_intact1)']*fold_redn;
Single_num=Single_num_defect+Single_num_intact;
Single_perc=Single_num./Clone_num*100;
Single_perc_intact=Single_num_intact./Clone_num_intact*100;

% half life of singelton percentage of total clones
bsingle=regress(log(Single_perc)',[ones(size(times))' times']);
thalf_single_clone=-log(2)/bsingle(2);
iij=Single_perc_intact>0; % want to remove low numbers if not running large no. agents
bsingle_intact=regress(log(Single_perc_intact(iij))',[ones(size(times(iij)))' times(iij)']);
thalf_single_clone_intact=-log(2)/bsingle_intact(2);

% percentage singletons of sampled clones
% only from start of ART
istart=length(t_nonART)+1;
Single_perc_samples=cellfun(@(x) x(1,2),tbl_size(:,1))'./m_ART*100;


out_sims(6)={[Single_perc; Single_perc_intact]}; 
out_sims(7)={Single_perc_samples};

% store Res_size total and intact
out_sims(1)={[Res_size; Res_size_intact]};

% store growth and decay half-lives
out_sims(2)={[tgrowth, thalf_cd4, thalf_defect_cd4, thalf_intact_cd4]};

%     store the clones sampled at first time that appear in the last sample
for j=1:2 % 1 is pat 1 and 2 ispat3
    tbl=clones_sampled{1,j};
    tbl2=clones_sampled{5-j,j};
    [tbl7,ia,ib]=intersect(tbl(:,1),tbl2(:,1));
    if ~isempty(ia)
        george=[tbl(ia,2) tbl2(ib,2)]; % the clone sizes of rpts in each time
        tbl3=tabulate(tbl(ia,2));
        iij=tbl3(:,2)>0;
        tbl3=tbl3(iij,:);
        gu=unique(george(:,1));
        ng=length(gu);
        ggu=NaN(2,ng);
        for k=1:ng
            iij=george(:,1)==gu(k);
            ggu(1,k)=min(george(iij,2));
            ggu(2,k)=max(george(iij,2));
        end
        % the size of a cross-sampled clone in the first sample (gu), how
        % many there were of this size (tbl3) and the min and max sizes
        % when sampled at the last time
        out_sims(2+j)={[gu tbl3(:,2) ggu']};
    else
        out_sims(2+j)={[]};
    end
end
    
file_out=1; % if outputting these summary data
    if file_out
        out_filename=[file_name,num2str(myindex),'.csv'];
        writematrix(param0,out_filename)
        writematrix(out_sims{1},out_filename,'WriteMode','append')
         writematrix(out_sims{2},out_filename,'WriteMode','append')
         writematrix(out_sims{5},out_filename,'WriteMode','append')
         writematrix(out_sims{6},out_filename,'WriteMode','append')
         writematrix(out_sims{7},out_filename,'WriteMode','append')
         writematrix(I_lam_ART1',out_filename,'WriteMode','append')
         writematrix(I_act_ART1',out_filename,'WriteMode','append')
         writematrix('ART 0.2 pat1',out_filename,'WriteMode','append')
         writematrix(tbl_size{1,1},out_filename,'WriteMode','append')
         writematrix('ART 4.8',out_filename,'WriteMode','append')
         writematrix(tbl_size{2,1},out_filename,'WriteMode','append')
         writematrix('ART 11.4',out_filename,'WriteMode','append')
         writematrix(tbl_size{4,1},out_filename,'WriteMode','append')
         writematrix('ART 0.2 & 11.4',out_filename,'WriteMode','append')
         writematrix(out_sims{3},out_filename,'WriteMode','append')
         writematrix('ART 0.2 pat3',out_filename,'WriteMode','append')
         writematrix(tbl_size{1,2},out_filename,'WriteMode','append')
         writematrix('ART 7.2 pat3',out_filename,'WriteMode','append')
         writematrix(tbl_size{3,2},out_filename,'WriteMode','append')
         writematrix('ART 0.2 & 7.2',out_filename,'WriteMode','append')
         writematrix(out_sims{4},out_filename,'WriteMode','append')
    end
end 

