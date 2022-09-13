% May 1, 2020: revising code by S. Gavrilets on evolution of conformity 
% First draft by S. Gavrilets, modifications by A. Jacquet
% There are 2 diallelic genes one of which with alleles A,a is expressed in males 
% only and another with alleles B, b in females only; the recombination rate is r.
% The frequencies of alleles A and B are p and q, respectively. Male trait
% is neutral. Females b mate randomly. 'Copier' females B prefer males that had high mating
% success in the previous generation. I am using an S-shaped "conformity
% function" in which probablity of mating with males of a paticular type is
% p^beta/(p^beta+(1-p)^beta), where p is their mating success in the previous 
% generation and beta is the steepness parameter.
%
% The population size is split into N groups of 2n individuals each. After
% mating and reproduction, a fraction d of each group exits the group for
% migration. The migrants are pooled together by gender and redistributed
% randomly on the spots left vacant by other migrants (males take male spots
% and females take female spots to prevent gender bias).
%
% Initial conditions are such that alleles A and a are at 50:50 proportions
% whereas copier females are present at frequency q0 (overall, may differ between groups).
% Simulations stop whenever genetic variation is lost in one of the genes.
%


%% INTRO %%

clear
close all

if ismac
    cd '/Users/Antoine/Documents/GitHub/project-conformity/code Sergey'
end

if 1==1
    rng('shuffle');
    seed=rng;
    save seed;
else
    load seed;
    rng(seed)
end


% Parameters
ntot = 1000;           % total number of individuals
N = 10;                % number of groups
n = ntot/(2*N);        % number of females/males per group

T = 1000;              % number of generations
T0 = 1;                % generation when conformity is introduced (minimum 1)
Runs = 20;             % number of runs

Q0 = .2; %[.1 .2 .5];        % initial frequency of copier females
D = 0:0.005:.10;             % dispersal rate
Beta = 1:.05:3;              % conformity strength parameter
R = .2; %[0 .05 .1 .5 1];    % recombination rate (0 never recombined; 1 always recombined; .5 baseline)

group_invasion = 0;    % conformity appears as a cluster (1) or not (0) 
vs_setprefs = 0;       % competition against set preferences for a given male phenotype (1) or random copying (0)

save_graphs = 0;
display_figures = 'off';
I = 1:n ;                        % for recombination


% Setting up end successes matrix: [ntot N T T0 vs_setprefs group_invasion q0 d beta r Successes/Runs]
nComb = length(R)*length(Q0)*length(Beta)*length(D);
MatComb = nan(nComb,11);

i=1;
for i_q0 = 1:length(Q0)
    q0 = Q0(i_q0);
for i_d = 1:length(D)
    d = D(i_d);
for i_beta = 1:length(Beta)
    beta = Beta(i_beta);
for i_r = 1:length(R)
    r = R(i_r);
    
    MatComb(i, 1:10) = [ntot N T T0 vs_setprefs group_invasion q0 d beta r];
    
    i = i+1;
    
end
end
end
end

if ismac
    %parpool(2);
    ;
elseif ispc
    parpool(8);
end

progressSimulations = 0;



%% SIMULATIONS VS RANDOM COPYING %%

if vs_setprefs == 0

tic

%parfor iComb = 1:nComb
for iComb = 1:nComb
    
    currentComb = MatComb(iComb,:);
    
    q0 = currentComb(7);
    d = currentComb(8);
    beta = currentComb(9);
    r = currentComb(10);
    
    P = nan(1,Runs);                  % vector of end-run proportions of the A allele 
	Q = nan(1,Runs);                  % vector of end-run proportions of the B allele (copier)
	Tau = nan(1,Runs);                % time to fixation of the male trait
	p = nan(T,N,Runs);                % frequency of matings by males of type 1 in each group
	q = nan(T,N,Runs);                % 
    
    Successes = 0;                    % number of 'successful' runs ('success' is defined below)
    
for run=1:Runs
    
    if T0 == 1
        q(1,:,run) = q0;                   % initial frequency of 'copier' gene (might differ between groups)
    else
        q(1,:,run) = 0;
    end

    Z0 = nan(2*n,2,N);                % create population genotypes: 2n individuals with 2 genes each, repeated for N groups
                                      % first (1:n) are females, last (n+1:2n) are males
    Z0(:,1,:) = 0;                    % initially only bearers of allele b (random copiers)
    Z0(:,2,:) = randi(2,2*n,1,N)-1;   % initial bearers of allele A
     
    p(1,:,run) = mean(Z0(:,2,:), 1);  % initial proportion of male phenotype A in each group
    Z = Z0;                           % initial genotypes in the population
	m1 = ones(1,N);                   % per-group count of allele A
	m2 = ones(1,N);
	f1 = ones(1,N);                   % per-group count of allele B
	f2 = ones(1,N);
	gametes = nan(T,4);               % ???
    
    t = 1;
	while sum(m1)*sum(m2)>0 && t<T %&& sum(f1)*sum(f2)>0    % to stop simulations if there is fixation in male or female phenotype
        
        if t == T0
            n_mod = q0*ntot/2;
            if group_invasion == 1
                k = floor(n_mod/n);
                for i = 1:k
                    Z(:,1,i) = 1;
                end
                n_rem = n_mod - k*n;
                Z(1:n_rem,1,N) = 1;
                Z(n+1:n+n_rem,1,N) = 1;
            else
                Z(:,1,:) = rand(2*n,1,N)<q0;
            end
        end

        X = Z(1:n,:,:);                   % genotypes of females
        Y = Z(n+1:2*n,:,:);               % genotypes of males
        
        % Reproduction within groups
        for g = 1:N

            F1 = find(X(:,1,g)==1);             % indexes of type B females (copier)
                f1(g) = length(F1);
            F2 = find(X(:,1,g)==0);             % indexes of type b females (random mating)
                f2(g) = n - f1(g);                % test: f2(g)==length(F2)

            M1 = find(Y(:,2,g)==1);             % indexes of type A males
                m1(g) = length(M1);
            M2 = find(Y(:,2,g)==0);             % indexes of type a males
                m2(g) = n - m1(g);                % test: m2(g)==length(M2)
    
            if m1(g) && m2(g)
                u = double(rand(f1(g),1) < p(t,g,run)^beta/(p(t,g,run)^beta+(1-p(t,g,run))^beta));   % females of type B choose among males of type 1 with probability p
                fathers_1 = u.*M1(randi(m1(g),f1(g),1))+(1-u).*M2(randi(m2(g),f1(g),1));    % indexes of fathers for B females
            else
                fathers_1 = randi(n,f1(g),1);        % indexes of fathers for B females
            end
            fathers_2 = randi(n,f2(g),1);            % indexes of fathers for b females

            mothers = [F1;F2];
            fathers = [fathers_1;fathers_2];
            p(t+1,g,run) = mean(Y(fathers,2,g));
            q(t+1,g,run) = mean(X(mothers,1,g));

            %G= Z(:,1)+2*Z(:,2);
            %h=hist(G,edges);
            %gametes(t,:)=h;

            % Four types of offspring genotypes (vertically: both from mother; 1 from mother and 2 from father; etc.)
            First =  [X(mothers,1,g)'; X(mothers,1,g)'; Y(fathers,1,g)'; Y(fathers,1,g)'];
            Second = [X(mothers,2,g)'; Y(fathers,2,g)'; X(mothers,2,g)'; Y(fathers,2,g)'];

            % the type from recombination probabilities
            ind1 = randp([(1-r)/2, r/2, r/2, (1-r)/2], n, 1);
            Ind1 = 4*(I'-1) + ind1;
            ind2 = randp([(1-r)/2, r/2, r/2, (1-r)/2], n, 1);
            Ind2 = 4*(I'-1) + ind2;

            Z(1:n,:,g) = [First(Ind1) Second(Ind1)];          % new generation (child 1)
            Z(n+1:2*n,:,g) = [First(Ind2) Second(Ind2)];      % new generation (child 2)
            Z(:,:,g) = Z(randperm(2*n),:,g);                  % assign sex randomly
           
        end
        t = t+1;
        
        % Dispersal: shuffling of some individuals across groups; males replace males and females replace females
        Dispersal = rand(2*n,N)<d;
        
        [DFrow,DFcol] = find(Dispersal(1:n,:)==1);       % indexes of female migrants
        DF = [DFrow,DFcol];
        df = size(DF,1);
        shuffleF = DF;
        idx = randperm(df);
        shuffleF(idx,:) = DF(:,:);               % shuffling of migrant females over all vacant spots
        
        [DMrow,DMcol] = find(Dispersal(n+1:2*n,:)==1);   % indexes of male migrants
        DM = [DMrow,DMcol];
        dm = size(DM,1);
        shuffleM = DM;
        idx = randperm(dm);
        shuffleM(idx,:) = DM(:,:);               % shuffling of migrant males over all vacant spots
        
        Zpost = Z;
        
        for i=1:df
           Zpost(shuffleF(i,1),1,shuffleF(i,2)) = Z(DF(i,1),1,DF(i,2)) ;        % female shuffleF(i) replaces female i (first gene)
           Zpost(shuffleF(i,1),2,shuffleF(i,2)) = Z(DF(i,1),2,DF(i,2)) ;        % (second gene)
        end
        for i=1:dm
           Zpost(n+shuffleM(i,1),1,shuffleM(i,2)) = Z(n+DM(i,1),1,DM(i,2)) ;    % male shuffleM(i) replaces male i (first gene)
           Zpost(n+shuffleM(i,1),2,shuffleM(i,2)) = Z(n+DM(i,1),2,DM(i,2)) ;    % (second gene)
        end
        
        Z=Zpost;
        
	end
    P(run)=mean(p(t,:,run));
    Q(run)=mean(q(t,:,run));
    Tau(run)=t;
    if P(run)*(1-P(run))>0 && Q(run) > q0      % just one way to define success: male heterogeneity and more conformists than at first
        Successes = Successes+1;
    end
end

maxTau=max(Tau(:));

MatComb(iComb, :) = [ntot N T T0 vs_setprefs group_invasion q0 d beta r Successes/Runs];
%if sum(ismember(Success_rate(:, 1:10), [ntot N T T0 vs_setprefs group_invasion q0 d beta r], 'rows'))==0     % only add data point if previously inexistant
%    Success_rate = [Success_rate; ntot N T T0 vs_setprefs group_invasion q0 d beta r Successes/Runs];
%    if ispc
%        save('Success_rate.mat', 'Success_rate')   % run only on the server
%    end
%end

figure('visible',display_figures);
subplot(2,1,1)
datap=squeeze(mean(p(1:maxTau,:,:),2));
plot(1:maxTau,datap);
%ylabel('p')
ylim([0 1])
%title([num2str(N), ' groups of ', num2str(2*n), ' individuals each, q_0=', num2str(q0),', d=',num2str(d),', \beta=',num2str(beta),', r=',num2str(r)]);

subplot(2,1,2)
dataq=squeeze(mean(q(1:maxTau,:,:),2));
plot(1:maxTau,dataq);
%ylabel('q')
ylim([0 1])
%xlabel('time')

png_file = '';
if save_graphs == 1
    if ismac
        png_file = sprintf('/Users/Antoine/Documents/GitHub/project-conformity/code Sergey/Figs/randmatch q%2.1f d%2.2f b%3.2f  r%2.2f.png',q0,d,beta,r);
    elseif ispc
        png_file = sprintf('Z:/Conformity/MATLAB/Figs/randmatch q%2.1f d%2.2f b%3.2f r%2.2f.png',q0,d,beta,r);
    end
	print(png_file,'-dpng','-r150');
%   print('-depsc2',eps_file);
%   print('-dpdf',pdf_file,'-bestfit');
end

progressSimulations = progressSimulations + 1

disp(['Progress: ', num2str(100*progressSimulations/nComb, '%2.2f') ,'% completed'])

end

toc

end



%% SIMULATIONS VS SET PREFERENCES %%

if vs_setprefs == 1

tic

parfor iComb = 1:nComb
    
    currentComb = MatComb(iComb,:);
    
    q0 = currentComb(7);
    d = currentComb(8);
    beta = currentComb(9);
    r = currentComb(10);
    
    P = nan(1,Runs);                  % vector of end-run proportions of the A allele 
	Q = nan(1,Runs);                  % vector of end-run proportions of the C allele (copier)
	Tau = nan(1,Runs);                % time to fixation of the male trait
	p = nan(T,N,Runs);                % frequency of matings by males of type A in each group
	q = nan(T,N,Runs);                % 
    
    Successes = 0;                    % number of 'successful' runs ('success' is defined below)
    
for run=1:Runs
    
    if T0 == 1
        q(1,:,run) = q0;                   % initial frequency of 'copier' gene (might differ between groups)
    else
        q(1,:,run) = 0;
    end

    Z0 = nan(2*n,2,N);                % create population genotypes: 2n individuals with 2 genes each, repeated for N groups
                                      % first (1:n) are females, last (n+1:2n) are males
    Z0(:,1,:) = randi(2,2*n,1,N)-1;   % initial bearers of allele ca (0), or cA (1) [C introduced later as 2]
    Z0(:,2,:) = randi(2,2*n,1,N)-1;   % initial bearers of allele A
     
    p(1,:,run) = mean(Z0(:,2,:), 1);  % initial proportion of male phenotype A in each group
    Z = Z0;                           % initial genotypes in the population
	m1 = ones(1,N);                   % per-group count of allele A
	m2 = ones(1,N);
	f1 = ones(1,N);                   % per-group count of allele ca
	f2 = ones(1,N);                   % cA
    f3 = ones(1,N);                   % C
	gametes = nan(T,4);               % ???
    
    t = 1;
	while sum(m1)*sum(m2)>0 && sum(f1)*sum(f2+f3)>0 && t<T    % to stop simulations if there is fixation in male phenotype
        
        if t == T0
            n_mod = q0*ntot/2;
            if group_invasion == 1
                k = floor(n_mod/n);
                for i = 1:k
                    Z(:,1,i) = 2;
                end
                n_rem = n_mod - k*n;
                Z(1:n_rem,1,N) = 2;
                Z(n+1:n+n_rem,1,N) = 2;
            else
                Mutants = rand(2*n,1,N)<q0;
                for g = 1:N
                    for i = 1:2*n
                        if Mutants(i,1,g) == 1
                            Z(i,1,g) = 2;
                        end
                    end
                end
            end
        end

        X = Z(1:n,:,:);                   % genotypes of females
        Y = Z(n+1:2*n,:,:);               % genotypes of males
        
        for g = 1:N

            F1 = find(X(:,1,g)==0);             % indexes of type ca females
                f1(g) = length(F1);
            F2 = find(X(:,1,g)==1);             % indexes of type cA females
                f2(g) = length(F2);
            F3 = find(X(:,1,g)==2);             % indexes of type C females
                f3(g) = n - f1(g) - f2(g);      % test: f3(g)==length(F3)

            M1 = find(Y(:,2,g)==0);             % indexes of type a males
                m1(g) = length(M1);
            M2 = find(Y(:,2,g)==1);             % indexes of type A males
                m2(g) = n - m1(g);              % test: m2(g)==length(M2)
    
            if m1(g)&&m2(g)
                fathers_1 = M1(randi(m1(g),f1(g),1));            % indexes of fathers for ca females
                fathers_2 = M2(randi(m2(g),f2(g),1));            % indexes of fathers for cA females
                u = double(rand(f3(g),1) < p(t,g,run)^beta/(p(t,g,run)^beta+(1-p(t,g,run))^beta));   % females of type C choose among males of type 1 with probability p
                fathers_3 = u.*M1(randi(m1(g),f3(g),1))+(1-u).*M2(randi(m2(g),f3(g),1));             % indexes of fathers for C females
            else
                fathers_1 = randi(n,f1(g),1);
                fathers_2 = randi(n,f2(g),1);
                fathers_3 = randi(n,f3(g),1);
            end

            mothers = [F1;F2;F3];
            fathers = [fathers_1;fathers_2;fathers_3];
            p(t+1,g,run) = mean(Y(fathers,2,g));
            q(t+1,g,run) = sum(X(mothers,1,g)==2, 'all')/n;

            %G= Z(:,1)+2*Z(:,2);            % ???
            %h=hist(G,edges);
            %gametes(t,:)=h;

            % Four types of offspring genotypes (vertically: both from mother; 1 from mother and 2 from father; etc.)
            First =  [X(mothers,1,g)'; X(mothers,1,g)'; Y(fathers,1,g)'; Y(fathers,1,g)'];
            Second = [X(mothers,2,g)'; Y(fathers,2,g)'; X(mothers,2,g)'; Y(fathers,2,g)'];

            % the type from recombination probabilities
            ind1 = randp([(1-r)/2, r/2, r/2, (1-r)/2], n,1);
            Ind1 = 4*(I'-1) + ind1;
            ind2 = randp([(1-r)/2, r/2, r/2, (1-r)/2], n,1);
            Ind2 = 4*(I'-1) + ind2;

            Z(1:n,:,g) = [First(Ind1) Second(Ind1)];          % new generation (child 1)
            Z(n+1:2*n,:,g) = [First(Ind2) Second(Ind2)];      % new generation (child 2)
            Z(:,:,g) = Z(randperm(2*n),:,g);                  % assign sex randomly
           
        end
        t = t+1;
        
        % Dispersal: shuffling of some individuals across groups; males replace males and females replace females
        Dispersal = rand(2*n,N)<d;
        
        [DFrow,DFcol] = find(Dispersal(1:n,:)==1);       % indexes of female migrants
        DF = [DFrow,DFcol];
        df = size(DF,1);
        shuffleF = DF;
        idx = randperm(df);
        shuffleF(idx,:) = DF(:,:);                       % shuffling of migrant females over all vacant spots
        
        [DMrow,DMcol] = find(Dispersal(n+1:2*n,:)==1);   % indexes of female migrants
        DM = [DMrow,DMcol];
        dm = size(DM,1);
        shuffleM = DM;
        idx = randperm(dm);
        shuffleM(idx,:) = DM(:,:);                       % shuffling of migrant males over all vacant spots
        
        Zpost = Z;
        
        for i=1:df
           Zpost(shuffleF(i,1),1,shuffleF(i,2)) = Z(DF(i,1),1,DF(i,2)) ;        % female shuffleF(i) replaces female i (first gene)
           Zpost(shuffleF(i,1),2,shuffleF(i,2)) = Z(DF(i,1),2,DF(i,2)) ;        % (second gene)
        end
        for i=1:dm
           Zpost(n+shuffleM(i,1),1,shuffleM(i,2)) = Z(n+DM(i,1),1,DM(i,2)) ;    % male shuffleM(i) replaces male i (first gene)
           Zpost(n+shuffleM(i,1),2,shuffleM(i,2)) = Z(n+DM(i,1),2,DM(i,2)) ;    % (second gene)
        end
        
        Z=Zpost;
        
	end
    P(run)=mean(p(t,:,run));
    Q(run)=mean(q(t,:,run));
    Tau(run)=t;
    if P(run)*(1-P(run))>0 && Q(run) > q0      % just one way to define success: male heterogeneity and more conformists than at first
        Successes = Successes+1;
    end
end

maxTau=max(Tau(:));

MatComb(iComb, :) = [ntot N T T0 vs_setprefs group_invasion q0 d beta r Successes/Runs];

%if sum(ismember(Success_rate(:, 1:10), [ntot N T T0 vs_setprefs group_invasion q0 d beta r], 'rows'))==0     % only add data point if previously inexistant
%    Success_rate = [Success_rate; ntot N T T0 vs_setprefs group_invasion q0 d beta r Successes/Runs];
%    if ispc
%        save('Success_rate.mat', 'Success_rate')   % run only on the server
%    end
%end

figure
subplot(2,1,1)
datap=squeeze(mean(p(1:maxTau,:,:),2));
plot(1:maxTau,datap);
ylabel('p')
ylim([0 1])
title([num2str(N), ' groups of ', num2str(2*n), ' individuals each, q_0=', num2str(q0),', d=',num2str(d),', \beta=',num2str(beta),', r=',num2str(r)]);

subplot(2,1,2)
dataq=squeeze(mean(q(1:maxTau,:,:),2));
plot(1:maxTau,dataq);
ylabel('q')
ylim([0 1])
xlabel('time')

png_file = '';
if save_graphs == 1
    if ismac
        png_file = sprintf('/Users/Antoine/Documents/GitHub/project-conformity/code Sergey/Figs/setprefs q%2.1f d%2.2f b%3.2f  r%2.2f.png',q0,d,beta,r);
    elseif ispc
        png_file = sprintf('Z:/Conformity/MATLAB/Figs/setprefs q%2.1f d%2.2f b%3.2f r%2.2f.png',q0,d,beta,r);
    end
	print(png_file,'-dpng');
%   print('-depsc2',eps_file);
%   print('-dpdf',pdf_file,'-bestfit');
end

progressSimulations = progressSimulations + 1

disp(['Progress: ', num2str(100*progressSimulations/nComb, '%2.2f') ,'% completed'])

end

toc

end


%delete(gcp('nocreate'));


%% SAVE SIMULATION DATA %%

load('ConformitySimulations.mat', 'ConformitySimulations')
ConformitySimulations = [ConformitySimulations; MatComb];
save('ConformitySimulations.mat', 'ConformitySimulations')





