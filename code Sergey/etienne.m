% Apr.18-19, 2020: checking Etienne's idea on sex and conformity
% S. Gavrilets
% There are 2 diallelic genes one of which with allels A,a is expressed in males 
% only and another with alleles B, b in females only; the recombinatio rate is r.
% The frequencies of alleels A and B are pa nd 1, respectively. Male trait
% is neutral. Females b mate randomly. Copier females B prefer males that had high mating
% success in the previous generation. I am using an S-shaped "conformity
% function" in which probablity of mating with males of a paticular type is
% p^beta/(p^beta+(1-p)^beta), where p is their mating success in the previous 
% generation and beta is the steepnes parameter. The population size is
% 2000 individuals. Initial conditions are such that alleles B and b are at 50:50 proportions
% whereis copier females are present at frequency q0. 100 independent runs
% for each parameter combination (as shown in the graphs). Simulations stop 
% whenever genetic varition is lost in one of the genes.
%
% The results show that while copier females can rapidly drive the most common 
% male trait to fixation, there is no detectable selection on female alleles.
% I think this result is actually compatible with Kirkpatrick (1982).
%
% Three possibilities: your scenario doesn't work, I misunderstood your
% scenario and implemented it in a wrong way, my code has some bugs.


clear all
close all

if 1
    rng('shuffle');
    seed=rng;
    save seed;
else
    load seed;
    rng(seed)
end
ones(10)*ones(10);  


% Parameters
n=1000;             % number of females and of males
T=500;             % nunber of generations
Runs=100;           % number of runs

save_graphs=0;

R=0.05; %[0 0.01 0.05 0.1 0.5];    % recombination rate
Q0= 0.5; %[0.1, 0.5];             % initial frequency of copier females
Beta=1.5; % [1.25, 1.50, 2.0];    % conformity strength parameter


% Useful
I=1:2*n;                    % for recombination
edges=-0.5:1:2.5;           % for histogram

tic

for q0=Q0
    for beta=Beta
        for i=1:length(R)
            r=R(i);
            [r q0 beta]
            
            P=nan(length(R),Runs);
            Q=nan(length(R),Runs);
            Tau=nan(length(R),Runs);        % time to fixation of the male trait
            p=nan(T,Runs);                  % frequency of matings by males of type 1
            q=nan(T,Runs);  
            
            for run=1:Runs
                p(1,run)=0.5;
                q(1,run)=q0;

                Z0=zeros(2*n,1);
                Z0(:,1)=rand(2*n,1)<q0;
                Z0(:,2)=randi(2,2*n,1)-1;                
                Z=Z0;
                m1=1;
                m2=1;
                f1=1;
                f2=1;
                gametes=nan(T,4);


        %         for t=1:T
                t=1;
                while m1*m2>0 && f1*f2>0 && t<T     % to stop simulations if there is fixation in ne gene or another

                    X=repmat(Z(1:n,:),2,1);         % genotypes of mothers of 2n offspring
                    Y=Z(n+1:2*n,:);                 % genotype of males

                    F1=find(X(:,1)==1);   % indexes of type 1 females (copier)
                    f1=length(F1);
                    F2=find(X(:,1)==0);   % indexes of type 0 females (random mating)
                    f2=2*n-f1;

                    Y1=find(Y(1:n,2)==1);   % indexes of type 1 males
                    m1=length(Y1);
                    Y2=find(Y(1:n,2)==0);   % indexes of type 0 males
                    m2=n-m1;

        %             u=double(rand(f1,1)<p(t)+e*p(t)*(1-p(t))*(2*p(t)-1));      % females of type 1 choose among males of type 1 with probability p
                    u=double(rand(f1,1)<p(t,run)^beta/(p(t,run)^beta+(1-p(t,run))^beta));
                    if m1&&m2
                        fathers_1=u.*Y1(randi(m1,f1,1))+(1-u).*Y2(randi(m2,f1,1));  
                    else
                        fathers_1=randi(n,f1,1);        % indexes of fathers
                    end
                    fathers_2=randi(n,f2,1);            % indexes of fathers

                    mothers=[F1;F2];
                    fathers=[fathers_1;fathers_2];
                    p(t+1,run)=mean(Y(fathers,2));
                    q(t+1,run)=mean(X(mothers,1));

                    G= Z(:,1)+2*Z(:,2);
                    h=hist(G,edges);
                    gametes(t,:)=h;

                    % For types of offspring genotypes (vertically)
                    First= [X(mothers,1)'; X(mothers,1)'; Y(fathers,1)'; Y(fathers,1)'];
                    Second=[X(mothers,2)'; Y(fathers,2)'; X(mothers,2)'; Y(fathers,2)'];

                    % the type from recombination probabilities
                    ind=randp([(1-r)/2, r/2, r/2, (1-r)/2], 2*n,1);
                    Ind=4*(I'-1)+ind;

                    Z=[First(Ind) Second(Ind)];     % new generation
                    Z=Z(randperm(2*n),:);           % assign sex randomly
                    t=t+1;
                end
                    P(i,run)=p(t,run);
                    Q(i,run)=q(t,run);
                    Tau(i,run)=t;
            end
            maxTau=max(Tau(i,:));
            
            figure
            subplot(2,1,1)
            plot(1:maxTau,p(1:maxTau,:));
            ylabel('p')
            ylim([0 1])
            title(['r=', num2str(r),', q_0=',num2str(q0),', \beta=',num2str(beta)]);

            subplot(2,1,2)
            plot(1:maxTau,q(1:maxTau,:));
            ylabel('q')
            ylim([0 1])
            xlabel('time')

            if save_graphs 
                png_file=sprintf('Figs/r%2.2fq%2.1f.b%3.2f.png',r,q0,beta);
                print('-dpng',png_file);
        %         print('-depsc2',eps_file);
        %         print('-dpdf',pdf_file,'-bestfit');
            end

        end
    end
end
toc

[mean(P,2) std(P')']
[mean(Q,2) std(Q')']
[mean(Tau,2) std(Tau')']
