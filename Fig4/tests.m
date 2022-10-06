% tests
n = 200 ;
N = 10 ;
q0 = .5 ;
r = 0 ;
threshold = q0*2*n*N ;
Runs = 100000 ;
Z0 = nan(2*n,2,N) ;
count = 0 ;

Beta = 1:1:200;
step = .005;

plotbeta = nan(1,1+1/step);

for beta = Beta
for p = 0:step:1
    i = int8(1+p/step);
    plotbeta(beta, i) = p^beta/(p^beta+(1-p)^beta);
end
end

figure
plot(0:step:1,plotbeta);


