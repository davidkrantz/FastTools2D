% This is a test script for the doubly periodic Stokes single-layer
% potential

close all
clearvars
clc

initewald

%% Set up data
Nsrc = 10000;
Ntar = 10000;

Lx = 2*pi;
Ly = 2*pi;

% Two components of the density function
f1 = 10*rand(Nsrc,1);
f2 = 10*rand(Nsrc,1);

% Source and target loccations
xsrc = Lx*rand(Nsrc,1);
ysrc = Ly*rand(Nsrc,1);

xtar = Lx*rand(Ntar,1);
ytar = Ly*rand(Ntar,1);

%% Check that timings scale as O(N Log N)

Nsrc = 10.^[3, 4, 5, 6];
tol = 10.^[-4, -8, -12, -16];

% NB: Choosing too low a tolerance can lead to problems...

Lx = 1;
Ly = 2;

times = zeros(length(tol), length(Nsrc));

loglog(Nsrc, Nsrc.*log(Nsrc));
xlabel('Number of points');
ylabel('time (s)');
drawnow;

hold on

for i = 1:length(tol)
    for j = 1:length(Nsrc)
        % Two components of the density function
        f1 = 10*rand(Nsrc(j),1);
        f2 = 10*rand(Nsrc(j),1);
        
        % Source and target loccations
        xsrc = Lx*rand(Nsrc(j),1);
        ysrc = Ly*rand(Nsrc(j),1);
        
        xtar = Lx*rand(Nsrc(j),1);
        ytar = Ly*rand(Nsrc(j),1);
        
        tic
        [u1, u2] = StokesSLP_ewald_2p(xsrc, ysrc, xtar, ytar, f1, f2, Lx, Ly, 'tol', tol(i), 'verbose', 1);
        
        times(i,j) = toc;
                
    end
    
    loglog(Nsrc, times(i,:), '-o');
    drawnow
    
end

%% Check that timings scale as O(N Log N)

Nsrc = 10.^[3, 4, 5, 6];
tol = 10.^[-4, -8, -12, -16];

Lx = 1;
Ly = 2;

times = zeros(length(tol), length(Nsrc));

loglog(Nsrc, Nsrc.*log(Nsrc));
xlabel('Number of points');
ylabel('time (s)');
drawnow;

hold on

for i = 1:length(tol)
    for j = 1:length(Nsrc)
        %Two components of the density function
        f1 = 10*rand(Nsrc(j),1);
        f2 = 10*rand(Nsrc(j),1);
        
        % normal vector at source points
        n1 = rand(Nsrc(j),1);
        n2 = sqrt(1 - n1.^2);
        
        % Source and target loccations
        xsrc = Lx*rand(Nsrc(j),1);
        ysrc = Ly*rand(Nsrc(j),1);
        
        xtar = Lx*rand(Nsrc(j),1);
        ytar = Ly*rand(Nsrc(j),1);
        
        tic
        [u1, u2] = StokesDLP_ewald_2p(xsrc, ysrc, xtar, ytar, n1, n2, ...
            f1, f2, Lx, Ly, 'tol', tol(i), 'verbose', 1);
        
        times(i,j) = toc;
              
    end
    
    loglog(Nsrc, times(i,:), '-o');
    drawnow
    
end
        

