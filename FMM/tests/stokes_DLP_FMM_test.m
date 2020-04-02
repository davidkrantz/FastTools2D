% Test the FMM for stokes DLP
clearvars
close all
addpath('../src/StokesDLP');

%% Set up data
Nsrc = 10000;

%Two components of the density function
f1 = 10*rand(Nsrc,1);
f2 = 10*rand(Nsrc,1);

n1 = rand(Nsrc,1);
n2 = sqrt(1 - n1.^2);

% Source and target loccations
xsrc = 2*rand(Nsrc,1) - 1;
ysrc = 2*rand(Nsrc,1) - 1;


%% Compute direct sum
udirect = zeros(numel(f1),1);
vdirect = zeros(numel(f1),1);

tic
for k = 1:numel(f1)
    
    % skip self interaction term
    ind = [(1:k-1) (k+1:numel(f1))];

    rx = xsrc(k) - xsrc(ind);
    ry = ysrc(k) - ysrc(ind);
    rho4 = (rx.^2 + ry.^2).^2;
    
    rdotf = rx.*f1(ind) + ry.*f2(ind);
    rdotn = rx.*n1(ind) + ry.*n2(ind);
    
    udirect(k) = 4*sum(rdotn.*rdotf./rho4.*rx);
    vdirect(k) = 4*sum(rdotn.*rdotf./rho4.*ry);
end

udirect = udirect/4/pi;
vdirect = vdirect/4/pi;
time_direct = toc;

disp(['Direct sum for ', num2str(Nsrc) ' sources took ',...
        num2str(time_direct), ' seconds.']);

%% Compute FMM for a range of tolerances

tic
[ufmm,vfmm] = stokesDLPfmm(f1(:),f2(:),xsrc(:),ysrc(:),n1(:),n2(:));

t = toc;
e = max(abs(ufmm + 1i*vfmm - (udirect + 1i*vdirect)));

disp(['FMM took ', num2str(t),...
    ' seconds. Maximum error: ', num2str(e)]);


%% Check time scaling
Npts = [1e3, 1e4, 1e5, 1e6];
times_fmm = zeros(size(Nsrc));

for i = 1:length(Npts)
    %Two components of the density function
    f1 = 10*rand(Npts(i),1);
    f2 = 10*rand(Npts(i),1);
    
    n1 = rand(Npts(i),1);
    n2 = sqrt(1 - n1.^2);

    % Source and target loccations
    xsrc = 2*rand(Npts(i),1) - 1;
    ysrc = 2*rand(Npts(i),1) - 1;
    
    tic
    [ufmm,vfmm] = stokesDLPfmm(f1(:),f2(:),xsrc(:),ysrc(:),n1(:),n2(:));
    
    times_fmm(i) = toc;
    
end

figure()
loglog(Npts, times_fmm, '-o');
hold on
loglog(Npts, times_fmm(1)*Npts/Npts(1));

xlabel('number of points');
ylabel('time');
legend({'FMM', '$\mathcal{O}(N)$'}, 'interpreter', 'latex');




