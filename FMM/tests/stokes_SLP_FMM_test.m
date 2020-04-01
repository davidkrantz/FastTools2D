% Test the FMM for stokes SLP
clearvars
close all
addpath('../src/StokesSLP');

%% Set up data
Nsrc = 10000;
Ntar = 10000;

%Two components of the density function
f1 = 10*rand(Nsrc,1);
f2 = 10*rand(Nsrc,1);

% Source and target loccations
xsrc = 2*rand(Nsrc,1) - 1;
ysrc = 2*rand(Nsrc,1) - 1;

xtar = 2*rand(Ntar,1) - 1;
ytar = 2*rand(Ntar,1) - 1;

srcEqualsTar = 0;

if srcEqualsTar
    xtar = xsrc;
    ytar = ysrc;
end


%% Compute direct sum
udirect = zeros(numel(f1),1);
vdirect = zeros(numel(f1),1);

tic
for k = 1:numel(f1)
    
    if srcEqualsTar
        ind = [(1:k-1) (k+1:numel(f1))];
    else
        ind = 1:Nsrc;
    end
    rx = xtar(k) - xsrc(ind);
    ry = ytar(k) - ysrc(ind);
    rho2 = rx.^2 + ry.^2;
    rdotf = rx.*f1(ind) + ry.*f2(ind);
    udirect(k) = sum(-0.5*log(rho2).*f1(ind) + rdotf./rho2.*rx);
    vdirect(k) = sum(-0.5*log(rho2).*f2(ind) + rdotf./rho2.*ry);
end

udirect = udirect/4/pi;
vdirect = vdirect/4/pi;
time_direct = toc;

disp(['Direct sum for ', num2str(Nsrc) ' sources and ', num2str(length(xtar))...
    ' targets took ', num2str(time_direct), ' seconds.']);

%% Compute FMM for a range of tolerances
iprec = -2:5;

for i = 1:length(iprec)
    tic
    [ufmm,vfmm] = stokesSLPfmm(f1(:),f2(:),xsrc(:),ysrc(:),xtar(:),ytar(:),...
        srcEqualsTar,iprec(i));
    
    t = toc;
    e = max(abs(ufmm + 1i*vfmm - (udirect+1i*vdirect)));
    
    disp(['FMM with iprec ', num2str(iprec(i)), ' took ', num2str(t),...
        ' seconds. Maximum error: ', num2str(e)]);
    
end

%% Check time scaling
Npts = [1e3, 1e4, 1e5, 1e6];
times_fmm = zeros(size(Nsrc));

for i = 1:length(Npts)
    %Two components of the density function
    f1 = 10*rand(Npts(i),1);
    f2 = 10*rand(Npts(i),1);
    
    % Source and target loccations
    xsrc = 2*rand(Npts(i),1) - 1;
    ysrc = 2*rand(Npts(i),1) - 1;
    
    xtar = 2*rand(Npts(i),1) - 1;
    ytar = 2*rand(Npts(i),1) - 1;
    
    tic
    [ufmm,vfmm] = stokesSLPfmm(f1(:),f2(:),xsrc(:),ysrc(:),xtar(:),ytar(:),...
        0,5);
    
    times_fmm(i) = toc;
    
end

figure()
loglog(Npts, times_fmm, '-o');
hold on
loglog(Npts, times_fmm(1)*Npts/Npts(1));

xlabel('number of points');
ylabel('time');
legend({'FMM', '$\mathcal{O}(N)$'}, 'interpreter', 'latex');




