% Check stresslet identity. To do this we set up a circle of radius r, and
% prescribe a density f = e_i (i = 1,2). Then if u(x) is given by the
% double-layer potential, the following identity holds, even for periodic
% sums:
% u(x) = e_i, x inside circle, 
%         0, x outside circle,
%        e_i/2, x on the circle.
% We will consider only the first two cases, and look at the average u at 
% a circle of points inside and another one outside the original circle. 

clearvars
close all
addpath('../src/StokesDLP');

fprintf("*********************************************************\n");
fprintf("Checking stresslet identity\n")
fprintf("*********************************************************\n\n");

r = 1;
r_in = 0.5*r;
r_out = 1.5*r;

%% set up source and targets
Nsrc = 512;
Ntar = 16;

% use trapezoid rule
hsrc = 2*pi/Nsrc;
theta = (hsrc : hsrc : 2*pi)';

xsrc = r*cos(theta);
ysrc = r*sin(theta);

% x component = 1
f1 = hsrc*ones(Nsrc,1);
f2 = zeros(Nsrc,1);

n1 = cos(theta);
n2 = sin(theta);

htar = 2*pi/Ntar;
theta = (htar : htar : 2*pi)';

xtar = zeros(2*Ntar, 1);
ytar = zeros(2*Ntar, 1);

xtar(1:end/2) = r_in*cos(theta);
xtar(end/2+1:end) = r_out*cos(theta);
ytar(1:end/2) = r_in*sin(theta);
ytar(end/2+1:end) = r_out*sin(theta);

%% Call FMM and check average values

% NOTE double-layer FMM assumes identical source and target points, so we
% pass in the target points with zero weights

xsrc = [xsrc; xtar];
ysrc = [ysrc; ytar];
f1 = [f1; zeros(2*Ntar,1)];
f2 = [f2; zeros(2*Ntar,1)];

[u1, u2] = stokesDLPfmm(f1(:),f2(:),xsrc(:),ysrc(:),n1(:),n2(:));

% FMM not scaled by -1
u1 = -u1;
u2 = -u2;
fprintf("f1 = 1, f2 = 0\n");
fprintf("error in mean interior solution: %3.3e\n",...
    norm(u1(Nsrc+1:Nsrc+Ntar) + 1i*u2(Nsrc+1:Nsrc+Ntar) - 1));
fprintf("error in mean exterior solution %3.3e\n",...
    norm(u1(Nsrc+Ntar+1:end) + 1i*u2(Nsrc+Ntar+1:end))); 

% y component = 1
f1 = zeros(Nsrc,1);
f2 = hsrc*ones(Nsrc,1);
f1 = [f1; zeros(2*Ntar,1)];
f2 = [f2; zeros(2*Ntar,1)];

[u1, u2] = stokesDLPfmm(f1(:),f2(:),xsrc(:),ysrc(:),n1(:),n2(:));

% FMM not scaled by -1
u1 = -u1;
u2 = -u2;

fprintf("f1 = 0, f2 = 1\n");
fprintf("error in mean interior solution: %3.3e\n",...
    norm(u1(Nsrc+1:Nsrc+Ntar) + 1i*u2(Nsrc+1:Nsrc+Ntar) - 1i));
fprintf("error in mean exterior solution %3.3e\n",...
    norm(u1(Nsrc+Ntar+1:end) + 1i*u2(Nsrc+Ntar+1:end))); 




