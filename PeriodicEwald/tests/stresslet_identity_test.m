% Check stresslet identity

fprintf("*********************************************************\n");
fprintf("Checking stresslet identity\n")
fprintf("*********************************************************\n\n");

initewald

r = 1;
r_in = 0.5*r;
r_out = 1.5*r;

Lx = 4*r;
Ly = 4*r;

%% set up source and targets
Nsrc = 256;
Ntar = 16;

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

%% Call spectral Ewald and check average values
[u1, u2] = StokesDLP_ewald_2p(xsrc, ysrc, xtar, ytar, n1, n2, f1, f2, Lx, Ly);

fprintf("f1 = 1, f2 = 0\n");
fprintf("error in mean interior solution: (%3.3e, %3.3e)\n",...
    mean(u1(1:end/2) + 1), mean(u2(1:end/2)))
fprintf("error in mean exterior solution (%3.3e, %3.3e)\n",...
    mean(u1(end/2+1:end)), mean(u2(end/2+1:end))); 

% y component = 1
f1 = zeros(Nsrc,1);
f2 = hsrc*ones(Nsrc,1);

[u1, u2] = StokesDLP_ewald_2p(xsrc, ysrc, xtar, ytar, n1, n2, f1, f2, Lx, Ly);

fprintf("error in mean interior solution: (%3.3e, %3.3e)\n",...
    mean(u1(1:end/2)), mean(u2(1:end/2) + 1))
fprintf("error in mean exterior solution (%3.3e, %3.3e)\n",...
    mean(u1(end/2+1:end)), mean(u2(end/2+1:end))); 




