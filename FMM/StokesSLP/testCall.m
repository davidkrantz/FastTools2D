%addpath ../src
%Test the FMM for stokes SLP.
N = 2^13;
theta = (0:N-1)'*2*pi/N;

f1 = [ones(N,1);ones(N,1)]*2*pi/64;
f2 = [ones(N,1);ones(N,1)]*2*pi/64;
%Two components of the density function

xs = [cos(theta)-1.02 ;cos(theta) + 1.02];
ys = [2/3*sin(theta) ;2/3*sin(theta)];
%source/target locations

%f1 = ones(N,1);
%xs = cos(theta)*2*pi/N;
%ys = sin(theta)*2*pi/N;

tic
[u1,v1] = stokesSLPfmm(f1(:),f2(:),xs(:),ys(:));
toc


u2 = zeros(numel(f1),1);
v2 = zeros(numel(f1),1);
tic
for k = 1:numel(f1)
  ind = [(1:k-1) (k+1:numel(f1))];
  rx = xs(k) - xs(ind);
  ry = ys(k) - ys(ind);
  rho2 = rx.^2 + ry.^2;
  rdots = rx.*f1(ind) + ry.*f2(ind);
  u2(k) = sum(-0.5*log(rho2).*f1(ind) + rdots./rho2.*rx);
  v2(k) = sum(-0.5*log(rho2).*f2(ind) + rdots./rho2.*ry);
end
u2 = u2/4/pi;
v2 = v2/4/pi;
toc

