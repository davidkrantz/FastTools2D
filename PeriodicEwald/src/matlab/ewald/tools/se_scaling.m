function Hhat = se_scaling(ghat,ewaldparam,alpha)
% Compute scaling of k-space K0 term for the spectral Ewald method.
% Periodic case, Gaussian window function

kvecx = -ewaldparam.kinfx + (0:ewaldparam.Mx-1)/ewaldparam.Mx*2*ewaldparam.kinfx;
kvecy = -ewaldparam.kinfy + (0:ewaldparam.My-1)/ewaldparam.My*2*ewaldparam.kinfy;
[KJ,KL] = meshgrid(kvecx, kvecy);
K2 = (KJ.^2 + KL.^2);
Kmod = sqrt(K2);

% Scaling factor
AF = 2*pi./(alpha^2+K2).*exp(-alpha^2/(4*ewaldparam.xi^2));

% Window  function
B = exp(-(1-ewaldparam.eta)*K2/(4*ewaldparam.xi^2));


Hhat = B.*AF.*ghat;

end