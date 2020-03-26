function domain = domSetup(param)


% SOURCES -> Place random sources
rng('default')
domain.psrc = zeros(param.Nsrc,2);
domain.psrc(:,1) = param.Lx*rand(param.Nsrc,1) - param.Lx/2;
domain.psrc(:,2) = param.Ly*rand(param.Nsrc,1) - param.Ly/2;

% Randomise force at source
domain.fsrc = zeros(param.Nsrc,param.nbrDim);
domain.Q = 0;
for j=1:param.nbrDim
    domain.fsrc(:,j) = rand(param.Nsrc,1);
    domain.fsrc(:,j) = domain.fsrc(:,j)/sum(domain.fsrc(:,j));
    domain.Q = domain.Q + sum(domain.fsrc(:,j).^2);
end
% domain.fsrc = domain.fsrc - sum(domain.fsrc)/param.Nsrc; % This is only needed if removing the k=0
% term from the Fourier space sum

% TARGETS
if param.Nsrc == param.Ntar                               % -> Assign targets = sources
    domain.ptar = domain.psrc;
else
    domain.ptar = zeros(param.Ntar,2);                     % -> Place random targers
    domain.ptar(:,1) = param.Lx.*rand(param.Ntar,1) - param.Lx/2;
    domain.ptar(:,2) = param.Ly.*rand(param.Ntar,1) - param.Ly/2;
end

if param.Ntar == 1
%         domain.ptar = [0 0];
    domain.ptar = [0.1 3];
end

if param.Ntar == 2
    domain.ptar(1,:) = [0 0];
    domain.ptar(2,:) = [1 1];
end

if param.Ntar == 4
    domain.ptar(1,:) = [0 0];
    domain.ptar(2,:) = [1 1];
    domain.ptar(3,:) = [2 2];
    domain.ptar(4,:) = [3 3];
end

end