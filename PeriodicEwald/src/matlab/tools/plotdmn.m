function plotdmn(pt,ps,Lx,Ly,figind)
% This function plots the domain sources (ps) and targets (pt) with the reference
% cell boundaries
figure(figind);clf;
p1= plot(ps(:,1),ps(:,2),'k.','MarkerSize',20);
hold on
p2 = plot(pt(:,1),pt(:,2),'rx');
plot(Lx/2*[-1 -1],Ly/2*[-1 1],'k-');
plot(Lx/2*[1 1],Ly/2*[-1 1],'k-');
plot(Lx/2*[-1 1],Ly/2*[1 1],'k-');
plot(Lx/2*[-1 1],Ly/2*[-1 -1],'k-');

legend([p1 p2],'sources','targets')

axis equal
end