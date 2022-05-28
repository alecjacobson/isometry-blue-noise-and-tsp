%X = readDMAT('secord-X-4096-tsp-greedy-insertion.dmat');
%% upsample
%X = reshape(cat(1,X',0.5*(X+X([2:end 1],:))'),2,[])';
%X = [0 0;1 0;1 0.5;0.5 0.5;0.5 1;0 1];
X = [0 0;0.5 0.5;1 0];
% κ-curve to smoothly interpolate
tic;
[V,Q] = kappa_curve(X,'Closed',false);
toc
% Sample regularly. Could use _flat to reduce count
stride = 10;
t = linspace(0,1,stride);
Y = cell2mat(arrayfun(@(q) quadratic_eval(V(Q(q,:),:),t),1:size(Q,1),'UniformOutput',false)');

clf;
hold on;
plt(Y);
hold off;
axis equal;
set(gca,'Ydir','reverse');
