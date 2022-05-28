% outward angle
%rng(0);
%r = 0.005;
%X = blue_noise([0 0;1 1],r*100);
X = readDMAT('secord-X-4096-tsp-greedy-insertion.dmat');

Y = add_swirls(X,[]);


clf;
hold on;
%plt(X([1:end 1],:),'-k');
%sct(X,'k','MarkerFaceColor','k');
plt(Y([1:end 1],:),'-k','LineWidth',1);
hold off;
set(gca,'Ydir','reverse');
axis equal;
 
