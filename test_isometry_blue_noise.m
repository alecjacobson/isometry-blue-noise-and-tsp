
%im = imresize(im2double(imread('~/Downloads/alek-wek')),2^8*[1 1]);
%im = imresize(im2double(imread('~/Downloads/alek-wek-repeat.jpg')),2^8*[1 1]);
im = imresize(im2double(imread('~/Downloads/alek-wek-cropped.jpg')),2^8*[1 1]);
im = im(:,:,1);
[V,F] = create_regular_grid(size(im,2)+1,size(im,1)+1);
% image as scaling signal
S = reshape(repmat(im(:),1,2)',[],1);
% Conformal doesn't work... exactly because it _doesn't_ care about scale...
T = repmat( (1-S).^2,3,1);
T = T/mean(T);
U = V(F,:).*T;
G = reshape(1:numel(F),size(F));
P = sparse(G,F,1);

animated = true;
if animated
  VS = V;
  tsh = tsurf(F,VS,'CData',S,falpha(1,0));
  colormap(repmat(linspace(0,1,256),3,1)');
  axis equal;
  set(gca,'Ydir','reverse')
  drawnow;
  
  data = [];
  method = 'arap';max_iter = 20;
  method = 'ss';max_iter=1000;
  tic;
  VV = V;
  for iter = 1:max_iter
    switch method
    case 'arap'
      %[VS,data] = arap(VS,F,[],[],'V0',VS,'Energy','elements','Ref',{U,G,P'},'Data',data,'MaxIter',1);
      [VS,data] = arap(VS,F,[],[], ...
        'V0',VS,'Energy','elements','Ref',{U,G,P'},'Data',data,'MaxIter',1, ...
        'Dynamic',zeros(size(V)),'TimeStep',5e-1,'Vm1',VV(:,:,max(iter-1,1)));
    case 'ss'
      [VS,data] = smith_and_schaefer(U,F,VS,'Data',data,'MaxIters',1,'Tol',avgedge(V,F)/4,'PreventBoundaryOverlaps',false);
      if data.converged 
        break;
      end
    end
    tsh.Vertices = VS;
    VV = cat(3,VV,VS);
    drawnow;
  end
  toc
  
  
  r = 0.002;
  tic;
  %[bP,I,B] = blue_noise(VS,F,r);
  [BV,BF] = create_regular_grid(2,2);
  % Could use Bridson's method here
  BV = BV.*(max(VS)-min(VS))+min(VS);
  bP = blue_noise(VS,r*norm(max(VS)-min(VS)),'K',100);
  
  % This is slower than it needs to be:
  tic;
  [~,I,J,B] = in_element(VS,F,bP,'Method','ray-stabbing');
  
  D = squeeze(sqrt(sum(sum(diff(VV,[],3).^2,1),2)));
  C = cumsum([0;D]);
  C = C/C(end);
  ease = @(t) 3.*t.^2-2.*t.^3;ease = @(t) ease(min(max(t,0),1));
  tsh = tsurf(F,VS,'CData',S,falpha(1,0));
  colormap(repmat(linspace(0,1,256),3,1)');
  axis equal;
  set(gca,'Ydir','reverse')
  set(gca,'Visible','off','Position',[0 0 1 1]);
  set(gcf,'Color','w');
  
  for t = ease(linspace(0,1,30))
    Vt = permute(interp1(permute(C,[2 1 3]),permute(VV,[3 1 2]),t),[2 3 1]);
    tsh.Vertices = Vt;
    drawnow;
    %figgif('alek-wek-sample.gif','Delay',(t==0)*0.5);
  end
  sh = sct([nan nan],'.k','SizeData',1)
  set(gca,'Ydir','reverse')
  axis equal;
  set(gca,'Visible','off','Position',[0 0 1 1]);
  set(gcf,'Color','w');
  for t = ease(linspace(1,0,30))
    Vt = permute(interp1(permute(C,[2 1 3]),permute(VV,[3 1 2]),t),[2 3 1]);
    X = Vt(F(J,1),:).*B(:,1)+Vt(F(J,2),:).*B(:,2)+Vt(F(J,3),:).*B(:,3);
    set(sh,'Xdata',X(:,1),'YData',X(:,2));
    drawnow;
    %figgif('alek-wek-sample.gif','Delay',(t==1)*0.5 + (t==0)*0.5);
  end
else
  tic;
  VS = arap(V,F,[],[],'V0',V,'Energy','elements','Ref',{U,G,P'},'MaxIter',15);
  r = 0.001;
  tic;
  bP = blue_noise(VS,r*norm(max(VS)-min(VS)),'K',100);
  toc
  tic;
  [~,I,J,B] = in_element(VS,F,bP,'Method','ray-stabbing');
  X = V(F(J,1),:).*B(:,1)+V(F(J,2),:).*B(:,2)+V(F(J,3),:).*B(:,3);
  toc

  sct(X,'.k','SizeData',1);
  set(gca,'Ydir','reverse')
  axis equal;
  set(gca,'Visible','off','Position',[0 0 1 1]);
  set(gcf,'Color','w');
end
