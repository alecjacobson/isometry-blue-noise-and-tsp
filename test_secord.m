oim = im2double(imread('~/Downloads/alek-wek-cropped.jpg'));
oim = oim(:,:,1);
oim = imresize(oim,2);
%oim = repmat(linspace(0,1,2^9),2^9,1);
%ease = @(t) 3.*t.^2-2.*t.^3;ease = @(t) ease(min(max(t,0),1));
%oim = ease(ease(oim))  *0.99+0.005;

nearest_power_of_two = @(x) 2.^ceil(log(x)/log(2));
N = numel(oim)/50
N = nearest_power_of_two(N)
rng(0);

levels = log(N)/log(2)-log(512)/log(2);
max_iter = 2000;
tol = 0.01;

F = inf;

tic;
for level = levels-1:-1:0
  im = imresize(oim,1/(2^(level)));

  [pX,pY] = meshgrid(linspace(0,1,size(im,2)),linspace(0,1*(size(im,1)-1)/(size(im,2)-1),size(im,1)));
  Y = [pX(:) pY(:)];
  rho = max((1-im(:)),0);
  % Precompute
  
  % min ∑ⱼ minᵢ ρⱼ ‖ Xᵢ - xⱼ ‖²
  %
  % or 
  %
  % min ∑ⱼ ρⱼ ‖ X_ℓⱼ - xⱼ ‖²
  %
  % Suppose Xᵢ are fixed. The only variables are ℓⱼ for each pixel, choose the
  % closest one regardless of ρ...
  %
  
  if level == levels-1
    %n = nearest_power_of_two(numel(pY)/50)
    n = nearest_power_of_two(numel(oim)/100)/(4^(levels-1));
    %% white noise initialization
    %X = rand(n,2).*(max(Y)-min(Y))+min(Y);
    % Sampled from rho
    X = datasample(Y,n,'Replace',false,'Weights',rho);
  else
    Xold = X;
    bad = any(Xold<min(Y),2) | any(Xold>max(Y),2);
    fprintf('  #bad: %d\n',sum(bad));
    Xold(bad,:) = datasample(Y,sum(bad),'Replace',false,'Weights',rho);
    X = repmat(Xold,4,1);
    counts = accumarray(I,1,[size(Xold,1) 1]);
    %X = X+(rand(size(X))-0.5)*(1/size(im,2))*1e-3;
    %X = X+(rand(size(X))-0.5)*(1/size(im,2)).*repmat(max(sqrt(counts),1),4,1);
    th = reshape(rand(size(Xold,1),1)*2*pi + repmat([0 0.5 1 1.5]*pi,size(Xold,1),1),size(X,1),1);
    X = X+(1/size(im,2)).*0.5*repmat(max(sqrt(counts),1),4,1).*[cos(th) sin(th)];

    %R = randperm(size(Y,1));
    %I = knnsearch(Xold,Y(R,:));
    %L = sparse(1:numel(I),I,1);
    %J = zeros(size(L,2),4);
    %Z = zeros(size(L,2),1);
    %for l = 1:size(J,2)
    %  [Zl,J(:,l)] = max(L);
    %  Z = max(Z,Zl');
    %  if l>1
    %    J(Zl==0,l) = J(Zl==0,l-1);
    %  end
    %  L(sub2ind(size(L),J(:,l),(1:size(L,2))')) = 0;
    %end
    %J = R(J);
    %X = Y(J,:);
    X = min(max(X,min(Y)-0.1),max(Y)+0.1);

    I = knnsearch(Xold,Y);
    %ish = surf(pX,pY,pX*0,'CData',reshape(I,size(pX)),falpha(0.1,0));
    %set(gca,'Ydir','reverse');
    %hold on;
    %%qvr(X,repmat(Xold,size(J,2),1)-X,0);
    %sct(Xold,'.r');
    %sct(X,'.k');
    %hold off;
    %view(2);
    %axis equal;
    %drawnow
  end
    
  method = 'gradient-descent'; alpha=1;
  method = 'adam'; alpha = 2e0;
  method = 'sobolev'; alpha=1e-2;
  method = 'lloyds';
  method = 'lbfgs'; alpha = 1;
  % lbfgs state
  M = 8;
  k = 1;
  % rolling tape of M steps, 
  x = nan([size(X) M]);
  % rolling tape of M gradients
  g = nan([size(X) M]);
  % s(:,:, mod(k-1,M)+1 ) = x(:,:,mod(k+1-1,M)+1)-x(:,:,mod(k-1,M)+1)
  s = nan([size(X) M]);
  % rolling tape of M changes in gradient
  % y(:,:, mod(k-1,M)+1 ) = g(:,:,mod(k+1-1,M)+1) - g(:,:,mod(k-1,M)+1)
  y = nan([size(X) M]);
  % rolling tape of p (ρ) and gamma
  p = nan(M,1);
  vec = @(X) X(:);
  inner = @(X,Y) sum(X(:).*Y(:));
  %gamma = nan(M,1);
  
  with_line_search = false;
  with_clamping = true;
  % adam state
  m = 0;
  v = 0;
  beta1 = 0.9;
  beta2 = 0.999;
  epsilon = 1e-8;
  t = 0;
  
  for iter = 1:max_iter
    switch method
    case 'lloyds'
      X0 = X;
      [I,D] = knnsearch(X0,Y,'K',1);
      C = [accumarray(I,rho.*Y(:,1),[n 1]) accumarray(I,rho.*Y(:,2),[n 1])]./accumarray(I,rho,[n 1]);
      X = C;
      F0 = sum(rho.*D.^2);
      fprintf('%s: $%g @\n',method,F0,toc);
    case {'adam','gradient-descent','lbfgs','sobolev'}
      [F,dFdX,I,D,C] = voronoi_objective(X,Y,rho);
      X0 = X;
      F0 = F;
      dFdX0 = dFdX;
      dFdX0(isnan(dFdX0)) = 0;
      %dFdX0 = (X-C);
  
      switch method
      case 'adam'
        t = t+1;
        m0 = m;
        v0 = v;
        g = dFdX0;
        m = beta1 * m0 + (1-beta1)*g;
        v = beta2 * v0 + (1-beta2)*g.^2;
        mhat = m./(1-beta1^t);
        vhat = v./(1-beta2^t);
        dX = -mhat./(sqrt(vhat) + epsilon);
      case 'gradient-descent'
        dX = -dFdX0;
      case 'lbfgs'
        x(:,:,mod(k-1,M)+1) = X0;
        g(:,:,mod(k-1,M)+1) = dFdX0;
        if k>1
          % now we have enough info to update sk-1, yk-1, pk-1
          s(:,:, mod(k-1-1,M)+1 ) = x(:,:,mod(k-1,M)+1)-x(:,:,mod(k-1-1,M)+1);
          y(:,:, mod(k-1-1,M)+1 ) = g(:,:,mod(k-1,M)+1)-g(:,:,mod(k-1-1,M)+1);
          p(mod(k-1-1,M)+1) = 1./inner(y(:,:,mod(k-1-1,M)+1),s(:,:,mod(k-1-1,M)+1));
          % initialize
          r = -g(:,:,mod(k-1,M)+1);
          % 1st LBFGS update
          gamma = nan(M,1);
          for i = mod((k-1:-1:max(k-M,1))-1,M)+1
            gamma(i) = p(i)*inner(s(:,:,i),r);
            r = r-gamma(i)*y(:,:,i);
          end
          dk = inner(s(:,:,mod(k-1-1,M)+1),y(:,:,mod(k-1-1,M)+1)) / ...
            inner(y(:,:,mod(k-1-1,M)+1),y(:,:,mod(k-1-1,M)+1)) * r;
          % 2nd LBFGS update
          for i = flip(mod((k-1:-1:max(k-M,1))-1,M)+1)
            dk = dk + s(:,:,i)*(gamma(i)-p(i)*inner(y(:,:,i),dk));
          end
          dX = dk;
        else
          dX = -1e-3*g(:,:,mod(k-1,M)+1);
        end
        k = k + 1;
      case 'sobolev'
        [XI,DX] = knnsearch(X,X,'K',6+1);
        LJ = XI(:,2:end);
        LI = repmat(1:size(X,1),size(LJ,2),1)';
        LV = -1./DX(:,2:end);
        L = sparse([LI LJ LI LJ],[LJ LI LI LJ],[LV LV -LV -LV]);
        dX = (L+1e-3*speye(size(L)))\-dFdX0;
      end
      if with_line_search
        % This is quite expensive
        [alpha,X,F] = backtracking_line_search( ...
          @(X) voronoi_objective(X,Y,rho), ...
          X0,dFdX0,dX,0.3,0.5);
        if alpha < eps
          fprintf('Line search stalled.\n');
          break;
        end
        fprintf('$%g → $%g via %g\n',F0,F,alpha);
      else
        %fprintf('%s: $%g\n',method,F0);
        X = X0+alpha*dX;
      end

      if norm(dFdX0(:)) < tol
        break;
      end

      if with_clamping
        X = min(max(X,min(Y)-0.1),max(Y)+0.1);
      end
    end
  
    %clf;
    %%ish = surf(pX,pY,pX*0,'CData',reshape(I,size(pX)),falpha(0.1,0));
    %set(gca,'Ydir','reverse');
    %hold on;
    %sct(X0,'.r','SizeData',1);
    %sct(X,'.k','SizeData',3);
    %hold off;
    %view(2);
    %set(gca,'Visible','off','Position',[0 0 1 1])
    %set(gcf,'Color','w');
    %axis equal;
    %drawnow
    %figgif('secord.gif');
  end
  fprintf('%s: %d: %d: $%g @%g secs\n',method,level,iter,F0,toc);

end
toc
