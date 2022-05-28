function P = kappa_curve_demo(im)

  function P0 = current_point()
    P0 = get(gca,'currentpoint');
    P0 = P0(1,1:2);
  end
  % Callback for mouse press
  function onmove(src,ev)
    if ~closed && ~isempty(P)
      P0 = current_point();
      tmp_closed = false;
      if normrow(P0-P(1,:))<tol1
        tmp_closed = true;
        Q = P;
      else
        if size(P,1) >1
          if size(P,1)==2
            dir = P(2,:)-P(1,:);
          else
            [kV,kQ] = kappa_curve(P,'Closed',false);
            dir = (kV(kQ(end,3),:)-kV(kQ(end,2),:));
          end
          f = ease(normrow(P0-P(end,:))/tol);
          Psafe = P(end,:) + max(normrow(P0-P(end,:)),0.001)*dir;
          P0 = f*P0+(1-f)*Psafe;
        end
        Q = [P;P0];
      end
      draw_curve(Q,tmp_closed);
    end
  end
  function draw_curve(Q,closed)
    switch size(Q,1)
    case {0,1}
      % do nothing
      return;
    case 2
      % just draw straight line
      t = linspace(0,1)';
      X = Q(1,:)+t.*(Q(2,:)-Q(1,:));
    otherwise
      [kV,kQ] = kappa_curve(Q,'Closed',closed);
      t = linspace(0,1)';
      X = cell2mat(arrayfun(@(q) quadratic_eval(kV(kQ(q,:),:),t),1:size(kQ,1),'UniformOutput',false)');
    end
    set(xsh,'XData',X(:,1),'YData',X(:,2));
  end
  function ondown(src,ev)
    % Tell window that we'll handle drag and up events
    set(gcf,'windowbuttonmotionfcn', @ondrag);
    set(gcf,'windowbuttonupfcn',     @onup);
    P0 = current_point();
    if closed
      [~,sel] = min(normrow(P-P0));
    else
      if size(P,1)>=3 && normrow(P0-P(1,:))<tol1
        closed = true;
        set(gcf,'WindowButtonDownFcn', '');
        set(psh,'ButtonDownFcn', @ondown);
        sel = 1;
      else
        P = [P;P0];
        sel = size(P,1);
      end
    end
    draw_points();
  end
  % Callback for mouse press
  function ondrag(src,ev)
    P(sel,:) = current_point();
    draw_points();
    draw_curve(P,closed);
  end
  function draw_points()
    set(psh,'XData',P(:,1),'YData',P(:,2));
  end
  % Callback for mouse release
  function onup(src,ev)
    % Tell window to handle drag and up events itself
    set(gcf,'windowbuttonmotionfcn',@onmove);
    set(gcf,'windowbuttonupfcn','');
  end
  function onkeypress(src,ev)
    % escape character id
    ESC = {char(27),char(3)};
    ENTER = char(13);
    switch ev.Character
    case ESC
      finish();
    case ENTER
      finish();
    otherwise
      if ~isempty(ev.Character)
        warning(['Unknown key: ' ev.Character ...
          ' (' num2str(uint8(ev.Character(1))) ')']);
      end
    end
  end
  function finish()
    done = true;
    set(gca,'ButtonDownFcn','');
    set(gcf,'keypressfcn','');
    set(gcf,'windowbuttondownfcn','');
    set(gcf,'windowbuttonmotionfcn','');
    set(gcf,'windowbuttonupfcn','');
    delete(dsh);
  end
  
  P = [];
  sel = 0;
  closed = false;
  ease = @(t) 3.*t.^2-2.*t.^3;ease = @(t) ease(min(max(t,0),1));
  tol = 0.02;
  tol1 = 0.02;

  clf;
  axis([0 1 0 1]);
  axis equal;
  axis manual;
  set(gcf,'WindowButtonDownFcn',@ondown);
  set(gcf,'windowbuttonmotionfcn',@onmove);
  set(gcf,'keypressfcn',@onkeypress);

  % https://www.mathworks.com/matlabcentral/answers/397491-wait-a-boolean-variable-becomes-false#answer_374431
  hold on;
  if nargin>=1 
    [GX,GY] = meshgrid(0:1,[size(im,1)/size(im,2)]*[-0.5 0.5]+0.5);
    surf(GX,GY,0*GX,'CData',im,'FaceColor','texture','EdgeColor','none');
  end
  green = hex2rgb('0FFF50');
  xsh = plot(nan,nan,'-','Color',green,'LineWidth',2);
  psh = scatter(nan,nan,'MarkerFaceColor',0.0*[1 1 1],'MarkerEdgeColor',green,'SizeData',100,'LineWidth',2);
  dsh = plot(nan,nan);
  set(gcf,'Color',0.15*[1 1 1]);
  set(gca,'Visible','off','Position',[0 0 1 1]);
  hold off;
  waitfor(dsh)

end
