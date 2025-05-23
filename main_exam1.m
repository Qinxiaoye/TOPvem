%%%% A 99 LINE TOPOLOGY OPTIMIZATION CODE WITH VEM  %%%
% see reference: Sigmund, O. A 99 line topology optimization code written in Matlab. https://doi.org/10.1007/s001580050176

clear;
volfrac = 0.5;
penal = 3;
rmin = 0.5;
loop = 0; 
change = 1.;

figure('color',[1 1 1]);
%-------------------------------parameter------------------------------------
ndim =2;
load mesh2.mat; % input mesh

sumNode = size(node,1);
sumElem = size(elem,1);
EX = 200000; mu = 0.3;
%boundary conditions
nodeL = find(node(:,1)<0.001); % nodes on the left side
nodeR = find(node(:,1)>20-0.001); % nodes on the right side
nodeRD = nodeR(node(nodeR,2)<0.001); % bottom right node
sL = size(nodeL,1);
fixNode = [nodeL,ones(sL,1),zeros(sL,1);...
        nodeRD,2,0];

% load
nodeLT = nodeL(node(nodeL,2)>10-0.001); % Apply load at the upper left corner
nodeForce = [nodeLT,2,-100]; % nodeID, load direction(1-x, 2-y), load value
%----------------calculate stiffness matrix with sparse form---------------
[GK_u,GK_v,GK_a] = globalK2D(EX,mu,elem,node);
%
Emin=1e-9;
E0=1;
%
[neibor,weightE]=filterE(elem,node,sumElem,rmin,ndim);
neiborNum=neibor(:,1);
neibor(:,1)=[];
% initial value
x=zeros(sumElem,1);
x(:)=volfrac;

% loop

while change > 0.05||loop<=100
    loop = loop + 1;
    xold = x;
    xa=x.^penal;


    [GK,force] = boundary_chang_one(xa,E0,Emin,GK_u,GK_v,GK_a,elem,fixNode,nodeForce,sumNode,sumElem,ndim);
    % solve equation
    u = GK\force;
    ce=evalE(u,GK_a,sumElem,elem);
    c = sum(sum((Emin+x.^penal*(E0-Emin)).*ce));
    dc = -penal*(E0-Emin)*x.^(penal-1).*ce;
    %
    cr_dc=zeros(sumElem,1);
    for i=1:sumElem
        cr_dcA=0;
        cr_dcB=0;
        for j=1:neiborNum(i)
            cr_dcA=cr_dcA+dc(neibor(i,j))*x(neibor(i,j))*weightE(i,j);
            cr_dcB=cr_dcB+x(i)*weightE(i,j);
        end
        cr_dc(i)=cr_dcA/cr_dcB;
    end
    % oc
    [x]    = OC(sumElem,x,volfrac,cr_dc);
    % show solution
    change = max(max(abs(x-xold)));
    disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
        ' Vol.: ' sprintf('%6.3f',sum(sum(x))/sumElem) ...
        ' ch.: ' sprintf('%6.3f',change )])

    colormap(flip(gray));
    max_n_vertices = max(cellfun(@length, elem));
    padding_func = @(vertex_ind) [vertex_ind,...
        NaN(1,max_n_vertices-length(vertex_ind))];  % function to pad the vacancies
    tpad = cellfun(padding_func, elem, 'UniformOutput', false);
    tpad = vertcat(tpad{:});
    if loop>1
        delete(h)
    end
    h=patch('Faces', tpad,...
        'Vertices', node,'EdgeColor','k',...
        'FaceColor', 'flat',...
        'FaceVertexCData', x);
    drawnow
    axis equal; axis off;


    Y(loop)=c; Vol(loop)=sum(sum(x))/(sumElem);
end

%--------------------------------Show iteration history------------------------------

figure;[AX,h1,h2]=plotyy(1:loop,Y,1:loop,Vol);legend([h1,h2],'Objective Function','Volume fraction','Location','Best');
set(h1,'linestyle','-','LineWidth',1.5,'color','r');
set(h2,'linestyle','-.','LineWidth',1.5,'color','b');
set(get(AX(1),'ylabel'),'string','Objective Function');set(get(AX(2),'ylabel'),'string','Volume fraction');xlabel('Iter.');
set(AX(2),'ylim',[0 1]);set(AX(2),'ytick',0:.25:1);


%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xnew]=OC(sumElem,x,volfrac,dc)  
l1 = 0; l2 = 100000; move = 0.05;
while (l2-l1 > 1e-4)
  lmid = 0.5*(l2+l1);
  xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*sqrt(-dc./lmid)))));  
  if sum(sum(xnew)) - volfrac*sumElem > 0
    l1 = lmid;
  else
    l2 = lmid;
  end
end
end
%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcn]=check(nelx,nely,rmin,x,dc)
dcn=zeros(nely,nelx);
for i = 1:nelx
  for j = 1:nely
    sum=0.0; 
    for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
      for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
        fac = rmin-sqrt((i-k)^2+(j-l)^2);
        sum = sum+max(0,fac);
        dcn(j,i) = dcn(j,i) + max(0,fac)*x(l,k)*dc(l,k);
      end
    end
    dcn(j,i) = dcn(j,i)/(x(j,i)*sum);
  end
end
end