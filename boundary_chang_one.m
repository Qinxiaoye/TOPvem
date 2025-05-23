function [GK,force] = boundary_chang_one(xa,E0,Emin,GK_u,GK_v,GK_a,elem,fixNode,nodeForce,sumNode,sumElem,ndim)
% first type boundary condition
% You can rewrite this code in a more concise form

fixNodeNew   = (fixNode(:,1)-1)*ndim+fixNode(:,2);
if ~isempty(nodeForce)
    nodeForceNew = (nodeForce(:,1)-1)*ndim+nodeForce(:,2);
else
    nodeForceNew = [];
end


sumForce = size(nodeForceNew,1);

a = (1:sumNode*ndim)';
c = ones(sumNode*ndim,1);
c(fixNodeNew) = 0;
b = sparse(a,a,c,sumNode*ndim,sumNode*ndim);

elemLen = cellfun('length',elem); 
nnz = sum((2*elemLen).^2);
GK_a2 = zeros(nnz,1);
ia=0;
for n = 1:sumElem    
    nodeID = elem{n};
    Ndof = 2*length(nodeID);
    GK_a2(ia+1:ia+Ndof^2)=GK_a(ia+1:ia+Ndof^2).*(Emin+xa(n)*(E0-Emin));
    ia = ia + Ndof^2;
end


GK = sparse(GK_u,GK_v,GK_a2,sumNode*ndim,sumNode*ndim);
GK = b*GK*b;  
GK = GK  -(b-speye(sumNode*ndim));


if size(nodeForce,1) > 0
    force = sparse(nodeForceNew,ones(sumForce,1),nodeForce(:,3),sumNode*ndim,1);%+sparse(bf);
    force = b*force;
else
    force = sparse(sumNode*ndim,1);%+sparse(bf);
    force = b*force;
end

