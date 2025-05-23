function [u,v,a] = globalK2D(EX,mu,elem,node)
% globel stiffness matrix
% for 2D problems
% based on VEM
% u,v,a for sparse matrix in MATLAB

sumElem = size(elem,1);
elemLen = cellfun('length',elem); 
nnz = sum((2*elemLen).^2);
u = zeros(nnz,1); % 
v = zeros(nnz,1);
a = zeros(nnz,1);
ia = 0;
for n = 1:sumElem
    
    coor = node(elem{n},:);
    K = elemKVEM(EX,mu,coor);  % element stiffness matrix 

    nodeID = elem{n};
    Ndof = 2*length(nodeID);
    ndfID = zeros(1,Ndof);
    ndfID(1:2:Ndof-1) = nodeID*2-1;
    ndfID(2:2:Ndof) = nodeID*2;
    ddd = repmat(ndfID,Ndof,1); 
    u(ia+1:ia+Ndof^2) = ddd(:);
    v(ia+1:ia+Ndof^2) = repmat(ndfID,1,Ndof);
    a(ia+1:ia+Ndof^2) = K(:);
    ia = ia + Ndof^2;
end










