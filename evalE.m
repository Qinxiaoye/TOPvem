function ce=evalE(u,GK_a,sumElem,elem)
ce=zeros(sumElem,1);
ia=0;
for n = 1:sumElem
    nodeID = elem{n};
    Ndof = 2*length(nodeID);
    ue=zeros(Ndof,1);
    Ke0=GK_a(ia+1:ia+Ndof^2);
    Ke=reshape(Ke0,Ndof,[]);
    ia = ia + Ndof^2;
    for in=1:Ndof/2
        ue(in*2-1:2*in)=u(nodeID(in)*2-1:nodeID(in)*2);
    end
    ce(n)=ue'*Ke*ue;
end