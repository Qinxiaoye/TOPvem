function [neibor,weightE]=filterE(elem,node,sumElem,rmin,ndim)

%---------------------------各单元中心坐标----------------------------------
elemCD=zeros(sumElem,ndim);
for i=1:sumElem
    nodeID = elem{i};
    mnode = length(nodeID);
    for n=1:mnode
        elemCD(i,:)=elemCD(i,:)+node(nodeID(n),:)/mnode;
    end
end
%--------------------------各单元附近的单元---------------------------------
neibor=zeros(sumElem,sumElem+1);
weightE=zeros(sumElem,sumElem);
for i=1:sumElem
    tt=0;
    for j=1:sumElem
        dist=norm(elemCD(j,:)-elemCD(i,:));
        if dist<=rmin
            tt=tt+1;
            neibor(i,tt+1)=j;
            weightE(i,tt)=rmin-dist;
        end
    end
    neibor(i,1)=tt;
end
            
    