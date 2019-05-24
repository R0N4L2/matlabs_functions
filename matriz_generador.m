function [K,D,A,B,E]=matriz_generador(M,carne,peso)
P=M(2:end,:);
P(isnan(P))=0;
n=size(P,1);
P=sparse(P);
R=P;
R(R<0)=0;
Q=R;
a=find(P(1,:)==0,1,'first');
for j=2:n-3*strcmp(carne,'pA')
    A=nchoosek(2:n,j);
    if strcmp(carne,'pA')
        %restriccion que no se pueden combinar tomahauk con rib steak o
        %club steak
        C=A((all(A~=9,2)|all(A<10|A>11,2))&(all(A~=7,2)|all(A<8|A>11,2)),:);
    else
        C=A;
    end  
    if ~isempty(C)
        i1=arrayfun(@(x)M(1,find(P(x,a:end)~=0)+a-1),C,'UniformOutput',0);
        S=mat2cell(C,ones(1,size(C,1)),size(C,2));
        i2=cell2mat(cellfun(@(x)sum(sum(Q(x,a:end),1)),S,'UniformOutput',0));
        i3=cell2mat(cellfun(@(x)sum(Q(1,ismember(M(1,1:a-1),setdiff(M(1,find(any(P(x,a:end)~=0,1))+a-1),[15273,15274])))),S,'UniformOutput',0));
        ie=abs(i3-i2)<1e-5;
        C=C(ie,:);
        if ~isempty(C)
            i1=i1(ie,:);
            S=mat2cell(C,ones(1,size(C,1)),size(C,2));
            Z=sparse(cell2mat(cellfun(@(x)[full(Q(1,1:a-1)),sum(full(Q(x,a:end)),1)],S,'UniformOutput',0)));
            ue=arrayfun(@(x)find(ismember(M(1,1:a-1),uniongen(i1(x,:)))),1:size(i1,1),'UniformOutput',0)';
            row=find(cellfun(@(x)~isempty(x),ue));
            for p=row'
                Z(p,ue{p})=0;
            end
            Z(Z<0)=0;
            R=[R;Z];
            R=unique(R(abs(sum(R,2)-peso)<1.e-4,:),'rows');
        end
    end
end
f=find(any(P>0,1));
R=R(:,f);
b=P(1,:)==0;
[n,m]=size(R);
D=zeros(n,m*n);
E=zeros(n-1,m*n);
for k=1:n
    D(k,(1:m)+(k-1)*m)=full(R(k,:));
    E(k,(1:m)+(k-1)*m)=[zeros(1,a-1),R(k,a:end)];
end
t=sum(E,2);
E=sparse([E,-diag(t)]);
E=E(t>0,:);
D=sparse([D,-peso*eye(n)]);
id=any(D~=0,1);
A=[repmat(b(1,f),[1,n]),zeros(1,n)];
K=[repmat(M(1,f),[1,n]),-(1:n)-(strcmp(carne,'pA')+2*strcmp(carne,'pB')+3*strcmp(carne,'bA'))*10^(2+ceil(log10(n))+(floor(log10(n))==ceil(log10(n))))];
A=A(id);
K=K(id);
E=E(:,id);
D=D(:,id);
B=[M(1,f);R];
end


