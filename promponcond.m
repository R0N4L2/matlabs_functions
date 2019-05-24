function ue=promponcond(ss,tt)
    ut=union(ss(:,1),tt(:,1))';
    ue=[];
    for u=ut
        a=ss(ss(:,1)==u,2);
        b=tt(tt(:,1)==u,2);
        if size(a,1)<size(b,1)
            a=[a;zeros(size(b,1)-size(a,1),1)];
        elseif size(a,1)>size(b,1)
            b=[b;zeros(size(a,1)-size(b,1),1)];
        end
        ue=[ue;[u,sum(a.*b)/sum(b)]];
    end
end