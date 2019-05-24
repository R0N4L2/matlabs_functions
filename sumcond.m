function ue=sumcond(ss)
    ut=unique(ss(:,1))';
    ue=[];
    for u=ut
        ue=[ue;[u,sum(ss(ss(:,1)==u,2))]];
    end
end