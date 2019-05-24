function ue=promcond(ss)
    ut=unique(ss(:,1))';
    ue=[];
    for u=ut
        ue=[ue;[u,mean(ss(ss(:,1)==u,2))]];
    end
end