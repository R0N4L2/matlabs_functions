function inter=intergen(cel)
    k=numel(cel);
    inter=cel{1};
    if k>1
        for j=2:k
            inter=intersect(inter,cel{j});
        end
    end
end