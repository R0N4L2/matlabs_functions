function uni=uniongen(cel)
    k=numel(cel);
    uni=cel{1};
    if k>1
        for j=2:k
            uni=union(uni,cel{j});
        end
    end
end