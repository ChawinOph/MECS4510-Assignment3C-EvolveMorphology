function mutant = mutate2(gene)
%This function takes a gene sequence and increments a row of values by + or - 0.1
    mutant= gene;
    row = randi(size(gene,1));
    for i = 1:size(gene,2)
        flip = rand;
        if flip>0.5
            mutant(row,i) = mutant(row,i)+ 0.1;
        else
            mutant(row,i) = mutant(row,i)- 0.1;
        end
    end
end