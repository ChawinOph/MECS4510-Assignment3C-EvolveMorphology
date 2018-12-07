function mutant = mutate1(gene)
%This function takes a gene sequence and changes one of the values to a new
%random value
    mutant= gene;
    mutant(randi(numel(gene))) = rand;
end