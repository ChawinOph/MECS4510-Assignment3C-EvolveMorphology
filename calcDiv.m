function div = calcDiv(genes)
% This function calculates the sum of standard deviation of the genes of a
% population.
divMat = std(genes, 0, 3);
div = sum(divMat, 'all');

end