function [childA, childB] = crossoverMat(parentA, parentB)
cut = [randperm(5,2); randperm(9,2)]
% cut = [4,1; 7,2]
% childA = zeros(4,3,5);
% childB = childA;
if cut(1,1) < cut(1,2)
    if cut(2,1) < cut(2,2)
        'case 1'
        childA = parentA;
        childA(cut(1,1):cut(1,2), cut(2,1):cut(2,2)) = parentB(cut(1,1):cut(1,2), cut(2,1):cut(2,2));
        childB = parentB;
        childB(cut(1,1):cut(1,2), cut(2,1):cut(2,2)) = parentA(cut(1,1):cut(1,2), cut(2,1):cut(2,2));
    else
        'case 2'
        childA = parentA;
        childA(cut(1,1):cut(1,2), 1:cut(2,2)) = parentB(cut(1,1):cut(1,2), 1:cut(2,2));
        childA(cut(1,1):cut(1,2), cut(2,1):end) = parentB(cut(1,1):cut(1,2), cut(2,1):end);
        childB = parentB;
        childB(cut(1,1):cut(1,2), 1:cut(2,2)) = parentA(cut(1,1):cut(1,2), 1:cut(2,2));
        childB(cut(1,1):cut(1,2), cut(2,1):end) = parentA(cut(1,1):cut(1,2), cut(2,1):end);
    end
else
    if cut(2,1) < cut(2,2)
        'case 3'
        childA = parentA;
        childA(1:cut(1,2), cut(2,1):cut(2,2)) = parentB(1:cut(1,2), cut(2,1):cut(2,2));
        childA(cut(1,1):end, cut(2,1):cut(2,2)) = parentB(cut(1,1):end, cut(2,1):cut(2,2));
        childB = parentB;
        childB(1:cut(1,2), cut(2,1):cut(2,2)) = parentA(1:cut(1,2), cut(2,1):cut(2,2));
        childB(cut(1,1):end, cut(2,1):cut(2,2)) = parentA(cut(1,1):end, cut(2,1):cut(2,2));
    else
        'case 4'
        childA = parentA;
        childA(1:cut(1,2), 1:cut(2,2)) = parentB(1:cut(1,2), 1:cut(2,2));
        childA(cut(1,1):end, cut(2,1):end) = parentB(cut(1,1):end, cut(2,1):end);
        childA(cut(1,1):end, 1:cut(2,2)) = parentB(cut(1,1):end, 1:cut(2,2));
        childA(1:cut(1,2), cut(2,1):end) = parentB(1:cut(1,2), cut(2,1):end);
        childB = parentB;
        childB(1:cut(1,2), 1:cut(2,2)) = parentA(1:cut(1,2), 1:cut(2,2));
        childB(cut(1,1):end, cut(2,1):end) = parentA(cut(1,1):end, cut(2,1):end);
        childB(cut(1,1):end, 1:cut(2,2)) = parentA(cut(1,1):end, 1:cut(2,2));
        childB(1:cut(1,2), cut(2,1):end) = parentA(1:cut(1,2), cut(2,1):end);
    end

end