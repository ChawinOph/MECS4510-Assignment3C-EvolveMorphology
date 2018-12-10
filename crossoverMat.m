function [childA, childB] = crossoverMat(parentA, parentB)
cut = [randperm(5,2); randperm(9,2)];
% cut = [4,1; 7,2]
% childA = zeros(4,3,5);
% childB = childA;
if cut(1,1) < cut(1,2)
    if cut(2,1) < cut(2,2)
        'case 1';
        childC = parentA;
        childC(cut(1,1):cut(1,2), cut(2,1):cut(2,2)) = parentB(cut(1,1):cut(1,2), cut(2,1):cut(2,2));
        childD = parentB;
        childD(cut(1,1):cut(1,2), cut(2,1):cut(2,2)) = parentA(cut(1,1):cut(1,2), cut(2,1):cut(2,2));
    else
        'case 2';
        childC = parentA;
        childC(cut(1,1):cut(1,2), 1:cut(2,2)) = parentB(cut(1,1):cut(1,2), 1:cut(2,2));
        childC(cut(1,1):cut(1,2), cut(2,1):end) = parentB(cut(1,1):cut(1,2), cut(2,1):end);
        childD = parentB;
        childD(cut(1,1):cut(1,2), 1:cut(2,2)) = parentA(cut(1,1):cut(1,2), 1:cut(2,2));
        childD(cut(1,1):cut(1,2), cut(2,1):end) = parentA(cut(1,1):cut(1,2), cut(2,1):end);
    end
else
    if cut(2,1) < cut(2,2)
        'case 3';
        childC = parentA;
        childC(1:cut(1,2), cut(2,1):cut(2,2)) = parentB(1:cut(1,2), cut(2,1):cut(2,2));
        childC(cut(1,1):end, cut(2,1):cut(2,2)) = parentB(cut(1,1):end, cut(2,1):cut(2,2));
        childD = parentB;
        childD(1:cut(1,2), cut(2,1):cut(2,2)) = parentA(1:cut(1,2), cut(2,1):cut(2,2));
        childD(cut(1,1):end, cut(2,1):cut(2,2)) = parentA(cut(1,1):end, cut(2,1):cut(2,2));
    else
        'case 4';
        childC = parentA;
        childC(1:cut(1,2), 1:cut(2,2)) = parentB(1:cut(1,2), 1:cut(2,2));
        childC(cut(1,1):end, cut(2,1):end) = parentB(cut(1,1):end, cut(2,1):end);
        childC(cut(1,1):end, 1:cut(2,2)) = parentB(cut(1,1):end, 1:cut(2,2));
        childC(1:cut(1,2), cut(2,1):end) = parentB(1:cut(1,2), cut(2,1):end);
        childD = parentB;
        childD(1:cut(1,2), 1:cut(2,2)) = parentA(1:cut(1,2), 1:cut(2,2));
        childD(cut(1,1):end, cut(2,1):end) = parentA(cut(1,1):end, cut(2,1):end);
        childD(cut(1,1):end, 1:cut(2,2)) = parentA(cut(1,1):end, 1:cut(2,2));
        childD(1:cut(1,2), cut(2,1):end) = parentA(1:cut(1,2), cut(2,1):end);
    end
end

if sum(sum(parentA==childC))>23
    childA = childC;
    childB = childD;
else
    childA = childD;
    childB = childC;
end

end