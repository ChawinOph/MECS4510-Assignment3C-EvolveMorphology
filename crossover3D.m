function [childA childB] = crossover3D(parentA, parentB)
cut = [randi(4), randi(4); randi(3), randi(3); randi(5), randi(5)];
indices = 
childA = zeros(4,3,5);
childB = childA;
if cut(1,1) > cut(1,2)
    if cut(2,1) > cut(2,2)
        if cut(3,1) > cut(3,2)
            childA = 
            childB = 
        else
            childA = 
            childB =
        end
    else
        if cut(3,1) > cut(3,2)
            childA = 
            childB = 
        else
            childA = 
            childB =
        end
    end
else
    if cut(2,1) > cut(2,2)
        if cut(3,1) > cut(3,2)
            childA = 
            childB = 
        else
            childA = 
            childB =
        end
    else
        if cut(3,1) > cut(3,2)
            childA = 
            childB = 
        else
            childA = [parent
            childB =
        end
    end

end