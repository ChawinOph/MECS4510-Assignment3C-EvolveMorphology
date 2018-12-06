classdef MorphCube < handle
    %MorphCube 
    properties       
        chromosome = zeros(4,3,5)  % double 4 x 3 x 5 array: indirect encoding 
        masses    % masses of the robot
        springs   % springs of the robot 
        age  =  0 % age of the robot (no. of variation parameters gone)
    end
    
    properties (Constant)
         % global constants
         mass = 0.1; % m
         omega = pi; % (0.5 Hz of breathing);
         cube_length = 0.1 % double: length of each individual cube strcture
         voxel_dim = [5 5 5]   % int 3 x 1: no. of cubes in each dimension (default 5 x 5 x 5)
    end
    
    %%
    methods
        function obj = MorphCube(chromosomes, sorted_indcs, p_init_offset, v_init, a_init)
            %MorphCube Construct an instance of this class
            % test nargin > 0 every time you create ane object array
            % otherwise it won't work
            % https://www.mathworks.com/help/matlab/matlab_oop/initialize-object-arrays.html
            if nargin ~= 0
                switch nargin
                    % set the kinematic variables of all masses in each
                    % robot to zero there are missing inputs
                    case 1
                        sorted_indcs = 1:5;
                        p_init_offset = repmat([0 0 0], size(chromosomes, 4), 1);
                        v_init = repmat([0 0 0], size(chromosomes, 4), 1);
                        a_init = repmat([0 0 0], size(chromosomes, 4), 1);
                    case 2
                        p_init_offset = repmat([0 0 0], size(chromosomes, 4), 1);
                        v_init = repmat([0 0 0], size(chromosomes, 4), 1);
                        a_init = repmat([0 0 0], size(chromosomes, 4), 1);
                    case 3
                        v_init = repmat([0 0 0], size(chromosomes, 4), 1);
                        a_init = repmat([0 0 0], size(chromosomes, 4), 1);
                    case 4
                        a_init = repmat([0 0 0], size(chromosomes, 4), 1);
                end
                
                [Y, X, Z] = meshgrid(cube_length*(1: voxel_dim(1)), ...
                    cube_length*(1 : voxel_dim(2)), ...
                    cube_length*(1 : voxel_dim(3)));
                voxel_pos = [X(:), Y(:), Z(:)];
                
                % create an object arrays of robots
                for n = size(chromosomes, 4):-1:1
                    
                    obj(n).chromosome = chromosomes(:, :, :, n);
                    if nargin > 1
                        % sort rows of gene in case we do linkage tightening
                        % do the inverse operation of the sorting
                        % note: sorted_indcs is a 1d array
                        unsorting_indcs = sorted_indcs;
                        for i = 1:length(sorted_indcs)
                            unsorting_indcs(i) = find(sorted_indcs == i);
                        end
                        obj(n).chromosome = obj(n).chromosome(:, :, unsorting_indcs);
                    end
                  
                    % create the type of robot based on the pdf functions 
                    n_voxel = numel(X);
                    
                    v1 = mvnpdf(voxel_pos, [0.3 0.3 0.3], 0.25*eye(3));
                    v2 = mvnpdf(voxel_pos, [0.1 0.1 0.1], 0.25*eye(3));
                    
                    % visualize the output function
                    figure;
                    pcolor3(Y, X, Z, reshape(v1 + v2, 5,5,5)); hold on;
                    scatter3(X(:), Y(:), Z(:),[], v)
                    
                    % change the position and orientation of the robot after
                    % constructing the springs
                    % R = obj.rotationAxisAngle([1 0 0], pi/6); % tile around x axis by 30 degree
                    R = eye(3);
                    p = R*p'; % tilt all masses
                    p = p' + p_init_offset(n,:); % add the offset (each row of 2d array);
                    obj(n).masses = PointMass(repmat(mass, size(p,1), 1), p, v_init, a_init);
                    
                    
                end
            end
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
    
end

