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
         voxel_dim = [4 4 4]   % int 3 x 1: no. of cubes in each dimension (default 5 x 5 x 5)
         var_range = [0.01 0.25; 0.01 0.25; 0.01 0.25];  % [Var_x; Var_y; Var_z]
         % Cov_xy = [0,1]*sqrt(Var_x)*sqrt(Var_y)
         % Cov_xz = [0,1]*sqrt(Var_x)*sqrt(Var_z)
         % Cov_yz = [0,1]*sqrt(Var_y)*sqrt(Var_z)
    end
    
    %%
    methods
        function obj = MorphCube(chromosomes, mat_sorted_indcs, p_init_offset, v_init, a_init)
            %MorphCube Construct an instance of this class
            % test nargin > 0 every time you create ane object array
            % otherwise it won't work
            % https://www.mathworks.com/help/matlab/matlab_oop/initialize-object-arrays.html
            if nargin ~= 0
                switch nargin
                    % set the kinematic variables of all masses in each
                    % robot to zero there are missing inputs
                    case 1
                        mat_sorted_indcs = 1:5;
                        p_init_offset = repmat([0 0 0], size(chromosomes, 3), 1);
                        v_init = repmat([0 0 0], size(chromosomes, 3), 1);
                        a_init = repmat([0 0 0], size(chromosomes, 3), 1);
                    case 2
                        p_init_offset = repmat([0 0 0], size(chromosomes, 3), 1);
                        v_init = repmat([0 0 0], size(chromosomes, 3), 1);
                        a_init = repmat([0 0 0], size(chromosomes, 3), 1);
                    case 3
                        v_init = repmat([0 0 0], size(chromosomes, 3), 1);
                        a_init = repmat([0 0 0], size(chromosomes, 3), 1);
                    case 4
                        a_init = repmat([0 0 0], size(chromosomes, 3), 1);
                end
                
                % create meshgrid of voxel centers
                [X, Y, Z] = meshgrid(obj.cube_length*(1: obj.voxel_dim(1)), ...
                    obj.cube_length*(1 : obj.voxel_dim(2)), ...
                    obj.cube_length*(1 : obj.voxel_dim(3)));
                % convert to 2d array for pdf calculation
                voxel_pos = [X(:), Y(:), Z(:)];
                % move the origin to the voxel center
                voxel_pos = voxel_pos - min(voxel_pos) - range(voxel_pos)/2;
                
                % create an object arrays of robots
                for n = size(chromosomes, 3):-1:1
                    
                    obj(n).chromosome = chromosomes(:, :, n);
                    if nargin > 1
                        % sort rows of gene in case we do linkage tightening
                        % do the inverse operation of the sorting
                        % note: sorted_indcs is a 1d array
                        unsorting_indcs = mat_sorted_indcs;
                        for i = 1:length(mat_sorted_indcs)
                            unsorting_indcs(i) = find(mat_sorted_indcs == i);
                        end
                        obj(n).chromosome = obj(n).chromosome(:, unsorting_indcs);
                    end       
                    
                    % initializa the prob density PD values as a n_voxel x
                    % n_matrl
                    pd = zeros(size(obj(n).chromosome, 1), size(voxel_pos, 1));
                                      
                    % convert normalized values of mu to the xyz positions                    
                    MU = obj(n).chromosome(:, 1:3).*range(voxel_pos) + min(voxel_pos);
                    
                    % calculate the PD for each material
                    for i = 1:size(obj(n).chromosome, 1)
                        % get variance values from normalized value
                        var = obj(n).chromosome(i, 4:6)'.*diff(obj(n).var_range, 1, 2) + obj(n).var_range(:, 1);
                        covar = diag(var);
                        % calculate covariance values from normalized value
                        covar(2,1) = obj(n).chromosome(i, 7)*sqrt(covar(1,1))*sqrt(covar(2,2));
                        covar(1,2) = covar(2,1);
                        covar(3,1) = obj(n).chromosome(i, 8)*sqrt(covar(1,1))*sqrt(covar(3,3));
                        covar(1,3) = covar(3,1);
                        covar(3,2) = obj(n).chromosome(i, 9)*sqrt(covar(2,2))*sqrt(covar(3,3));
                        covar(2,3) = covar(3,2);                      
                        % since covar(i,j) < max(var_range) by construction and a symmetric diagonally dominant matrix
                        % is symmetric positive definite, which can be ensured by adding nI
                        % https://math.stackexchange.com/questions/357980/how-to-generate-random-symmetric-positive-definite-matrices-using-matlab
                        covar = covar + diag(max(obj(n).var_range, [], 2));
                        pd(i, :) = mvnpdf(voxel_pos, MU(i, :), covar);
                    end
                    
                    % select material with prob in proportion to weights (PD) of
                    % material at each columbn of voxel coordinate (or just find the max value) 
%                     matrl = zeros(1, size(pd, 2));
%                     for i = 1:size(pd, 2)
%                         matrl(i) = randsample(1:size(pd, 1), 1, true,rand(1,5));
%                     end

                    [~, matrl] = max(pd); % max value case
                    % reshape back to voxel dimension
                    voxel_matrl = reshape(matrl, obj(n).voxel_dim(1), obj(n).voxel_dim(2), obj(n).voxel_dim(3));  
                    figure;
                    scatter3(voxel_pos(:,1), voxel_pos(:,2), voxel_pos(:,3), [], matrl, '.');
                    colorbar; colormap('jet')
                    axis equal;
                    
                    %% contruct the robot based on the choice of material                                       
                    
                    % change the position and orientation of the robot after
                    % constructing the springs
                    % R = obj.rotationAxisAngle([1 0 0], pi/6); % tile around x axis by 30 degree
%                     R = eye(3);
%                     p = R*p'; % tilt all masses
%                     p = p' + p_init_offset(n,:); % add the offset (each row of 2d array);
%                     obj(n).masses = PointMass(repmat(mass, size(p,1), 1), p, v_init, a_init);
                    
                    
                end
                
            end
        end
        
        function plotMaterialPDF(obj)
            %METHOD1 Summary of this method goes here
                                % visualize the output function
%                     figure;
%                     pcolor3(X, Y, Z, reshape(v1 + v2, 5,5,5));
%                     colorbar
%                     figure;
%                     pcolor3(X, Y, Z, reshape(v1, 5,5,5));
%                     colorbar
%                     figure;
%                     pcolor3(X, Y, Z, reshape(v2, 5,5,5));
%                     colorbar
%                     scatter3(X(:), Y(:), Z(:),[], v1); hold on;
%                     scatter3(X(:), Y(:), Z(:),[], matrl, 'x');
%                     colorbar; colormap('jet')
                    
%                     figure;
%                     v1_2d = reshape(v1, 5,5,5); v1_2d = v1_2d(:,:,1);
%                     surf(X(:,:,1), Y(:,:,1), v1_2d)
        end
    end
    
end

