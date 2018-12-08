classdef MorphCube < handle
    %MorphCube
    properties
        chromosome = zeros(5,9)  % double 5 x 9 array: indirect encoding
        masses    % masses of the robot
        springs   % springs of the robot
        age  % age of the robot (no. of variation parameters gone)
        voxel_pos %
        voxel_matrl %
        pd %
    end
    
    properties (Constant)
        % global constants
        mass = 0.1; % m
        omega = pi; % (0.5 Hz of breathing);
        cube_length = 0.1 % double: length of each individual cube strcture
        voxel_dim = [2 2 2]   % int 3 x 1: no. of cubes in each dimension (default 5 x 5 x 5)
        var_range = [0.01 0.04; 0.01 0.04; 0.01 0.04];  % [Var_x; Var_y; Var_z]
        % Cov_xy = [0,1]*sqrt(Var_x)*sqrt(Var_y)
        % Cov_xz = [0,1]*sqrt(Var_x)*sqrt(Var_z)
        % Cov_yz = [0,1]*sqrt(Var_y)*sqrt(Var_z)
    end
    
    %% member functions
    methods
        function obj = MorphCube(chromosomes, age, mat_sorted_indcs, p_init_offset, v_init, a_init)
            %MorphCube Construct an instance of this class
            % test nargin > 0 every time you create ane object array
            % otherwise it won't work
            % https://www.mathworks.com/help/matlab/matlab_oop/initialize-object-arrays.html
            if nargin ~= 0
                switch nargin
                    % set the kinematic variables of all masses in each
                    % robot to zero there are missing inputs
                    case 1
                        age = 0;
                        mat_sorted_indcs = 1:5;
                        p_init_offset = repmat([0 0 obj.voxel_dim(3)*obj.cube_length/2], size(chromosomes, 3), 1);
                        v_init = repmat([0 0 0], size(chromosomes, 3), 1);
                        a_init = repmat([0 0 0], size(chromosomes, 3), 1);
                    case 2
                        mat_sorted_indcs = 1:5;
                        p_init_offset = repmat([0 0 obj.voxel_dim(3)*obj.cube_length/2], size(chromosomes, 3), 1);
                        v_init = repmat([0 0 0], size(chromosomes, 3), 1);
                        a_init = repmat([0 0 0], size(chromosomes, 3), 1);
                    case 3
                        p_init_offset = repmat([0 0 0], size(chromosomes, 3), 1);
                        v_init = repmat([0 0 0], size(chromosomes, 3), 1);
                        a_init = repmat([0 0 0], size(chromosomes, 3), 1);
                    case 4
                        v_init = repmat([0 0 0], size(chromosomes, 3), 1);
                        a_init = repmat([0 0 0], size(chromosomes, 3), 1);
                    case 5
                        a_init = repmat([0 0 0], size(chromosomes, 3), 1);
                end
                
                % create meshgrid of voxel centers
                [X, Y, Z] = meshgrid(obj.cube_length*(1: obj.voxel_dim(1)), ...
                    obj.cube_length*(1 : obj.voxel_dim(2)), ...
                    obj.cube_length*(1 : obj.voxel_dim(3)));
                
                % create an object arrays of robots
                for n_bot = size(chromosomes, 3):-1:1
                    
                    % convert to 2d array for pdf calculation
                    obj(n_bot).voxel_pos = [X(:), Y(:), Z(:)];
                    % move the origin to the voxel center
                    obj(n_bot).voxel_pos = obj(n_bot).voxel_pos - min(obj(n_bot).voxel_pos) - range(obj(n_bot).voxel_pos)/2;
                    obj(n_bot).chromosome = chromosomes(:, :, n_bot);
                    
                    if nargin > 1
                        % sort rows of gene in case we do linkage tightening
                        % do the inverse operation of the sorting
                        % note: sorted_indcs is a 1d array
                        unsorting_indcs = mat_sorted_indcs;
                        for i = 1:length(mat_sorted_indcs)
                            unsorting_indcs(i) = find(mat_sorted_indcs == i);
                        end
                        obj(n_bot).chromosome = obj(n_bot).chromosome(:, unsorting_indcs);
                    end
                    
                    % initializa the prob density PD values as a n_voxel x
                    % n_matrl
                    pd = zeros(size(obj(n_bot).chromosome, 1), size(obj(n_bot).voxel_pos, 1));
                    
                    % convert normalized values of mu to the xyz positions
                    MU = obj(n_bot).chromosome(:, 1:3).*range(obj(n_bot).voxel_pos) + min(obj(n_bot).voxel_pos);
                    
                    % calculate the PD for each material
                    for i = 1:size(obj(n_bot).chromosome, 1)
                        
                        % get variance values from normalized value
                        var = obj(n_bot).chromosome(i, 4:6)'.*diff(obj(n_bot).var_range, 1, 2) + obj(n_bot).var_range(:, 1);
                        covar = diag(var);
                        
                        % calculate covariance values from normalized value
                        covar(2,1) = obj(n_bot).chromosome(i, 7)*sqrt(covar(1,1))*sqrt(covar(2,2));
                        covar(1,2) = covar(2,1);
                        covar(3,1) = obj(n_bot).chromosome(i, 8)*sqrt(covar(1,1))*sqrt(covar(3,3));
                        covar(1,3) = covar(3,1);
                        covar(3,2) = obj(n_bot).chromosome(i, 9)*sqrt(covar(2,2))*sqrt(covar(3,3));
                        covar(2,3) = covar(3,2);
                        
                        % since covar(i,j) < max(var_range) by construction
                        % and a symmetric diagonally dominant matrix is
                        % symmetric positive definite, which can be ensured
                        % by adding nI
                        % https://math.stackexchange.com/questions/357980/
                        % how-to-generate-random-symmetric-positive-definite-matrices-using-matlab
                        covar = covar + diag(max(obj(n_bot).var_range, [], 2));
                        pd(i, :) = mvnpdf(obj(n_bot).voxel_pos, MU(i, :), covar);
                        
                    end
                    
                    obj(n_bot).pd = pd;
                    
                    % select material with prob in proportion to weights (PD) of
                    % material at each columbn of voxel coordinate (or just find the max value)
                    % matrl = zeros(1, size(pd, 2));
                    % for i = 1:size(pd, 2)
                    %    matrl(i) = randsample(1:size(pd, 1), 1, true,rand(1,5));
                    % end
                    
                    [~, matrl] = max(pd); % max value case
                    
                    % reshape back to voxel dimension and store in
                    % properties
                    obj(n_bot).voxel_matrl = reshape(matrl, obj(n_bot).voxel_dim(1), obj(n_bot).voxel_dim(2), obj(n_bot).voxel_dim(3));
                    
                    %% check connectivity of the cubes (ignore it for now)
                    % if there are separate bodies, choose the body with
                    % most masses.
                    % If the number of masses are equal, keep only the
                    % masse that is found first strting from (0,0,0) pos
                    
                    % check the matrial cube if it could go to every other
                    % material cube in the voxel space
                    % find all checking pairs of the voxel
                    
                    % construct the adjacency matrix of 6 connected
                    % neighbors (top-bottom-N-W-S-E)
                    
                    %% construct the cube
                    
                    % create half diagonal vectors from the vox center
                    cube_corner = obj(n_bot).cube_length*ones(1,3)/2;
                    corner_signs =  [1, 1, 1;   % octant 1
                        -1  1  1;   % octant 2
                        -1 -1  1;   % octant 3
                        1 -1  1;   % octant 4
                        1  1 -1;   % octant 5
                        -1  1 -1;   % octant 6
                        -1 -1 -1;   % octant 7
                        1 -1 -1];  % octant 8
                    
                    non_empty_indcs = find(obj(n_bot).voxel_matrl(:) < 5); % array of non-empty material
                    p = []; % unique position of masses
                    spring_connect_indcs = []; % unique pair index of p for constructing the spring
                    spring_type = [];
                    L_0 = [];
                    acts = [];
                    K = [];
                    
                    for i = 1:length(non_empty_indcs)
                        vox_center = obj(n_bot).voxel_pos(non_empty_indcs(i),:);
                        current_matrl = matrl(non_empty_indcs(i));
                        p_candidate = zeros(8, 3);
                        for j = 1:size(corner_signs,1)
                            p_candidate(j, :) = vox_center + corner_signs(j, :).*cube_corner;
                            % check if the position already exists
                            % putting isempty as the first expression will
                            % prevent the error in ismember()
                            if isempty(p) || ~ismember(p_candidate(j, :), p, 'rows')
                                % store in the p
                                p = [p;  p_candidate(j, :)]; %#ok<AGROW>
                            end
                        end
                        % all canndidates are used to construct the springs
                        for j = 1: size(p_candidate)
                            % get the index of the postiion in p
                            [~, node_indx] = ismember(p_candidate(j, :), p, 'rows');
                            
                            % get all posibility of neighbor positions
                            % within the same voxel at each of the eight pos candidates
                            neighber_rel_pos = -diag(2*corner_signs(j, :).*cube_corner);
                            
                            for k = 1: size(neighber_rel_pos, 1)
                                neighbor_pos = p_candidate(j, :) + neighber_rel_pos(k, :);
                                
                                % find the index of the neighbor pos in p
                                [~, neighbor_indx] = ismember(neighbor_pos, p, 'rows');
                                
                                % pair up the node and the current neighbor
                                % indx and sort in ascending order
                                candidate_connect_indcs = sort([node_indx, neighbor_indx]);
                                
                                % check if the connecting indcs already exist
                                if isempty(spring_connect_indcs) || ~ismember(candidate_connect_indcs, spring_connect_indcs, 'rows')
                                    % if it's a new pair, add it to the
                                    % list spring_connect_indcs
                                    spring_connect_indcs = [spring_connect_indcs; candidate_connect_indcs];      %#ok<AGROW>
                                    new_L0 = vecnorm(p(candidate_connect_indcs(1), :) - p(candidate_connect_indcs(2), :));
                                    
                                    % current_matrl (propeties of spring is based on the rolling basis of the voxel pos)
                                    L_0 = [L_0; new_L0]; %#ok<AGROW>
                                    switch(current_matrl)
                                        case 1 % hard
                                            b_breathing = 0;
                                            phase = 0;
                                            k_spring = 1000;
                                        case 2 % soft
                                            b_breathing = 0;
                                            phase = 0;
                                            k_spring = 200;
                                        case 3 % sine med
                                            b_breathing = 1;
                                            phase = 0;
                                            k_spring = 600;
                                        case 4 % cosine med
                                            b_breathing = 1;
                                            phase = pi/2;
                                            k_spring = 600;
                                    end
                                    acts = [acts; b_breathing*[new_L0/2, obj(n_bot).omega, phase]]; %#ok<AGROW>
                                    K = [K; k_spring]; %#ok<AGROW>
                                    spring_type = [spring_type; current_matrl]; %#ok<AGROW>
                                end
                            end
                        end
                        
                        
                    end
                    
                    % change the position and orientation of the robot after
                    % constructing the springs
                    % R = obj.rotationAxisAngle([1 0 0], pi/6); % tile around x axis by 30 degree
                    if ~isempty(p)
                        R = eye(3);
                        p = R*p'; % tilt all masses
                        p = p' + p_init_offset(n_bot,:); % add the offset (each row of 2d array);
                    end
                    obj(n_bot).masses = PointMass(repmat(obj(n_bot).mass, size(p, 1), 1), p, repmat(v_init(n_bot ,:), size(p,1), 1), repmat(a_init(n_bot, :), size(p,1), 1));
                    obj(n_bot).springs = Spring(L_0, K, spring_connect_indcs, acts, spring_type);
                    
                end
                
            end
        end
        
        %% setters
        function obj = updateP(obj, p)
            %UPDATEP Updates the position of all masses
            %   p is an array of position vectors (each row is a mass, the
            %   columns are x,y,z)
            P = num2cell(p,2);
            [obj.masses.p] = P{:};
        end
        
        function obj = updateV(obj, v)
            %UPDATEP Updates the velocity of all masses
            %   v is an array of velocity vectors (each row is a mass, the
            %   columns are x,y,z)
            V = num2cell(v,2);
            [obj.masses.v] = V{:};
        end
        
        function obj = updateA(obj, a)
            %UPDATEP Updates the acceleration of all masses
            %   a is an array of acceleration vectors (each row is a mass, the
            %   columns are x,y,z)
            A = num2cell(a,2);
            [obj.masses.a] = A{:};
        end
 
        %% energy functions
        function pe = calcPE(obj, g, t)
            my_masses = obj.masses;
            my_springs = obj.springs;
            pe = 0;
            % gravitational potential energy
            for i = 1:length(my_masses)
                pe = pe + my_masses(i).mass * abs(g(3)) * my_masses(i).p(3);
            end
            % spring energy
            for i = 1:length(my_springs)
                pair_indcs = my_springs(i).m;
                vector = my_masses(pair_indcs(1)).p - my_masses(pair_indcs(2)).p;
                L = vecnorm(vector);
                act = my_springs(i).act;
                L_act = my_springs(i).L_0 + act(1)*sin(act(2)*t + act(3));
                pe = pe + 0.5*my_springs(i).k*(L - L_act).^2;
            end
            % GRF energy will be summed in the simulation object
        end
        
        function ke = calcKE(obj)
            my_masses = obj.masses;
            ke = 0;
            for i = 1:length(my_masses)
                ke = ke + 0.5*my_masses(i).mass*(my_masses(i).v(1)^2 + ...
                    my_masses(i).v(2)^2 + my_masses(i).v(3)^2);
            end
        end
        
        function com_pos = calcCOM(obj)
            mass_pos = reshape([obj.masses.p], 3, []);
            com_pos = mean(mass_pos, 2);
        end
        
        
        %% visualization functions
        
        function plotPDF(obj)
            figure;
            [X, Y, Z] = obj.getMeshGrid();
            lim = zeros(size(obj.chromosome, 1), 2);
            
            labels = {'Hard','Soft','Sine','Cosine','Empty'};
            
            for i = size(obj.chromosome, 1): -1 : 1
                ax(i) = subplot(size(obj.chromosome, 1), 1, i);
                pcolor3(X, Y, Z, reshape(obj.pd(i, :), ...
                    obj.voxel_dim(1), obj.voxel_dim(2), ...
                    obj.voxel_dim(3)));
                colorbar;
                lim(i, :) = caxis;
                title(labels{i})
            end
            
            min_lim = min(min(lim));
            max_lim = max(max(lim));
            
            for i = 1:size(obj.chromosome, 1)
                caxis(ax(i),[min_lim max_lim])
            end
        end
        
        function [X, Y, Z] = getMeshGrid(obj)
            x_range = min(obj.voxel_pos(:,1)): obj.cube_length: max(obj.voxel_pos(:,1));
            y_range = min(obj.voxel_pos(:,2)): obj.cube_length: max(obj.voxel_pos(:,2));
            z_range = min(obj.voxel_pos(:,3)): obj.cube_length: max(obj.voxel_pos(:,3));
            [X, Y, Z] = meshgrid(x_range, y_range, z_range);
        end
        
        function plotMaterial(obj)
            %PlotMaterial Summary of this method goes here
            % visualize the output function
            matrl = obj.voxel_matrl(:);
            matrl_color = flipud(jet(4));
            
            figure;
            %             scatter3(obj.voxel_pos(:,1), obj.voxel_pos(:,2), obj.voxel_pos(:,3), [], matrl, 'x'); hold on;
            %             colormap(jet(4)); caxis([1 4]);
            
            cube_dim = obj.cube_length*ones(1,3);
            
            % create offset to see the scatter plots at the center of voxels
            vox_plot_pos = obj.voxel_pos - cube_dim/2;
            for i = 1:length(obj.voxel_pos)
                if matrl(i) ~= 5
                    voxel(vox_plot_pos(i, :), cube_dim, matrl_color(matrl(i), :), 0.25);
                end
            end
            
            colormap(jet(4));
            labels = {'Cosine','Sine','Soft','Hard'};
            lcolorbar(labels,'fontweight','bold');
            xlabel('X'); ylabel('Y'); zlabel('Z');
            
            view(3); axis equal;
        end
        
        %% force functions
        function forces = calcForces(obj, g, f_ext, t)
            %CALCFORCES Calculates the vector forces on each mass
            %   g is the gravitational constant (1x3 vector)
            %   f are additional forces on the nodes (num_masses x 3 array)
            my_masses = obj.masses;
            my_springs = obj.springs;
            forces = zeros(length(my_masses), 3);
            for i = 1:length(my_masses)
                % add gravitational force
                forces(i,:) = my_masses(i).mass * g + f_ext(i,:);
                % go through all springs in the robot
                for j = 1:length(my_springs)
                    % check if the current mass index is attahced to the current
                    % spring
                    if ismember(i, my_springs(j).m)
                        
                        % find the current spring length L
                        pair_indcs = my_springs(j).m;
                        vector = my_masses(pair_indcs(1)).p - my_masses(pair_indcs(2)).p;
                        L = vecnorm(vector);
                        act = my_springs(j).act;
                        L_act = my_springs(j).L_0 + act(1)*sin(act(2)*t + act(3));
                        spring_f = my_springs(j).k*(L - L_act);
                        
                        % create the force vector with correct direction
                        if my_springs(j).m(1) == i
                            spring_v = -spring_f*vector/L;
                        elseif my_springs(j).m(2) == i
                            spring_v = spring_f*vector/L;
                        end
                        
                        % add more forces to the mass
                        forces(i,:) = forces(i,:) + spring_v;
                    end
                end
            end
        end
        
        
        %% kinematics functions
        function [a, v, p] = calcKin(obj, f, dt, mu_s, mu_k)
            %CALCKIN Calculates the kinematics of the robot (acceleration,
            %velocity, position) based of the forces
            %   f is the array of force vectors on all the masses (not just
            %   external forces)
            %   a is the array of acceleration vectors for all masses
            %   v is the array of velocity vectors for all masses
            %   p is the array of position vectors for all masses
            my_masses = obj.masses;
            mass_pos = reshape([obj.masses.p], 3, []);
            mass_pos_z = mass_pos(3, :);
            
            a = zeros(size(f,1), size(f,2));
            v = a;
            p = v;
            
            % add friction
            fixed_indcs = [];
            slide_indc = [];
            if ~isempty(find(mass_pos_z < 0, 1))
                contact_indcs = find(mass_pos_z < 0);
                % calculate the magnitude of the horizontal forces
                f_h = vecnorm(f(contact_indcs, 1:2), 2, 2);
                f_z = abs(f(contact_indcs, 3));
                fixed_indcs = find(f_h <= f_z*mu_s);
                slide_indc = find(f_h > f_z*mu_s);
            end
            
            for i = 1:length(my_masses)
                if ismember(i, fixed_indcs)
                    % not position update
                    a(i,:) = 0;
                    v(i,:) = 0;
                    p(i,:) = my_masses(i).p;
                elseif ismember(i, slide_indc)
                    % add force from kinetic frictoin (opposite to horizontal v component)
                    f_z = abs(f(i, 3));
                    a(i,:) = (f(i,:) - mu_k*f_z*([my_masses(i).v(1:2), 0]))/ my_masses(i).mass;
                    v(i,:) = my_masses(i).v + a(i,:)*dt;
                    p(i,:) = my_masses(i).p + v(i,:)*dt;
                else % no contact with the ground
                    a(i,:) = f(i,:) / my_masses(i).mass;
                    v(i,:) = my_masses(i).v + a(i,:)*dt;
                    p(i,:) = my_masses(i).p + v(i,:)*dt;
                end
            end
        end
        
        
        
        function S = skew(~, v)
            % vec: 3 x 1 vector column
            S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
        end
        
        function R = rotationAxisAngle(obj, axis, angle)
            % axis:  double 3 x 1 vector
            % angle: double
            axis = [axis(1) axis(2), axis(3)]';
            axis = axis/vecnorm(axis);
            c = cos(angle); s = sin(angle); v = 1 - c;
            R = eye(3)*c + obj.skew(axis)*s + axis*axis'*v;
        end
        
    end
    
end

