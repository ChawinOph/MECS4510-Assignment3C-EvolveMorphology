classdef Simulator < handle
    %SIMULATOR Summary of this class goes here
    %   Detailed explanation goes here
    properties (Constant)
        g = [0, 0, -9.81]; % Double array. Gravitational acceleration (m/s^2)
        rho = 0.999; % Double. Global velocity damping parameter (0<p<=1)
        k_ground = 20000; % contact force constant (2500 default)
        mu_s = 1; % static friction coefficient (0.25 default)
        mu_k = 0.8; % kinetic friction coefficient (0.1 default)       
        run_time = 5; % seconds
    end
    
    properties
        bots % Array of robots Represents the bodies in the system
        t = 0 % Double. current time
        dt = 0.005; 
    end
    
    methods
        function obj = Simulator(bots, dt)
            %SIMULATOR Construct an instance of simulator class
            %   Detailed explanation goes here
            if nargin>0
                obj.bots = bots;
                if nargin>1
                    obj.dt = dt;
                    if nargin >2
                        disp('Too many arguments for simulator');
                    end
                end
            end            
        end
        
        function fitnesses = evaluate(obj, bots)
            % create array of robot objects the same number as no. of gene
            % columns (you have a choice to do one robot at a time by inputting just one gene)
            [~, ~, ~, fitnesses] = obj.simulate(bots);
        end
        
        function [K, V, COM, fitnesses] = simulate(obj, bots)
            % SIMULATE run the simulation without the plot'
            
            obj.bots = bots;
            
            T = 0: obj.dt : obj.run_time;
            
            K = zeros(length(T), length(obj.bots)); % kinetic energy
            V = zeros(length(T), length(obj.bots)); % potential energy
            COM = zeros(length(T), 3, length(obj.bots)); % potential energy
            
            % reset timer
            obj.t = 0;
            
            for i = 1:length(T)
                [ke, pe, com_pos] = obj.step();      
                V(i, :) = pe;
                K(i, :) = ke;
                COM(i, :, :) = reshape(com_pos, 1, 3, []);            
            end
            
            fitnesses = reshape(vecnorm(COM(end , 1:2, :) - COM(1, 1:2, :), 2, 2), 1, []);
            
        end
        
        function [frames, K, V, COM, fitness] = simulate_and_plot(obj, bots)
            % SIMULATE run the simulation over a period of 'time' in
            % seconds and record the animation for playback in 'frames'
            
            obj.bots = bots;
            
            fig = figure('pos',[10 10 900 600]);
            % crate slider for playback
            
            T = 0: obj.dt : obj.run_time;
            
            K = zeros(length(T), length(obj.bots)); % kinetic energy
            V = zeros(length(T), length(obj.bots)); % potential energy
            COM = zeros(length(T), 3, length(obj.bots)); % potential energy
            
            k = 0; % frame counter
            
            % reset timer
            obj.t = 0;
            
            for i = 1:length(T)
               
                t_step = T(i);     
                [ke, pe, com_pos] = obj.step();
                    
                V(i, :) = pe;
                K(i, :) = ke;
                COM(i, :, :) = reshape(com_pos, 1, 3, []);
                
                % desired frame rate is 25 frame/s meaning we need one
                % frame every other 0.04 sec (every other 0.04/dt frames)
                if mod(t_step, 0.04) == 0
                    k = k + 1;
                    clf;
                    obj.drawRobots();
                    
                    % draw trajectory of COM from the start
                    com_trajs_x = reshape(COM(1:i, 1, :), [], 1, 1);
                    com_trajs_y = reshape(COM(1:i, 2, :), [], 1, 1);
                    com_trajs_z = reshape(COM(1:i, 3, :), [], 1, 1);
                    scatter3(com_trajs_x, com_trajs_y, com_trajs_z, ...
                        'Marker', '.', ...
                        'MarkerEdgeColor', 'g', ...
                        'MarkerFaceColor', 'g');
                    
                    text(1, 1, 0.4, ['t = ' num2str(t_step) ' sec']);
                    
                    drawnow
                    frames(k) = getframe(fig);  %#ok<AGROW>
                end
                   
            end         
            
            fitness = reshape(vecnorm(COM(end , 1:2, :) - COM(1, 1:2, :), 2, 2), 1, []);
        end
        
        function [ke, pe, com] = step(obj)
            % loop through all robots
            ke = zeros(1,length(obj.bots));
            pe = zeros(1,length(obj.bots));
            com = zeros(3,length(obj.bots));
            
            for bot_no = 1:length(obj.bots)
                if isempty([obj.bots(bot_no).masses])
                    % if robot is empty, do nothing
                else
                    % calculate contact forces based on mass positions
                    f_contact = zeros(length(obj.bots(bot_no).masses), 3);
                    mass_pos = reshape([obj.bots(bot_no).masses.p], 3, []);
                    mass_pos_z = mass_pos(3, :);
                    % check if there are any masses underneath the ground
                    if ~isempty(find(mass_pos_z < 0, 1))
                        contact_inds = find(mass_pos_z < 0);
                        % calculate the restoration force
                        f_contact(contact_inds, 3) = -obj.k_ground*mass_pos_z(contact_inds);
                        pe_contact = 1/2*obj.k_ground*sum(abs(mass_pos_z(contact_inds).^2));
                    else
                        pe_contact = 0;
                    end
                    
                    f_ext = f_contact;
                                       
                    forces = obj.bots(bot_no).calcForces(obj.g, f_ext, obj.t);
                    
                    [a, v, p] = obj.bots(bot_no).calcKin(forces, obj.dt, obj.mu_s, obj.mu_k);
                    
                    obj.bots(bot_no).updateP(p);
                    obj.bots(bot_no).updateV(obj.rho*v);
                    obj.bots(bot_no).updateA(a);
                    
                    % get energy
                    % ke(:,bot_no) = obj.bots(bot_no).calcKE();
                    % pe(:,bot_no) = obj.bots(bot_no).calcPE(obj.g, obj.t) + pe_contact;
                    
                    com(:,bot_no) = obj.bots(bot_no).calcCOM();
                end
            end
            % update time
            obj.t = obj.t + obj.dt;
        end
        
        %% VISUALIZATION
        function drawRobots(obj)
            
            spring_color = flipud(jet(4));
            % loop through all robots
            for bot_no = 1:length(obj.bots)
                if ~isempty([obj.bots(bot_no).masses.p])
                    % get position of all point masses
                    mass_pos = reshape([obj.bots(bot_no).masses.p], 3, []);
                    scat = scatter3(mass_pos(1, :), mass_pos(2, :), mass_pos(3, :));
                    hold on;
                    scat.MarkerEdgeColor = 'k';
                    scat.MarkerFaceColor = 'k';                   
                    % obj.drawOctaSurface(bot_no)                    
                    % draw springs based on given pairs of mass indices
                    pair_indcs = reshape([obj.bots(bot_no).springs.m], 2, [])';
                    spring_types = [obj.bots(bot_no).springs.type];
                    for i = 1:size(pair_indcs, 1)                                               
                        pair_pos = reshape([obj.bots(bot_no).masses(pair_indcs(i, :)).p], 3, []);
                        plot3(pair_pos(1, :), pair_pos(2, :), pair_pos(3, :), 'Color', spring_color(spring_types(i), :), 'LineWidth', 1.5); hold on;
                    end
                end
            end
                       
            axis equal;  grid on;
%             view(-50, 25);
            view(3);
            xlim(1.5*[-1 1]);
            ylim(1.5*[-1 1]);
            zlim([-0.02 0.5]);
            xlm = xlim();
            ylm = ylim();
            % floor
            [X,Y] = meshgrid(xlm(1):0.1:xlm(2), ylm(1):0.1:(ylm(2)));
            Z = zeros(size(X));
            h = surf(X,Y,Z);
            h.FaceColor = 0.8*[1 1 1];
            h.FaceLighting = 'gouraud';
            
            % add light 
            lightangle(-45,30)
            
            xlabel('x (m)')
            ylabel('y (m)')
            zlabel('z (m)')
        end
        
        function drawSurface(~, vertex_pos)
            h2 = fill3(vertex_pos(1,:), vertex_pos(2,:), vertex_pos(3,:), 'g');
%             h2.FaceAlpha = 0.25;
            h2.EdgeColor = 'none';
            h2.FaceLighting = 'gouraud';
        end
        
         function drawOctaSurface(obj, bot_no)
            mass_pos = reshape([obj.bots(bot_no).masses.p], 3, []);
            % top surface
            obj.drawSurface(mass_pos(:, [1 2 3])); hold on;
            obj.drawSurface(mass_pos(:, [1 3 4]));
            obj.drawSurface(mass_pos(:, [1 4 5]));
            obj.drawSurface(mass_pos(:, [1 5 2]));
            obj.drawSurface(mass_pos(:, [6 2 3]));
            obj.drawSurface(mass_pos(:, [6 3 4]));
            obj.drawSurface(mass_pos(:, [6 4 5]));
            obj.drawSurface(mass_pos(:, [6 5 2]));
        end
        
    end
end

