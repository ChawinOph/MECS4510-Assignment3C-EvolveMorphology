classdef Simulator < handle
    %SIMULATOR Summary of this class goes here
    %   Detailed explanation goes here
    properties (Constant)
        g = [0, 0, -9.81]; % Double array. Gravitational acceleration (m/s^2)
        rho = 1; % Double. Global velocity damping parameter (0<p<=1)
        k_ground = 2500; % contact force constant (2500 default)
        mu_s = 0.5; % static friction coefficient (0.25 default)
        mu_k = 0.25; % kinetic friction coefficient (0.1 default)       
        run_time = 3; % seconds
    end
    
    properties
        bots % Array of robots Represents the bodies in the system
        t = 0 % Double. current time
        dt = 0.005; 
    end
    
    methods
        function obj = simulator(bots, dt)
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
            
            fitnesses = reshape(vecnorm(COM(end , 1:2, :), 2, 2), 1, []);
            
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
                    text(0.4, 0.4, 0.4, ['t = ' num2str(t_step) ' sec']);
                    drawnow
                    frames(k) = getframe(fig);  %#ok<AGROW>
                end
                   
            end         
            
            fitness = reshape(vecnorm(COM(end , 1:2, :), 2, 2), 1, []);
        end
        
        function [ke, pe, com] = step(obj)
            % loop through all robots
            ke = zeros(1,length(obj.bots));
            pe = zeros(1,length(obj.bots));
            com = zeros(3,length(obj.bots));
            
            for bot_no = 1:length(obj.bots)
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
               
                % calculate kinematic variables
%                 tic
%                 disp('calcForce')
                forces = obj.bots(bot_no).calcForces(obj.g, f_ext, obj.t);
%                 toc
                % check the friction on the masses in contact
%                 tic
%                 disp('calcKin')
                [a, v, p] = obj.bots(bot_no).calcKin(forces, obj.dt, obj.mu_s, obj.mu_k);
%                 toc
                % update all kinematic variables
%                 tic
                obj.bots(bot_no).updateP(p);
                obj.bots(bot_no).updateV(obj.rho*v);
                obj.bots(bot_no).updateA(a);        
%                 toc
                % get energy
                ke(:,bot_no) = obj.bots(bot_no).calcKE();
                pe(:,bot_no) = obj.bots(bot_no).calcPE(obj.g, obj.t) + pe_contact;
                com(:,bot_no) = obj.bots(bot_no).calcCOM();                        
            end
            % update time
            obj.t = obj.t + obj.dt;
        end
        
        %% VISUALIZATION
        function drawRobots(obj)
              
            % loop through all robots
            for bot_no = 1:length(obj.bots)
                % get position of all point masses
%                 mass_pos = reshape([obj.bots(bot_no).masses.p], 3, []);
%                 scat = scatter3(mass_pos(1, :), mass_pos(2, :), mass_pos(3, :));
%                 hold on;
%                 scat.MarkerEdgeColor = 'k';
%                 scat.MarkerFaceColor = 'b';

                obj.drawOctaSurface(bot_no)

%                 draw springs based on given pairs of mass indices
                pair_indcs = reshape([obj.bots(bot_no).springs.m], 2, [])';
                for i = 1:size(pair_indcs, 1)
                    pair_pos = reshape([obj.bots(bot_no).masses(pair_indcs(i, :)).p], 3, []);
                    plot3(pair_pos(1, :), pair_pos(2, :), pair_pos(3, :), 'k'); hold on;
                end
            end
                       
            axis equal;  grid on;
            view(-50, 25)
            xlim([-1 1]);
            ylim([-1 1]);
            zlim([-0.02 0.5]);
            xlm = xlim();
            ylm = ylim();
            % floor
            [X,Y] = meshgrid(xlm(1):0.1:xlm(2), ylm(1):0.1:(ylm(2)));
            Z = zeros(size(X));
            h = surf(X,Y,Z);
            h.FaceColor = 0.9*[1 1 1];
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
        
        function drawAllStarfishSurface(obj, bot_no)
            mass_pos = reshape([obj.bots(bot_no).masses.p], 3, []);
            % top surface
            obj.drawSurface(mass_pos(:, [14 15 16 17])); hold on;
            obj.drawSurface(mass_pos(:, [14 15 7]));
            obj.drawSurface(mass_pos(:, [16 15 9]));
            obj.drawSurface(mass_pos(:, [16 17 11]));
            obj.drawSurface(mass_pos(:, [14 17 13]));
            
            % east leg
            obj.drawSurface(mass_pos(:, [1 5 7])); 
            obj.drawSurface(mass_pos(:, [1 6 7])); 
            obj.drawSurface(mass_pos(:, [1 5 13])); 
            obj.drawSurface(mass_pos(:, [1 6 13])); 
            obj.drawSurface(mass_pos(:, [14 6 13])); 
            obj.drawSurface(mass_pos(:, [14 6 7])); 
            
            % north leg
            obj.drawSurface(mass_pos(:, [2 5 7])); 
            obj.drawSurface(mass_pos(:, [2 7 8])); 
            obj.drawSurface(mass_pos(:, [2 8 9])); 
            obj.drawSurface(mass_pos(:, [2 5 9])); 
            obj.drawSurface(mass_pos(:, [15 7 8])); 
            obj.drawSurface(mass_pos(:, [15 8 9])); 
            
            % west leg
            obj.drawSurface(mass_pos(:, [3 5 9])); 
            obj.drawSurface(mass_pos(:, [3 9 10])); 
            obj.drawSurface(mass_pos(:, [3 10 11])); 
            obj.drawSurface(mass_pos(:, [3 5 11])); 
            obj.drawSurface(mass_pos(:, [16 9 10])); 
            obj.drawSurface(mass_pos(:, [16 10 11])); 
            
            % south leg
            obj.drawSurface(mass_pos(:, [4 5 11])); 
            obj.drawSurface(mass_pos(:, [4 11 12])); 
            obj.drawSurface(mass_pos(:, [4 12 13])); 
            obj.drawSurface(mass_pos(:, [4 5 13])); 
            obj.drawSurface(mass_pos(:, [17 11 12])); 
            obj.drawSurface(mass_pos(:, [17 12 13])); 
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

