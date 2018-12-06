classdef PointMass < handle
    %point_mass This class models a point mass
    %   Detailed explanation goes here
    
    properties
        mass % Double. The mass of the point [kg]
        p % Double Array. The 3d position vector [meters]
        v % Double Array. The 3d velocity vector [meters/second]
        a % Double Array. The 3d acceleration vector [meters/second^2]
    end
    
    methods
        % Constructor
        function obj = PointMass(mass, p, v, a)
            %POINT_MASS Construct an instance of point_mass
            %   Detailed explanation goes here
            if nargin>1
                for i = length(mass):-1:1
                    obj(i).mass = mass(i);
                    obj(i).p = p(i,:);

                    switch nargin
                        case 2
                            obj(i).v = [0 0 0];
                            obj(i).a = [0 0 0];
                        case 3
                            obj(i).v = v(i,:);
                            obj(i).a = [0 0 0];
                        case 4
                            obj(i).v = v(i,:);
                            obj(i).a = a(i,:);
                    end
                end
            end
        end
    end
end

