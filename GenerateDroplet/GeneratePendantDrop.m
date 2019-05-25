classdef GeneratePendantDrop < handle
    properties(Access = private)
        
    end
    properties 
        CapillaryOn    % Boolean that determines if the capillary is present in the final profile
        Bo             % Vector of Bond numbers
        dropProfiles   % Cell of drop profiles
    end
    
    
    methods % Generic functions
        
        function obj = GeneratePendantDrop()
            %This function creates an instance of the class
        end                
    end    
    
    methods % Core functions
        function dropProfiles = generateDropProfiles(obj,Bo)
            %This function generates and returns the drop profile for a
            %given Bo (Bond number). Bo can be a vector
            
            %Create a cell to hold the dropProfiles
            dropProfiles = cell(2,length(Bo));
            
            for i = 1:length(Bo)
                %Solve the Young-Laplace equation for the given guess of parameters 
                system = @(s,y) obj.ShapeODE(s,y,Bo(i));
                Span = linspace(0,1.5*pi,360);
                [~,y] = ode45(system,Span,[0 0 0],[]);
                
                %Generate full drop shape from theoretical data
                xTheory = [-flipud(y(:,1));y(:,1)];
                yTheory = [flipud(y(:,2));y(:,2)];
            
                %Save profile 
                dropProfiles{1,i} =xTheory;
                dropProfiles{2,i} =yTheory;
                
            end
        end
    end
    
    
    methods (Access = private)
        %Special supporting functions
        function ydot = ShapeODE(~,~,y,Bo)
        % The Young - Laplace system of differential equations
        ydot(1) = cos(y(3));
        ydot(2) = sin(y(3));
        %Prevent Divide by Zero Issues
        if y(1) == 0
            ydot(3) = 2 - Bo*y(2) - 1;
        else
            ydot(3) = 2 - Bo*y(2) - sin(y(3))/y(1);
        end
        ydot = ydot';
        end
    end
    
    
end