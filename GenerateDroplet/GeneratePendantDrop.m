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
                
                %Update the properties of the class:
                obj.Bo  = Bo;
                obj.dropProfiles = dropProfiles;
                
            end
        end
        
        function dropImages = generateDropImages(obj,imgSize,imgAugment)
            % Generates dropImages for training a machine learning model
            % Inputs: imgSize -  The size of the output image (default: 64x64)
            %       : imgAugment - The number of translational and rotational augments (default: 0)
            
            % Check if dropProfiles are available
            if (isempty(obj.dropProfiles))
                error('Run generateDropProfiles() to generate the dropProfiles before running this function');
            end
            
            % Handle optional arguments:
            if (nargin<3)
                imgAugment = 0;
            end
            if (nargin<2)
                imgSize = 64;
            end
           
            % Preallocate memory for the arrays:
              dropImages = zeros(length(obj.Bo)*(imgAugment+1),imgSize,imgSize);
              
            for i =1:length(obj.Bo)
                dropImages(i,:,:)=obj.profile2Image(obj.dropProfiles{1,i},obj.dropProfiles{2,i},imgSize);
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
        
        function im = profile2Image(~,xProfile,yProfile,nsize)
        % This function generates an image corresponding to the supplied
        % profile
        
        figure('Units', 'pixels','Position',[500,500,300,300],'Color','white','visible','off')
        plot(xProfile,yProfile,'Linewidth',1.4,'Color',[0,0,1])
        ylim([-max(yProfile)*0.25 1.25*max(yProfile)]);
        xlim(1.1*[-max(xProfile) max(xProfile)]);
        axis equal
           
        set(gca, 'Units', 'pixels', 'Position', [10, 10, nsize, nsize]);
        fr = getframe(gca,[0,0,nsize,nsize]);
        close(gcf);
        im = frame2im(fr);
        im = im2double(im(round(end-nsize)/2:round(end-nsize)/2+nsize-1,round(end-nsize)/2:round(end-nsize)/2+nsize-1,2));           
        end
    end
    
    
end