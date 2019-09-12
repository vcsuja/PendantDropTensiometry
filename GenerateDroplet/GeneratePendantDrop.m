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
        
        function [dropImages, dropLabels] = generateDropImages(obj,imgSize,imgAugment)
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
              dropLabels = zeros(1,length(obj.Bo)*(imgAugment+1));
              
            % Initalize progress bar:
            f = waitbar(0,'1','Name','Generating Drop Profiles',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
        
            imgCount = 1;
            for i =1:length(obj.Bo)
                dropImages(imgCount,:,:)=obj.profile2Image(obj.dropProfiles{1,i},obj.dropProfiles{2,i},imgSize,0,0);
                imgCount = imgCount+1;
                
                %Randomly translate the drop:
                 for j=1: floor(imgAugment/3)+mod(imgAugment,3)
                     dropImages(imgCount,:,:)=obj.profile2Image(obj.dropProfiles{1,i},obj.dropProfiles{2,i},imgSize,1,0);
                     imgCount = imgCount+1;
                 end
                %Randomly rotate the drop:
                 for j=1: floor(imgAugment/3)
                     dropImages(imgCount,:,:)=obj.profile2Image(obj.dropProfiles{1,i},obj.dropProfiles{2,i},imgSize,0,1);
                     imgCount = imgCount+1;
                 end
                %Randomly do both:
                 for j=1: floor(imgAugment/3)
                     dropImages(imgCount,:,:)=obj.profile2Image(obj.dropProfiles{1,i},obj.dropProfiles{2,i},imgSize,1,1);
                     imgCount = imgCount+1;
                 end
                 
                %Display progress:
                % Check for clicked Cancel button
                if getappdata(f,'canceling')
                    break;
                end

                % Update waitbar and message
                waitbar(i/length(obj.Bo),f,'Generating Profiles');      
                % Update labels:
                dropLabels(1,(i-1)*(imgAugment+1)+1:i*(imgAugment+1)) = obj.Bo(i);
                
            end
              
            delete(f) %Close the progress bar. 
            
            %Shuffle labels and dropImages:
            randIndx = randperm(length(dropLabels));
            dropLabels = dropLabels(randIndx);
            dropImages = dropImages(randIndx,:,:);
            
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
        
        function im = profile2Image(~,xProfile,yProfile,nsize,translate,rotate)
        % This function generates an image corresponding to the supplied
        % profile
        
        figure('Units', 'pixels','Position',[500,500,300,300],'Color','white','visible','off')        
        ylim([-max(yProfile)*0.25 1.25*max(yProfile)]);
        hold on
        xlim(1.25*[-max(xProfile) max(xProfile)]);
        
        %Translate:
        if translate
            xtrans = rand()*max(xProfile)*0.05*(-1)*(rand()>0.5);
            ytrans = rand()*max(yProfile)*0.15*(-1)*(rand()>0.5);     
        else
            xtrans = 0;
            ytrans = 0;
        end
        %Rotae
        if rotate
            theta = rand()*5*(-1)*(rand()>0.5);    
        else
            theta = 0;
        end  
        
        %Apply the affine translation
        xProfile = (xProfile+xtrans)*cosd(theta) + (yProfile+ytrans)*sind(theta);
        yProfile = -(xProfile+xtrans)*sind(theta) + (yProfile+ytrans)*cosd(theta); 
        plot(xProfile+xtrans,yProfile+ytrans,'Linewidth',1.4,'Color',[0,0,1]);
        
        %Axis controls:
        axis equal          
        set(gca, 'Units', 'pixels', 'Position', [10, 10, nsize, nsize]);
        
        fr = getframe(gca,[0,0,nsize,nsize]);
        close(gcf);
        im = frame2im(fr);
        im = im2double(:,:,2));           
        
        %%%%Depracted due to an unresolved error%%
        %im = im2double(im(round(end-nsize)/2:round(end-nsize)/2+nsize-1,round(end-nsize)/2:round(end-nsize)/2+nsize-1,2));       
        %%%%Depracted due to an unresolved error%%
        end
    end
    
    
end