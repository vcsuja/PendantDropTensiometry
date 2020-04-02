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
                imgSize = 128;
            end
            
            % Limiters:
            curvaturethresh = 1.5;
            tresh_caplow =0.99;
            ncap = 20; %Capillary positions per image
            
            
            
        
            % Preallocate memory for the arrays:
              dropImages = zeros(length(obj.Bo)*(imgAugment+1)*ncap,imgSize,imgSize);
              dropLabels = zeros(2,length(obj.Bo)*(imgAugment+1)*ncap); % First input is the 1/Bond number and the nondimensional radius of the capillary
              
              
              
            % Initalize progress bar:
            f = waitbar(0,'1','Name','Generating Drop Profiles',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
        
            imgCount = 1;
            
   
            for i =1:length(obj.Bo)
                %Extract the drop profile: 
                X = obj.dropProfiles{1,i};
                Y=obj.dropProfiles{2,i};
                
                % Find curvature and limit searchspace to attach the capillary:
                [~,R,~] = curvature([X,Y]);
                ind = find((1./R(1:round(length(R)/2)))>curvaturethresh,1,'last'); 
                if isempty(ind)
                    ind=1;
                end
                Xrange = X(ind:end-ind+1);
                Yrange = Y(ind:end-ind+1);

                %Find the starting location of the capillary:
                capr1 = abs(find(abs(Xrange(1:round(length(Xrange)/2)))>1,1,'last'));
                capleftidx = abs(find(abs(Xrange(1:capr1))<tresh_caplow,1,'last'));

                    %Index positions of the capillary
                    capidxrange = round(1:capleftidx/ncap:capleftidx) + round(rand(1,ncap)*3);
                    capidxrangeend=length(Xrange)-capidxrange+1;
                    
              for kk = 1:length(capidxrange) 
                    %Core section where the droplet profiles are generated
                    
                    %Upate base case
                    dropImages(imgCount,:,:)=obj.profile2Image(Xrange(capidxrange(kk):capidxrangeend(kk)),Yrange(capidxrange(kk):capidxrangeend(kk)),imgSize,0,0);
                    dropLabels(1,imgCount) = 1/obj.Bo(i); dropLabels(2,imgCount) = 2*abs(Xrange(capidxrange(kk))); 
                    imgCount = imgCount+1;

                    %Note: Zoom will automatically be applied in all the cases:
                    %Randomly translate the drop:
                     for j=1: floor(imgAugment/3)+mod(imgAugment,3)
                        dropImages(imgCount,:,:)=obj.profile2Image(Xrange(capidxrange(kk):capidxrangeend(kk)),Yrange(capidxrange(kk):capidxrangeend(kk)),imgSize,1,0);
                        dropLabels(1,imgCount) = 1/obj.Bo(i); dropLabels(2,imgCount) = 2*abs(Xrange(capidxrange(kk))); 
                        imgCount = imgCount+1;
                     end
                    %Randomly rotate the drop:
                     for j=1: floor(imgAugment/3)
                        dropImages(imgCount,:,:)=obj.profile2Image(Xrange(capidxrange(kk):capidxrangeend(kk)),Yrange(capidxrange(kk):capidxrangeend(kk)),imgSize,0,1);
                        dropLabels(1,imgCount) = 1/obj.Bo(i); dropLabels(2,imgCount) = 2*abs(Xrange(capidxrange(kk))); 
                        imgCount = imgCount+1;
                     end
                    %Randomly do both:
                     for j=1: floor(imgAugment/3)
                        dropImages(imgCount,:,:)=obj.profile2Image(Xrange(capidxrange(kk):capidxrangeend(kk)),Yrange(capidxrange(kk):capidxrangeend(kk)),imgSize,1,1);
                        dropLabels(1,imgCount) = 1/obj.Bo(i); dropLabels(2,imgCount) = 2*abs(Xrange(capidxrange(kk))); 
                        imgCount = imgCount+1;
                     end

                    %Display progress:
                    % Check for clicked Cancel button
                    if getappdata(f,'canceling')
                        break;
                    end
              end
                % Update waitbar and message
                waitbar(i/length(obj.Bo),f,'Generating Profiles');                  
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

    function [L,R,k] = curvature(X)
    % Radius of curvature and curvature vector for 2D or 3D curve
    %  [L,R,Kappa] = curvature(X)
    %   X:   2 or 3 column array of x, y (and possibly z) coordiates
    %   L:   Cumulative arc length
    %   R:   Radius of curvature
    %   k:   Curvature vector
      N = size(X,1);
      dims = size(X,2);
      if dims == 2
        X = [X,zeros(N,1)];  % Do all calculations in 3D
      end
      L = zeros(N,1);
      R = NaN(N,1);
      k = NaN(N,3);
      for i = 2:N-1
        [R(i),~,k(i,:)] = circumcenter(X(i,:)',X(i-1,:)',X(i+1,:)');
        L(i) = L(i-1)+norm(X(i,:)-X(i-1,:));
      end
      i = N;
      L(i) = L(i-1)+norm(X(i,:)-X(i-1,:));
      if dims == 2
        k = k(:,1:2);
      end
   
      function [R,M,k] = circumcenter(A,B,C)
    % Center and radius of the circumscribed circle for the triangle ABC
    %  A,B,C  3D coordinate vectors for the triangle corners
    %  R      Radius
    %  M      3D coordinate vector for the center
    %  k      Vector of length 1/R in the direction from A towards M
    %         (Curvature vector)
      D = cross(B-A,C-A);
      b = norm(A-C);
      c = norm(A-B);
      if nargout == 1
        a = norm(B-C);     % slightly faster if only R is required
        R = a*b*c/2/norm(D);
        return
      end
      E = cross(D,B-A);
      F = cross(D,C-A);  
      G = (b^2*E-c^2*F)/norm(D)^2/2;
      M = A + G;
      R = norm(G);  % Radius of curvature
      if R == 0
        k = G;
      else
        k = G'/R^2;   % Curvature vector
      end
    end

    end

    function rotated_image = imrotate255(image, rot_angle_degree)
        %Rotates and fills new regions in the image with white space.
        tform = affine2d([cosd(rot_angle_degree)    -sind(rot_angle_degree)     0; ...
                          sind(rot_angle_degree)     cosd(rot_angle_degree)     0; ...
                          0                          0                          1]);

          rotated_image = imwarp(image, tform, 'interp', 'cubic', 'fillvalues', 255);
     end
        
        
        
        
        
        
        
        function im = profile2Image(~,xProfile,yProfile,nsize,translate,rotate)
        % This function generates an image corresponding to the supplied
        % profile
        

        set(gca, 'Units', 'pixels');
        hold on
        axis equal
        set(gcf,'Color',[1 1 1]);
        set(gca,'Position', [10, 10, nsize*1.2, nsize*1.2],'XColor', 'none', 'YColor', 'none');
        
        
        
        
%         figure('Units', 'pixels','Position',[500,500,300,300],'Color','white','visible','off')        
%         ylim([-max(yProfile)*0.25 1.25*max(yProfile)]);
%         hold on
%         xlim(1.25*[-max(xProfile) max(xProfile)]);
%         
%         %Translate:
%         if translate
%             xtrans = rand()*max(xProfile)*0.05*(-1)*(rand()>0.5);
%             ytrans = rand()*max(yProfile)*0.15*(-1)*(rand()>0.5);     
%         else
%             xtrans = 0;
%             ytrans = 0;
%         end
%         %Rotae
%         if rotate
%             theta = rand()*5*(-1)*(rand()>0.5);    
%         else
%             theta = 0;
%         end  
%         
%         %Apply the affine translation
%         xProfile = (xProfile+xtrans)*cosd(theta) + (yProfile+ytrans)*sind(theta);
%         yProfile = -(xProfile+xtrans)*sind(theta) + (yProfile+ytrans)*cosd(theta); 
%         plot(xProfile+xtrans,yProfile+ytrans,'Linewidth',1.4,'Color',[0,0,1]);
%         
%         %Axis controls:
%         axis equal
%         set(gca, 'Units', 'pixels');
%         set(gca,'Position', [10, 10, nsize, nsize],'XColor', 'none', 'YColor', 'none');
%         fr = getframe(gca,[0,0,nsize,nsize]);
%         close(gcf);
%         im = frame2im(fr);
%         im = im2double(im(:,:,2));           
        
        %%%%Depracted due to an unresolved error%%
        %im = im2double(im(round(end-nsize)/2:round(end-nsize)/2+nsize-1,round(end-nsize)/2:round(end-nsize)/2+nsize-1,2));       
        %%%%Depracted due to an unresolved error%%
        end
    end
    
    
end