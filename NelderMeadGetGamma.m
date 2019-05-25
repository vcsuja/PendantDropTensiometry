function [ST,predContours] = NelderMeadGetGamma(DropContour,DeltaRho,g,MaxReload)
            %Inputs:   DropContour - Contour of the drop (without the
            %          capillary) in scaled to meters.
            %          DeltaRho - the density difference between the suspended
            %          fluid and the ambient fluid.  
            %          g - accelaration due to gravity
            %          MaxReload - *optional* Number of times the
            %          algorithm reloads to reach the global minimum 
            %Outputs:  ST - Surface tension
            %          predContours - Numerically computed contours used 
            %          for calculating the surface tension.
            
            % Details: 2009, N. Alvarez - A non-gradient based algorithm for the 
            % determination of surface tension from a pendant drop: Application 
            % to low Bond number drop shapes
            
            % Author: Vinny (vineethcs.cet@gmail.com) (Based on orginal code by
            % J. Frostad (john.frostad@ubc.ca))
            
            % Handle optional arguments
            if nargin<4
                MaxReload = 15; %Max # of times to reload fitting algorithm
            end
            if nargin<3
                g = 9.81;
            end

            
            %Extract relavent geometrical points from DropContour
                %(x1,y1) = left point at location of maximum drop width
                %(x2,y2) = right point at location of maximum drop width
                %(x3,y3) = location of drop apex
                %(x4,y4) = left point at max drop width above apex
                %(x5,y5) = right point at max drop width above apex
                xData = DropContour(:,1);
                yData = DropContour(:,2);
                
                x1 = min(xData); x1Index = find(xData==x1); y1 = mean(yData(x1Index));
                x2 = max(xData); x2Index = find(xData==x2); y2 = mean(yData(x2Index));
                y3 = max(yData); y3Index = find(yData==y3); x3 = mean(xData(y3Index));
                D =(x2-x1);                
                y4Index = find(yData>=y3-D,1); x4 =  xData(y4Index);
                y5Index = find(yData>=y3-D,1,'last'); x5 =  xData(y5Index);
                S = mean(x5-x4)/D;
                
            %Make a guess for the fminsearch function optimizing
            %theoretical and predicted drop fits:
                STguess = AHTFitGetGamma(S,D,DeltaRho,g); % ST guess from AHT    
                Ro = (x2 - x1)/2;

            %Formulate Initial Guess for fitting parameters
            Guess = [STguess x3 y3 Ro 0];
            
                %Neldor-Mead Reloading Loop
                for j = 1:MaxReload+1
                    tic
                    [P,~] = fminsearch(@(P) FitLowBondShape(P),Guess);
                    ST = P(1); 
                    X0 = P(2);
                    Y0 = P(3);
                    R0 = P(4);
                    Theta = P(5);
                    Difference = max(abs(Guess-P));
                    Guess = P;
                    
                    toc
                    if Difference < 1e-5 
                        disp(['Image '...
                            ' analyzed after ' num2str(j-1) ' reloads.'])
                        break; 
                    end
                end
                
                
        % Return the fitted points after truncating and rotating,shifting and scaling them back 
        % to the original form
        
            % Truncate
            Cutoff = find(yTheory<=max(yProcessed));
            xTheory = xTheory(Cutoff,1);
            yTheory = yTheory(Cutoff,1);
            % Rotate
            predContours = [R0*xTheory+X0,-R0*yTheory+Y0]*[cos(Theta), -sin(Theta); sin(Theta),cos(Theta)]; %Inverse rotation  
                
        %% FitLowBondShape nested function
        
        function Error = FitLowBondShape(P)
            
            Gamma = P(1); %Interfacial Tension in N/m
            Xo = P(2); %Shift in x coord of drop apex in meters
            Yo = P(3); %Shift in z coord of drop apex in meters
            Ro = P(4); %Drop Radius at apex in meters
            Theta = P(5); %Camera rotation offset angle
            Bo = DeltaRho*g*Ro^2/Gamma; %Bond 

            %Shift and rotate coordinates
            xProcessed = (xData - Xo)*cos(Theta) + (yData - Yo)*sin(Theta);
            yProcessed = (yData - Yo)*cos(Theta) - (xData - Xo)*sin(Theta);
            %Nondimensionalize data for comparing to theoretical data
            xProcessed = xProcessed/Ro; yProcessed = -yProcessed/Ro; %Apex of experimental data is max(yData)
            
            %Refine # of Theoretical boundary points to be ~ twice the number of data
            N = length(xProcessed);
            
            %Solve the Young-Laplace equation for the given guess of parameters 
            system = @(s,y) ShapeODE(s,y,Bo);
            Span = linspace(0,1.5*pi,N);
            [~,y] = ode45(system,Span,[0 0 0],[]);

            %Generate full drop shape from theoretical data
            xTheory = [-flipud(y(:,1));y(:,1)];
            yTheory = [flipud(y(:,2));y(:,2)];
    

            %Calculate the Error in the fit by RMS difference between data and closest
            %theoretical value
            E = zeros(N,1); %Point by point RMS error
            Search = round(.06*2*N); % # of nearby points to check for closest theoretical point
            CloseIndex = E; %Index of closest points on theoretical curve
            
            for i = 1:N
                %Search for closest point to data on theoretical boundary.
                %Search entire curve on first point.
                %Search only nearby points on subsequent data points.
                %Limit search indices to values within bounds of theoretical points.
                if i == 1 
                    Range = 1:length(xTheory);
                elseif Index <= Search+1
                    Range = 1:Index+Search;
                elseif Index >= length(xTheory)-Search
                    Range = Index-Search:length(xTheory);
                else
                    Range = Index-Search:Index+Search;
                end
                
                %Calculate RMS distance and store closest distance and index of point
                Distance = sqrt((xTheory(Range) - xProcessed(i)).^2 + ...
                    (yTheory(Range) - yProcessed(i)).^2);
                [E(i),TempIndex] = min(Distance);
                Index = Range(1)-1+TempIndex;
                CloseIndex(i) = Index;
            end

            %Sum all of the RMS distances as the total error
            Error = sum(E);
            end
                         
            
              
end


function ydot = ShapeODE(~,y,Bo)
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