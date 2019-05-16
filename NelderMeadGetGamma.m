function [Gamma,S,x,y,xP,yP,scale] = GetGamma(image,FullMethod,NeedleD,...
    DropFluidRho,SuspendFluidRho,NPoints)
            %Refactor and clean code to split into the simple calculation of 
            %surface tension

            %User Input Parameters
            MaxReload = 15; %Max # of times to reload fitting algorithm
            flag2 = 1; %See location of boundary points used for calculating parameters
            flag3 = 1; %Save figure of fit.
            Skip = 999; 
            Skip2 = 999;% Image numbers to Skip
            Start = 1; %Image number ot start with
            Tolerance = 0.1;
            N = 1; %# of files to process

            %Initialize Data Arrays
            D = zeros(N,1); %Maximum drop width (for simplified calculation)
            Gamma = D; %Interfacial Tension
            d = D; % Drop width at D from bottom of drop (for simplified calc)
            Error = D; %Error in Young-Laplace fit
            scale = D; %conversion factor - meters/pixel
            R = D; %Value of fitted drop radius
            Xo = D; %Value of fitted shift in x coord of drop apex
            Zo = D; %Value of fitted shift in z coord of drop apex
            Theta = D; %Value of fitted rotation angle of drop coord system
            S = D; %ratio of d/D (for simplified calc)
            ST = D; %Interfacial tension by simple calculation
            Volume = D;
            Area = D;
            Bond = D;
            DeltaRho = DropFluidRho - SuspendFluidRho; %Relative Density in kg/m^3

            % Loop for processing images to fit Young-Laplace equation
            for i = Start:N
                if i == Skip
                    ST(i)=0;
                elseif i==Skip2
                    ST(i)= 0;
                else
                A = image;
                %Use Otsu's method to chose dividing line between black and white
                Intensity = graythresh(A);  
                BW = imcomplement(im2bw(A,Intensity)); %convert to a binary image
                %Trace drop boundary using bwtraceboundary, the outer boundary of the 
                %drop should be all dark with no glares, this can be achieved by
                %adjusting the backlighting in the experiment.
                row = 10; col = find(BW(row,:), 1);
                contour = bwtraceboundary(BW,[row, col], 'S',4,inf,'counterclockwise');
                x = contour(:,2); y = contour(:,1); %x and y valies of boundary points

            %Find Relevant Boundary Points for performing Young-Laplace fit
                %(x1,y1) = left point at location of maximum drop width
                %(x2,y2) = right point at location of maximum drop width
                %(x3,y3) = location of drop apex
                %(x4,y4) = left point at max drop width above apex
                %(x5,y5) = right point at max drop width above apex
                %(x6,y6) = array of points along left side of needle
                %(x7,y7) = array of points along right side of needle
                %(x8,x9) = first drop point used for fitting
                %(x9,y9) = last drop point used for fitting
                x1 = min(x); x1Index = find(x==x1); y1 = mean(y(x1Index));
                x2 = max(x); x2Index = find(x==x2); y2 = mean(y(x2Index));
                y3 = max(y); y3Index = find(y==y3); x3 = mean(x(y3Index));
                D(i) = x2-x1; %max drop width
                y4Index = find(y>=y3-D(i),1); y4 = y(y4Index); x4 = x(y4Index);
                y5Index = find(y>=y3-D(i),1,'last'); y5 = y(y5Index); x5 = x(y5Index);
                d(i) = mean(x5-x4);%drop width at max drop width above apex
                x6 = x(1:NPoints); y6 = y(1:NPoints);
                y7Index = find(y(1:end-5)>=y6(1),NPoints,'last'); 
                x7 = x(y7Index);
                y7 = y (y7Index);
                scale(i) = NeedleD/mean(x7-x6);
                xP = [x1;x2;x3;x4;x5;x6;x7];
                yP = [y1;y2;y3;y4;y5;y6;y7];

            %Run Simplified Calculation of Surface Tension
                d(i) = d(i)*scale(i); D(i) = D(i)*scale(i);
                S(i) = d(i)/D(i);
                ST(i) = SimpleSTCalculation(S(i),D(i),DeltaRho);

            %Run Full Calculation of Young-Laplace Fit
                if FullMethod
                    %Calculate offset angle of needle from perpendicular to image edge
                    x6 = x(1:NPoints); y6 = y(1:NPoints);
                    p = polyfit(y6,x6,1); slope = p(1); intercept = p(2);
                    Angle = atan(slope);
                    y7Index = find(y(1:end-5)>=y6(1),NPoints,'last'); 
                    y7 = y(y7Index); x7 = x(y7Index);
                    %Use pixel width of needle (corrected for offset angle) and known
                    %physical width of needle to get conversion factor.
                    scale(i) = NeedleD/(mean(x7-x6)*cos(Angle));
                    %Accounting for offset angle locate first and last  points
                    %to be used when fitting to Young-Laplace equation.
                   
                    % Strategy for finding drop contour start: Find the
                    % intercept of a straight line drawn over the capillary
                    % and the drop contour:
%                     x8IndexPre = find(abs(x-(slope*y+intercept))<1e-10);
%                     x8Index = x8IndexPre(find(diff(x8IndexPre)>2,1)+1);
%                     x8 = x(x8Index); y8 = y(x8Index);
%                     D1 = (slope*y(x8Index)+intercept)-x(x8Index);
                    
                    x8Index = find(abs((x-(slope*y+intercept))./(slope*y+intercept))>= Tolerance,1);
                    y8Index = find(y==y(x8Index)+7,1); % Shift up by 5 pixels (can cast as a Tol)
                    x8 = x(y8Index); y8 = y(y8Index);
                    
                    D1 = (slope*y(x8Index)+intercept)-x(x8Index);
                    y9Index = find(y>y8,1,'last'); % original code,
                    x9Index=y9Index;
                    
                    x9 = x(y9Index); y9 = y(y9Index);
                    %but MRH and XS ran into errors b/c of this line. So the next line
                    %was used.
%                     x9Index = y5Index - abs(x8Index - y4Index);
% 
%                     x9 = x(x9Index); y9 = y(x9Index);
%                     x8 = x4;
%                     x9 = x5;
%                     y8 = y4;
%                     y9 = y5;        
% 
                    plot(x,y,'k-'); hold on; axis equal
            plot(x1,y1,'y*');
            plot(x2,y2,'b*');
            plot(x3,y3,'r*');
            
                    plot(x6,y6,'y*');
                    plot(x7,y7,'c*');
                    
                            plot(x4,y4,'r*');
                    plot(x5,y5,'b*');
                    plot(x8,y8,'g*');
                    plot(x9,y9,'m*');
                    
                    %Use bondary points to make an initial guess for fitting parameters
                    Ro = (x2 - x1)*scale(i)/2;
                    if flag2
                        xs = [x1;x2;x3;x4;x5;x6;x7;x8;x9]; 
                        ys = [y1;y2;y3;y4;y5;y6;y7;y8;y9];
                        figure, imshow(A); hold on;
                        plot(x,y,'b',xs,ys,'gx','LineWidth',2); pause(1); close gcf;
                    end
                    %Truncate boundary points to remove nondrop points (due to needle) and 
                    %convert to meters
                    xData = scale(i)*(x(x8Index:x9Index)); 
                    zData = scale(i)*(y(x8Index:x9Index));
                    %Formulate Initial Guess for fitting parameters
                    Guess = [ST(i) x3*scale(i) -y3*scale(i) Ro Angle];
                    %Neldor-Mead Reloading Loop
                    for j = 1:MaxReload+1
                        tic
                        [P,Error(i)] = fminsearch(@(P) FitLowBondShape(P,xData,-zData,...
                            DeltaRho),Guess);
                        Gamma(i) = P(1); 
                        Xo(i) = P(2);
                        Zo(i) = P(3);
                        R(i) = P(4);
                        Theta(i) = P(5);
                        Difference = max(abs(Guess-P));
                        Guess = P;
                        if flag3
                            %Shift and rotate coordinates
                            xData1 = (xData - Xo(i))*cos(Theta(i)) + ...
                                (-zData - Zo(i))*sin(Theta(i));
                            zData1 = (-zData - Zo(i))*cos(Theta(i)) - ...
                                (xData - Xo(i))*sin(Theta(i));
                            g = 9.8; %Gravity in m/s^2
                            Bond(i) = DeltaRho*g*R(i)^2/Gamma(i); %Bond #
                            system = @(s,y) ShapeODE(s,y,Bond(i));
                            [s,Y] = ode45(system,[0 1.5*pi],[0 0 0]);
                            %Dimensionalize Data for Plotting
                            xTheory = Y(:,1)*R(i);
                            zTheory = Y(:,2)*R(i);
                            %Truncate Data for Plotting and Analysis
                            Cutoff = find(zTheory>=max(zData1));
                            xTheory = xTheory(1:Cutoff);
                            zTheory = zTheory(1:Cutoff);
                            s = s(1:Cutoff);
                            %Compute Volume and Surface Area of Drop
                            Volume(i) = trapz(zTheory,pi*xTheory.^2);
                            Area(i) = trapz(s,2*pi*xTheory);
                            %Generate full drop shape from theoretical data
                            xTheory = [-flipud(xTheory);xTheory];
                            zTheory = [flipud(zTheory);zTheory];
                %             imshow(A); hold on;
                %             xs = [x1;x2;x3;x4;x5;x6;x7;x8;x9]; 
                %             ys = [y1;y2;y3;y4;y5;y6;y7;y8;y9];
                            plot(xData1*1000,zData1*1000,'ko',xTheory*1000,...
                                zTheory*1000,'r','LineWidth',2); 
                            xlabel('X [mm]','FontSize',20);
                            ylabel('Z [mm]','FontSize',20);
                            axis equal; pause(1);
                        end
                        toc
                        if Difference < 1e-5; 
                            disp(['Image ' num2str(i) ' of ' num2str(N) ...
                                ' analyzed after ' num2str(j-1) ' reloads.'])
                %             saveas(gcf,[folder '/' 'Fit ' file num2str(i)],'fig')
                            break; 
                        end
                    end
                end
                end

            end

            Gamma = mean(ST(:))*1000;
            scale = mean(scale(:));
%             disp(['Sigma = ' num2str(mean(ST(:))*1000,3) ' +/- ' ...
%                 num2str(std(ST(:))*1000,1) ' mN/m'])
        end
    
    function ST = SimpleSTCalculation(S,D,DeltaRho)
        % Refer: M Misak, Equations for Detemining l/H Versus S Values for Interfacial Tension Calculations by the Pendant Drop Methods
        if S >= .9
            Hinv = (.30715/S^2.84636) + (-.69116*S^3)-(-1.08315*S^2)+...
                (-.18341*S)-(.20970) ;
        elseif S >= .68
            Hinv = (.31345/S^2.64267) - (.09155*S^2)+(.14701*S)-(.05877); 
        elseif S >= .59
            Hinv = (.31522/S^2.62435) - (.11714*S^2)+(.15756*S)-(.05285);
        elseif S >= .46
            Hinv = (.31968/S^2.59725) - (.46898*S^2)+(.50059*S)-(.13261);
        elseif S >= .401
            Hinv = (.32720/S^2.56651) - (.97553*S^2)+(.84059*S)-(.18069);
        else
            disp('Shape is too spherical');
            %Use formula for S > 0.401 even though it is wrong
            Hinv = (.32720/S^2.56651) - (.97553*S^2)+(.84059*S)-(.18069);
        end
    ST = (DeltaRho*9.8*(D)^2)*(Hinv);
    end