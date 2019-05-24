 function [S,D,DropContour,scale] = findDropContour(im,needleD,NPoints)
        %Inputs:   im - image containing the capillary and the drop
        %          needleDia - the diameter of the drop in *meters*  
        %          Npoints - *optional* the number of points to mark
        %          the sides of the capillary
        %Outputs:  D - The equitorial diameter of the drop
        %          S - Ratio of diameter at distance D from the apex of the
        %          drop to D
        %          Dropcontour - 2D array with drop boundary points
        %          excluding the capillary scaled to physical units
        %          scale - scale factor for Pixel to meters conversion 
        % Author: Vinny (vineethcs.cet@gmail.com) (Based on orginal code by
        % J. Frostad (john.frostad@ubc.ca))
        
        % Check validity of inputs
        if (~isnumeric(needleD))
            error("Enter numeric values (in meters) for the needle diameter");
        end
        
        % Handle optional arguments:
        if nargin<3
            NPoints=25;
        end
        TolDev = 0.1;
        TolPixShift = 7;
        
        % Extract contours
                %Use Otsu's method to chose dividing line between black and white
                Intensity = graythresh(im);  
                BW = imcomplement(im2bw(im,Intensity)); %convert to a binary image
                %Trace drop boundary using bwtraceboundary, the outer boundary of the 
                %drop should be all dark with no glares.
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
                Dpix =(x2-x1);
                
                y4Index = find(y>=y3-Dpix,1); y4 = y(y4Index); x4 = x(y4Index);
                y5Index = find(y>=y3-Dpix,1,'last'); y5 = y(y5Index); x5 = x(y5Index);
                
                x6 = x(1:NPoints); y6 = y(1:NPoints);
                y7Index = find(y(1:end-5)>=y6(1),NPoints,'last'); 
                x7 = x(y7Index);
                y7 = y (y7Index);
                                
                %Calculate offset angle of needle from perpendicular to image edge
                p = polyfit(y6,x6,1); slope = p(1); intercept = p(2);
                Angle = atan(slope);
                %Use pixel width of needle (corrected for offset angle) and known
                %physical width of needle to get conversion factor.
                scale = needleD/(mean(x7-x6)*cos(Angle));                 
                
                %Calculate drop co-ordinates
                 x8Index = find(abs((x-(slope*y+intercept))./(slope*y+intercept))>= TolDev,1);
                 y8Index = find(y==y(x8Index)+TolPixShift,1); % Shift up by 5 pixels (can cast as a Tol)
                 x8 = x(y8Index); y8 = y(y8Index);                   
                 y9Index = find(y>y8,1,'last'); % original code,
                 x9 = x(y9Index); y9 = y(y9Index);
                              
                % Parameters to return
                D = (x2-x1)*scale; %max drop width
                S = mean(x5-x4)/(x2-x1); %drop width at max drop width above apex normalized by D
                DropContour = scale*[x(y8Index:y9Index),y(y8Index:y9Index)];
                     
 end       