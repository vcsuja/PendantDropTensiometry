   function ST = AHTFitGetGamma(S,D,DeltaRho)
        %Params:  D - The equitorial diameter of the drop
        %         S - Ratio of diameter at distance D from the apex of the
        %         drop to D
        %         DeltaRho - Density difference between the suspended fluid
        %         and the ambient
        
        % Function utilizes the fits to the tabulated H-S function for the 
        % pendant drops by Andreas, Hauser and Tucker (1938)
        % Refer: (i) M Misak, Equations for Detemining l/H Versus S Values for Interfacial Tension Calculations by the Pendant Drop Methods
        %        (ii) J.M Andreas et.al Boundary tension by Pendant drops
        
        % Author: Vinny (vineethcs.cet@gmail.com)
        
        
        
        if S > .9
            Hinv = (.30715/S^2.84636) + (-.69116*S^3)-(-1.08315*S^2)+...
                (-.18341*S)-(.20970) ;
        elseif S > .68
            Hinv = (.31345/S^2.64267) - (.09155*S^2)+(.14701*S)-(.05877); 
        elseif S > .59
            Hinv = (.31522/S^2.62435) - (.11714*S^2)+(.15756*S)-(.05285);
        elseif S > .46
            Hinv = (.31968/S^2.59725) - (.46898*S^2)+(.50059*S)-(.13261);
        elseif S > .4
            Hinv = (.32720/S^2.56651) - (.97553*S^2)+(.84059*S)-(.18069);
        elseif S>=  .3 
            Hinv = (.34074/S^2.52303) + (123.9495*S^5)  - (72.82991*S^4) + (0.01320*S^3) - (3.38210*S^2) + (5.52969*S) - 1.07260;
        else
            disp('Shape is too spherical');
            %Use formula for S = 0.3 to S = 0.4
            Hinv = (.34074/S^2.52303) + (123.9495*S^5)  - (72.82991*S^4) + (0.01320*S^3) - (3.38210*S^2) + (5.52969*S) - 1.07260;
        end
    ST = (DeltaRho*9.8*(D)^2)*(Hinv);
    end