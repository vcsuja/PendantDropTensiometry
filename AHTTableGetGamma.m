function ST = AHTFitGetGamma(S,D,DeltaRho)
        %Params:  D - The equitorial diameter of the drop
        %         S - Ratio of diameter at distance D from the apex of the
        %         drop to D
        %         DeltaRho - Density difference between the suspended fluid
        %         and the ambient
        
        % Function utilizes the tabulated H-S function for the 
        % pendant drops developed by Andreas, Hauser and Tucker (1938) as
        % presented by Bidwell
        % Refer: (i) Bidwell et.al Tables for the determination of the
        % surface tensions of liquid metals by the pendant drop method
        %        (ii) J.M Andreas et.al Boundary tension by Pendant drops
        
        % Author: Vinny (vineethcs.cet@gmail.com)
  
        
        %Handle optional arguments:
        if nargin<4
            g = 9.8;
        end
        
% The S v/s 1/H data          
STable = 0.2:0.001:1.1;        
HinvTable = [19.3532500000000,19.1563300000000,18.9608900000000,18.7669200000000,18.5744300000000,18.3834300000000,18.1939000000000,18.0058400000000,17.8192700000000,17.6341800000000,17.4505600000000,17.2684200000000,17.0877600000000,16.9085800000000,16.7308700000000,16.5546500000000,16.3799000000000,16.2066300000000,16.0348400000000,15.8645200000000,15.6956900000000,15.5283300000000,15.3624500000000,15.1980500000000,15.0351300000000,14.8736900000000,14.7137200000000,14.5552400000000,14.3982300000000,14.2427000000000,14.0886400000000,13.9360700000000,13.7849700000000,13.6353500000000,13.4872100000000,13.3405500000000,13.1953700000000,13.0516700000000,12.9094400000000,12.7686900000000,12.6294200000000,12.4916300000000,12.3553100000000,12.2204800000000,12.0871200000000,11.9552400000000,11.8248400000000,11.6959200000000,11.5684800000000,11.4425100000000,11.3180200000000,11.1950100000000,11.0734800000000,10.9534300000000,10.8348500000000,10.7177600000000,10.6021400000000,10.4880000000000,10.3753400000000,10.2641500000000,10.1544500000000,10.0462200000000,9.93947000000000,9.83420000000000,9.73041000000000,9.62809000000000,9.52726000000000,9.42790000000000,9.33002000000000,9.23362000000000,9.13870000000000,9.04525000000000,8.95328000000000,8.86280000000000,8.77379000000000,8.68625000000000,8.60020000000000,8.51562000000000,8.43253000000000,8.35091000000000,8.27077000000000,8.19211000000000,8.11492000000000,8.03921000000000,7.96499000000000,7.89224000000000,8.05228000000000,7.98182000000000,7.91201000000000,7.84283000000000,7.77429000000000,7.70640000000000,7.63914000000000,7.57252000000000,7.50655000000000,7.44121000000000,7.37651000000000,7.31245000000000,7.24903000000000,7.18625000000000,7.12411000000000,7.06261000000000,7.00175000000000,6.94153000000000,6.88195000000000,6.82301000000000,6.76471000000000,6.70704000000000,6.65002000000000,6.59364000000000,6.53790000000000,6.48279000000000,6.42833000000000,6.37450000000000,6.32132000000000,6.26878000000000,6.21687000000000,6.16560000000000,6.11498000000000,6.06499000000000,6.01565000000000,5.96694000000000,5.91887000000000,5.87144000000000,5.82466000000000,5.77851000000000,5.73300000000000,5.68813000000000,5.64390000000000,5.64907000000000,5.60579000000000,5.56287000000000,5.52033000000000,5.47815000000000,5.43635000000000,5.39491000000000,5.35385000000000,5.31316000000000,5.27284000000000,5.23289000000000,5.19331000000000,5.15410000000000,5.11526000000000,5.07679000000000,5.03869000000000,5.00096000000000,4.96361000000000,4.92662000000000,4.89000000000000,4.85376000000000,4.81788000000000,4.78238000000000,4.74725000000000,4.71248000000000,4.67809000000000,4.64407000000000,4.61042000000000,4.57714000000000,4.54423000000000,4.51169000000000,4.47952000000000,4.44772000000000,4.41629000000000,4.38523000000000,4.37119000000000,4.34076000000000,4.31058000000000,4.28065000000000,4.25096000000000,4.22152000000000,4.19232000000000,4.16337000000000,4.13466000000000,4.10620000000000,4.07798000000000,4.05001000000000,4.02229000000000,3.99481000000000,3.96758000000000,3.94059000000000,3.91385000000000,3.88735000000000,3.86110000000000,3.83509000000000,3.80933000000000,3.78382000000000,3.75855000000000,3.73353000000000,3.70875000000000,3.68422000000000,3.65993000000000,3.63589000000000,3.61210000000000,3.58855000000000,3.56524000000000,3.54941000000000,3.52652000000000,3.50381000000000,3.48127000000000,3.45891000000000,3.43672000000000,3.41471000000000,3.39288000000000,3.37123000000000,3.34975000000000,3.32844000000000,3.30732000000000,3.28636000000000,3.26559000000000,3.24499000000000,3.22457000000000,3.20432000000000,3.18425000000000,3.16436000000000,3.14464000000000,3.12510000000000,3.10573000000000,3.08654000000000,3.06753000000000,3.04870000000000,3.03003000000000,3.01155000000000,2.99691000000000,2.97875000000000,2.96073000000000,2.94284000000000,2.92508000000000,2.90745000000000,2.88996000000000,2.87261000000000,2.85538000000000,2.83829000000000,2.82133000000000,2.80451000000000,2.78782000000000,2.77126000000000,2.75484000000000,2.73855000000000,2.72239000000000,2.70636000000000,2.69047000000000,2.67471000000000,2.65909000000000,2.64360000000000,2.62824000000000,2.61301000000000,2.59792000000000,2.58502000000000,2.57018000000000,2.55545000000000,2.54081000000000,2.52629000000000,2.51187000000000,2.49755000000000,2.48334000000000,2.46923000000000,2.45523000000000,2.44133000000000,2.42754000000000,2.41385000000000,2.40026000000000,2.38678000000000,2.37341000000000,2.36014000000000,2.34697000000000,2.33391000000000,2.32096000000000,2.30811000000000,2.29536000000000,2.28272000000000,2.27018000000000,2.25899000000000,2.24664000000000,2.23437000000000,2.22218000000000,2.21008000000000,2.19806000000000,2.18613000000000,2.17428000000000,2.16252000000000,2.15084000000000,2.13924000000000,2.12773000000000,2.11631000000000,2.10497000000000,2.09371000000000,2.08254000000000,2.07145000000000,2.06045000000000,2.04953000000000,2.03870000000000,2.02795000000000,2.01809000000000,2.00750000000000,1.99698000000000,1.98653000000000,1.97614000000000,1.96583000000000,1.95559000000000,1.94542000000000,1.93532000000000,1.92529000000000,1.91532000000000,1.90543000000000,1.89561000000000,1.88586000000000,1.87617000000000,1.86656000000000,1.85702000000000,1.84755000000000,1.83814000000000,1.82881000000000,1.81955000000000,1.81089000000000,1.80175000000000,1.79267000000000,1.78365000000000,1.77468000000000,1.76578000000000,1.75693000000000,1.74815000000000,1.73942000000000,1.73075000000000,1.72214000000000,1.71358000000000,1.70509000000000,1.69666000000000,1.68828000000000,1.67996000000000,1.67170000000000,1.66351000000000,1.65536000000000,1.64766000000000,1.63963000000000,1.63164000000000,1.62371000000000,1.61583000000000,1.60799000000000,1.60021000000000,1.59248000000000,1.58479000000000,1.57716000000000,1.56958000000000,1.56204000000000,1.55456000000000,1.54713000000000,1.53975000000000,1.53241000000000,1.52513000000000,1.51790000000000,1.51099000000000,1.50386000000000,1.49676000000000,1.48971000000000,1.48270000000000,1.47574000000000,1.46882000000000,1.46194000000000,1.45510000000000,1.44831000000000,1.44156000000000,1.43486000000000,1.42820000000000,1.42158000000000,1.41500000000000,1.40847000000000,1.40198000000000,1.39554000000000,1.38934000000000,1.38297000000000,1.37665000000000,1.37036000000000,1.36410000000000,1.35789000000000,1.35171000000000,1.34557000000000,1.33947000000000,1.33341000000000,1.32738000000000,1.32139000000000,1.31544000000000,1.30953000000000,1.30365000000000,1.29781000000000,1.29201000000000,1.28641000000000,1.28067000000000,1.27498000000000,1.26931000000000,1.26368000000000,1.25808000000000,1.25252000000000,1.24699000000000,1.24149000000000,1.23602000000000,1.23059000000000,1.22519000000000,1.21983000000000,1.21449000000000,1.20920000000000,1.20393000000000,1.19882000000000,1.19361000000000,1.18844000000000,1.18329000000000,1.17818000000000,1.17309000000000,1.16803000000000,1.16301000000000,1.15801000000000,1.15304000000000,1.14810000000000,1.14319000000000,1.13831000000000,1.13346000000000,1.12864000000000,1.12385000000000,1.11918000000000,1.11444000000000,1.10973000000000,1.10505000000000,1.10039000000000,1.09576000000000,1.09115000000000,1.08657000000000,1.08202000000000,1.07749000000000,1.07299000000000,1.06851000000000,1.06407000000000,1.05964000000000,1.05525000000000,1.05095000000000,1.04661000000000,1.04228000000000,1.03798000000000,1.03371000000000,1.02946000000000,1.02523000000000,1.02102000000000,1.01684000000000,1.01268000000000,1.00855000000000,1.00444000000000,1.00035000000000,0.996290000000000,0.992250000000000,0.988290000000000,0.984300000000000,0.980320000000000,0.976370000000000,0.972440000000000,0.968530000000000,0.964640000000000,0.960770000000000,0.956920000000000,0.953100000000000,0.949290000000000,0.945510000000000,0.941750000000000,0.938010000000000,0.934340000000000,0.930640000000000,0.926960000000000,0.923300000000000,0.919650000000000,0.916030000000000,0.912430000000000,0.908840000000000,0.905280000000000,0.901730000000000,0.898210000000000,0.894700000000000,0.891210000000000,0.887740000000000,0.884340000000000,0.880910000000000,0.877490000000000,0.874100000000000,0.870720000000000,0.867360000000000,0.864010000000000,0.860690000000000,0.857380000000000,0.854090000000000,0.850810000000000,0.847560000000000,0.844320000000000,0.841100000000000,0.837930000000000,0.834740000000000,0.831570000000000,0.828420000000000,0.825280000000000,0.822160000000000,0.819050000000000,0.815960000000000,0.812880000000000,0.809820000000000,0.806780000000000,0.803750000000000,0.800740000000000,0.797750000000000,0.794800000000000,0.791830000000000,0.788880000000000,0.785950000000000,0.783030000000000,0.780120000000000,0.777230000000000,0.774350000000000,0.771490000000000,0.768640000000000,0.765810000000000,0.762990000000000,0.760190000000000,0.757430000000000,0.754650000000000,0.751890000000000,0.749140000000000,0.746410000000000,0.743690000000000,0.740980000000000,0.738280000000000,0.735600000000000,0.732930000000000,0.730280000000000,0.727640000000000,0.725010000000000,0.722420000000000,0.719820000000000,0.717230000000000,0.714660000000000,0.712090000000000,0.709540000000000,0.707000000000000,0.704470000000000,0.701960000000000,0.699460000000000,0.696970000000000,0.694490000000000,0.692030000000000,0.689590000000000,0.687150000000000,0.684720000000000,0.682310000000000,0.679900000000000,0.677500000000000,0.675120000000000,0.672750000000000,0.670390000000000,0.668040000000000,0.665700000000000,0.663380000000000,0.661080000000000,0.658770000000000,0.656480000000000,0.654200000000000,0.651930000000000,0.649660000000000,0.647410000000000,0.645170000000000,0.642940000000000,0.640730000000000,0.638520000000000,0.636320000000000,0.634130000000000,0.631970000000000,0.629800000000000,0.627650000000000,0.625500000000000,0.623360000000000,0.621240000000000,0.619120000000000,0.617010000000000,0.614910000000000,0.612830000000000,0.610750000000000,0.608680000000000,0.606630000000000,0.604590000000000,0.602550000000000,0.600510000000000,0.598490000000000,0.596480000000000,0.594480000000000,0.592480000000000,0.590500000000000,0.588520000000000,0.586560000000000,0.584600000000000,0.582660000000000,0.580720000000000,0.578790000000000,0.576870000000000,0.574950000000000,0.573050000000000,0.571150000000000,0.569260000000000,0.567380000000000,0.565510000000000,0.563650000000000,0.561790000000000,0.559950000000000,0.558120000000000,0.556290000000000,0.554460000000000,0.552650000000000,0.550840000000000,0.549040000000000,0.547250000000000,0.545470000000000,0.543700000000000,0.541930000000000,0.540170000000000,0.538430000000000,0.536690000000000,0.534950000000000,0.533220000000000,0.531500000000000,0.529790000000000,0.528080000000000,0.526380000000000,0.524690000000000,0.523000000000000,0.521330000000000,0.519660000000000,0.518000000000000,0.516350000000000,0.514700000000000,0.513060000000000,0.511420000000000,0.509790000000000,0.508170000000000,0.506560000000000,0.504950000000000,0.503350000000000,0.501760000000000,0.500170000000000,0.498600000000000,0.497030000000000,0.495460000000000,0.493900000000000,0.492340000000000,0.490800000000000,0.489260000000000,0.487720000000000,0.486190000000000,0.484670000000000,0.483150000000000,0.481650000000000,0.480150000000000,0.478650000000000,0.477160000000000,0.475680000000000,0.474200000000000,0.472720000000000,0.471260000000000,0.469800000000000,0.468340000000000,0.466890000000000,0.465450000000000,0.464020000000000,0.462590000000000,0.461160000000000,0.459740000000000,0.458330000000000,0.456920000000000,0.455510000000000,0.454120000000000,0.452720000000000,0.451340000000000,0.449950000000000,0.448580000000000,0.447220000000000,0.445850000000000,0.444490000000000,0.443140000000000,0.441790000000000,0.440440000000000,0.439100000000000,0.437770000000000,0.436440000000000,0.435120000000000,0.433800000000000,0.432490000000000,0.431180000000000,0.429880000000000,0.428580000000000,0.427290000000000,0.426000000000000,0.424720000000000,0.423440000000000,0.422160000000000,0.420890000000000,0.419630000000000,0.418370000000000,0.417120000000000,0.415870000000000,0.414620000000000,0.413380000000000,0.412140000000000,0.410910000000000,0.409680000000000,0.408460000000000,0.407240000000000,0.406020000000000,0.404810000000000,0.403610000000000,0.402410000000000,0.401210000000000,0.400020000000000,0.398830000000000,0.397640000000000,0.396460000000000,0.395290000000000,0.394110000000000,0.392940000000000,0.391780000000000,0.390620000000000,0.389470000000000,0.388310000000000,0.387160000000000,0.386020000000000,0.384880000000000,0.383740000000000,0.382610000000000,0.381470000000000,0.380350000000000,0.379220000000000,0.378110000000000,0.377000000000000,0.375880000000000,0.374770000000000,0.373670000000000,0.372570000000000,0.371470000000000,0.370370000000000,0.369280000000000,0.368190000000000,0.367110000000000,0.366020000000000,0.364950000000000,0.363880000000000,0.362800000000000,0.361740000000000,0.360670000000000,0.359600000000000,0.358540000000000,0.357490000000000,0.356430000000000,0.355380000000000,0.354330000000000,0.353280000000000,0.352250000000000,0.351210000000000,0.350170000000000,0.349130000000000,0.348100000000000,0.347070000000000,0.346040000000000,0.345010000000000,0.343990000000000,0.342960000000000,0.341950000000000,0.340940000000000,0.339920000000000,0.338910000000000,0.337900000000000,0.336890000000000,0.335880000000000,0.334870000000000,0.333870000000000,0.332860000000000,0.331860000000000,0.330860000000000,0.329880000000000,0.328880000000000,0.327880000000000,0.326890000000000,0.325870000000000,0.324890000000000,0.323910000000000,0.322930000000000,0.321940000000000,0.320950000000000,0.319960000000000,0.318970000000000,0.317970000000000,0.316970000000000,0.315970000000000,0.314970000000000,0.313970000000000,0.312960000000000,0.311950000000000,0.310940000000000,0.309920000000000,0.308910000000000,0.308540000000000,0.307470000000000,0.306360000000000,0.305190000000000,0.303980000000000,0.302730000000000,0.301420000000000,0.300070000000000,0.298670000000000,0.297220000000000,0.295730000000000,0.294190000000000,0.292600000000000,0.290960000000000,0.285810000000000,0.284150000000000,0.282540000000000,0.280970000000000,0.279450000000000,0.277960000000000,0.276520000000000,0.275130000000000,0.273770000000000,0.272460000000000,0.271190000000000,0.269960000000000,0.268780000000000,0.267640000000000,0.266540000000000,0.265480000000000,0.264470000000000,0.263500000000000,0.264030000000000,0.263090000000000,0.262160000000000,0.261230000000000,0.260310000000000,0.259390000000000,0.258480000000000,0.257570000000000,0.256670000000000,0.255770000000000,0.254870000000000,0.253980000000000,0.253090000000000,0.252210000000000,0.251330000000000,0.250460000000000,0.249590000000000,0.248740000000000,0.247880000000000,0.247030000000000,0.246170000000000,0.245330000000000,0.244480000000000,0.243640000000000,0.242810000000000,0.241970000000000,0.241150000000000,0.240320000000000,0.239500000000000,0.238690000000000,0.237880000000000,0.237070000000000,0.236260000000000,0.235460000000000,0.234670000000000,0.233890000000000,0.233110000000000,0.232320000000000,0.231540000000000,0.230760000000000,0.229990000000000,0.229220000000000,0.228460000000000,0.227690000000000,0.226940000000000,0.226180000000000,0.225430000000000,0.224680000000000,0.223940000000000,0.223200000000000,0.222460000000000,0.221730000000000,0.221010000000000,0.220290000000000,0.219570000000000,0.218850000000000,0.218130000000000,0.217420000000000,0.216710000000000,0.216010000000000,0.215300000000000,0.214600000000000,0.213910000000000,0.213220000000000,0.212530000000000,0.211840000000000,0.211160000000000,0.210480000000000,0.209810000000000,0.209130000000000,0.208460000000000,0.207800000000000,0.207140000000000,0.206480000000000,0.205820000000000,0.205170000000000,0.204520000000000,0.203870000000000,0.203230000000000,0.202590000000000];       
        
%Round S to the nearest thousanth decimal place
Hinv = HinvTable(STable == round(S,3));
                
% Calculate tension        
ST = (DeltaRho*g*(D)^2)*(Hinv);
end        
        