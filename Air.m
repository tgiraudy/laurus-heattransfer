classdef Air < Fluid
   properties
    gbetanu2
   end
    
   methods
      
       function this = Air()
          
           this.Pr = 0.7;
           this.k = 23e-3; % [W/mK]
           this.gbetanu2 = 1.85e8; 
           this.U = 1; % [m/s]
           this.rho = 1.184; % [kg/m3] @25ºC
           this.mu = 18.37e-6; % [Pa s] @25ºC
           this.Cp = 1006; % [J / (Kg K)] 
           
       end
       
   end
    
end