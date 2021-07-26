classdef AirProperties
   properties
    Pr
    k
    gbetanu2
   end
    
   methods
      
       function this = AirProperties()
          
           this.Pr = 0.7;
           this.k = 23e-3; % [W/mK]
           this.gbetanu2 = 1.85e8; 
           
       end
       
   end
    
end