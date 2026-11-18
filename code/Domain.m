classdef Domain<handle
   properties
        uid=[];%unique id for isSame()
      region;
      phase;
      relVol;
      A_primM=[];
      B_primM=[];
   end
   methods
        function domain=Domain(phaseP)
            if nargin>0
                 domain.phase=phaseP;
                domain.A_primM=phaseP.A_primM;
              end
        end
   end
end
