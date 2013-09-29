classdef Av_class
   %AV_CLASS Summary of this class goes here
   %   Detailed explanation goes here
   
   properties
   end
   
   methods (Static)
      function sub = sub(obj)
         sub = obj.N - obj.M;
      end
      function sum = sum(obj)
         sum = obj.N + obj.M;
      end
   end
end

