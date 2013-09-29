classdef A_class
   properties 
      N
      M
      av
   end
   methods
      function obj = A_class(M,N)
         obj.M=M;
         obj.N=N;
         obj.av = obj.average(10);
      end
      function average=average(obj,C)
         average = (obj.N+obj.M+C)/2;
         fprintf('%i\n',Av_class.sum(obj));
         fprintf('%i',Av_class.sub(obj));
      end
   end
end
