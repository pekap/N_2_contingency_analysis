classdef Sz
   % size class makes comfortable to calculate row and columns of a matrix
   methods (Static)
      function r = r(matrix) % rows of a matrix
         D=size(matrix);
         r=D(1);
      end
      function c = c(matrix) % columns of a matrix
         D=size(matrix);
         c=D(1);
      end
   end
   
end

