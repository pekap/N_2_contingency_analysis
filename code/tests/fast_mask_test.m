% short script 

for i=1:Sz.r(Mask1)
   for j=i:Sz.r(Mask1) 
         if (Mask1(i,j) == 1)&&(Mask2(i,j) ~= 1)   
            fprintf('\n ooops %d %d',i,j); 
         end
   end
end


q=0;