a = cont_brute_force_algorithm;
b = cont_fast_algorithm;

q=[]; kq=0; % a-b
w=[]; kw=0; % b-a
e=[]; ke=0; % a&b

for i=1:Sz.r(a)
   found = (b(:,1)==a(i,1))&(b(:,2)==a(i,2));
   if sum(found)==0
      kq=kq+1;
      q(kq,1:3)=a(i,1:3);
   end
end

for i=1:Sz.r(b)
   found = (a(:,1)==b(i,1))&(a(:,2)==b(i,2));
   if sum(found)==0
      kw=kw+1;
      w(kw,1:3)=b(i,1:3);
   end
end

% q - array of excessive contingencies