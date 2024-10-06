function [x,w] = qnwlege1(n,a,b)

maxit = 100;
m = fix((n+1)/2);
xm = 0.5*(b+a);
xl = 0.5*(b-a);
x = zeros(n,1);
w = x;
i = (1:m)';
z = cos(pi*(i-0.25)./(n+0.5));
for its=1:maxit
   p1 = 1;
   p2 = 0;
   for j=1:n
      p3 = p2;
      p2 = p1;
      p1 = ((2*j-1)*z.*p2-(j-1)*p3)./j;
   end
   pp = n*(z.*p1-p2)./(z.*z-1);
   z1 = z;
   z = z1-p1./pp;
   if abs(z-z1)<1e-14
      break;
   end
end
if its==maxit
   error('Maximum iterations in qnwlege1')
end
x(i) = xm-xl*z;
x(n+1-i) = xm+xl*z;
w(i) = 2*xl./((1-z.*z).*pp.*pp);
w(n+1-i) = w(i);
end
