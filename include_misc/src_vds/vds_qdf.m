function [roots] = vds_qdf(a,b,c)

d = b^2 - 4*a*c;

roots(1) = (-b + sqrt(d))/(2*a);
roots(2) = (-b - sqrt(d))/(2*a);

end

