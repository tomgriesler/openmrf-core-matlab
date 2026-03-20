function dist=calcdist(arg1)

if nargin==0, arg1=gca;end
points=drawline(arg1);

dist=sqrt((points(2,1)-points(1,1)).^2+(points(2,2)-points(1,2)).^2);
a=line([points(1,1) points(2,1)],[points(1,2) points(2,2)])