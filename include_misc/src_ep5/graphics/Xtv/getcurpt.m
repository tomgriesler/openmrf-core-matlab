function [x,y]=getcurpt(arg1)

point=get(arg1,'currentpoint');

x=point(1,1);
y=point(1,2);
