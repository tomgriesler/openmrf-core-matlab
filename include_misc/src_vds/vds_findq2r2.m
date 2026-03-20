% =============== CVS Log Messages ==========================
%	This file is maintained in CVS version control.
%
%	$Log: vds.m,v $
%	Revision 1.5  2004/04/27 18:08:44  brian
%	Changed FOV to a polynomial of unlimited length,
%	and hopefully changed all comments accordingly.
%	Also moved sub-functions into vds.m so that
%	no other .m files are needed.
%	
%	Revision 1.2  2003/05/29 23:02:21  brian
%	minor edits
%	
%	Revision 1.1  2002/03/28 01:03:20  bah
%	Added to CVS
%	
%
% ===========================================================

function [q2,r2] = vds_findq2r2(smax,gmax,r,r1,T,Ts,N,Fcoeff,rmax)

gamma = 4258;			% Hz/G

F = 0;		% FOV function value for this r.
dFdr = 0;		% dFOV/dr for this value of r.
for rind = 1:length(Fcoeff)
	F = F+Fcoeff(rind)*(r/rmax)^(rind-1);
	if (rind>1)
		dFdr = dFdr + (rind-1)*Fcoeff(rind)*(r/rmax)^(rind-2)/rmax;
	end;
end;

GmaxFOV = 1/gamma /F/Ts;		% FOV limit on G
Gmax = min(GmaxFOV,gmax);	%

maxr1 = sqrt((gamma*Gmax)^2 / (1+(2*pi*F*r/N)^2));  


if (r1 > maxr1)			
			% Grad amplitude limited.  Here we
			% just run r upward as much as we can without
			% going over the max gradient.
	r2 = (maxr1-r1)/T; 
	%tt = sprintf('Grad-limited r=%5.2f, r1=%f',r,r1);
	%disp(tt);

else

	twopiFoN = 2*pi*F/N;
	twopiFoN2 = twopiFoN^2;

	%	A,B,C are coefficents of the equation which equates
	% 	the slew rate calculated from r,r1,r2 with the
	%	maximum gradient slew rate.
	%
	%	A*r2*r2 + B*r2 + C  =  0	
	%
	%	A,B,C are in terms of F,dF/dr,r,r1, N and smax.
	%


	A = 1+twopiFoN2*r*r;
	B = 2*twopiFoN2*r*r1*r1 + 2*twopiFoN2/F*dFdr*r*r*r1*r1;
	C = twopiFoN2^2*r*r*r1^4 + 4*twopiFoN2*r1^4 + (2*pi/N*dFdr)^2*r*r*r1^4 + 4*twopiFoN2/F*dFdr*r*r1^4 - (gamma)^2*smax^2;


	[rts] = vds_qdf(A,B,C);	% qdf = Quadratic Formula Solution.
	r2 = real(rts(1));	% Use bigger root.  The justification
				% for this is not entirely clear, but
				% in practice it seems to work, and 
				% does NOT work with the other root.




	% Calculate resulting slew rate and print an error 
	% message if it is too large.
	
	slew = 1/gamma*(r2-twopiFoN2*r*r1^2 + i*twopiFoN*(2*r1^2 + r*r2 + dFdr/F*r*r1^2));
	%tt = sprintf('Slew-limited r=%5.2d  SR=%f G/cm/s',r,abs(slew));
	%disp(tt);
	sr = abs(slew)/smax;

	if (abs(slew)/smax > 1.01)
		tt = sprintf('Slew violation, slew = %d, smax = %d, sr=%f, r=%f, r1=%f',round(abs(slew)),round(smax),sr,r,r1);
		disp(tt);
	end;

end;


%	Calculate q2 from other pararmeters.

q2 = 2*pi/N*dFdr*r1^2 + 2*pi*F/N*r2;

end

