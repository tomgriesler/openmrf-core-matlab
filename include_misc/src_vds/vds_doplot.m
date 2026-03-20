function vds_doplot(x,y,labs,col,acq)

axl={'x','y','z'};
ss = size(y);
if (min(size(acq))>0)
	mina = min(find(acq > 0));
	maxa = max(find(acq > 0));
	xsq = [x(mina) x(maxa) x(maxa) x(mina)];
	ysq = 1.1*sqrt(max(abs(sum(y'.*y'))))*[-1 -1 1 1];
	patchc = .8*[1 1 1];
end;
hold off;
if (min(size(acq))>0)
	h=patch(xsq,ysq,patchc);
	set(h,'EdgeColor',patchc);
	hold on;
end;
for q=1:ss(2)
	plot(x,y(:,q),col{q});
	hold on;
	leg{q} = sprintf('%s_%s',labs,axl{q});
end;
if (ss(2)>1)
	ab = sqrt(sum(y.'.*y.')).';
	plot(x,ab,col{ss(2)+1});
	leg{ss(2)+1}=sprintf('|%s|',labs);
end;
hold off;
legend(leg);
grid on;

end
