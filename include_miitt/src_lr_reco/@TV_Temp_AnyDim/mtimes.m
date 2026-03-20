function res = mtimes(a,b)


if a.adjoint
    res = adjDz(b,a.dim1);
else
    if a.dim1==6
        res = b(:,:,:,:,:,[2:end,end]) - b;
    elseif a.dim1==5
        res = b(:,:,:,:,[2:end,end],:) - b;
    elseif a.dim1==4
        res = b(:,:,:,[2:end,end],:,:) - b;
    elseif a.dim1==3
        res = b(:,:,[2:end,end],:,:,:) - b;
	elseif a.dim1==2
        res = b(:,[2:end,end],:,:,:,:) - b;
	elseif a.dim1==1
        res = b([2:end,end],:,:,:,:,:) - b;
	
    end
end

function y = adjDz(x,dim1)
if dim1==6
    y= x(:,:,:,:,:,[1,1:end-1]) - x;
    y(:,:,:,:,:,1) = -x(:,:,:,:,:,1);
    y(:,:,:,:,:,end) = x(:,:,:,:,:,end-1);
elseif dim1==5
    y= x(:,:,:,:,[1,1:end-1],:) - x;
    y(:,:,:,:,1,:) = -x(:,:,:,:,1,:);
    y(:,:,:,:,end,:) = x(:,:,:,:,end-1,:);
elseif dim1==4
    y= x(:,:,:,[1,1:end-1],:,:) - x;
    y(:,:,:,1,:,:) = -x(:,:,:,1,:,:);
    y(:,:,:,end,:,:) = x(:,:,:,end-1,:,:);
elseif dim1==3
    y= x(:,:,[1,1:end-1],:,:,:) - x;
    y(:,:,1,:,:,:) = -x(:,:,1,:,:,:);
    y(:,:,end,:,:,:) = x(:,:,end-1,:,:,:);
elseif dim1==2
    y= x(:,[1,1:end-1],:,:,:,:) - x;
    y(:,1,:,:,:,:) = -x(:,1,:,:,:,:);
    y(:,end,:,:,:,:) = x(:,end-1,:,:,:,:);
else % dim1==1
    y= x([1,1:end-1],:,:,:,:,:) - x;
    y(1,:,:,:,:,:) = -x(1,:,:,:,:,:);
    y(end,:,:,:,:,:) = x(end-1,:,:,:,:,:);
end

