function dat = ifft1d( dat, dim, zf )

if nargin == 3
    a = ones(1,ndims(dat));
    a(1,dim) = zf;
    dat = matresize(dat,round(size(dat).*a./2)*2);
end


dat = fftshift(dat,dim);
dat = ifft(dat,[],dim);
dat = fftshift(dat,dim);

end