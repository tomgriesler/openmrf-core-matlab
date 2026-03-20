function dat = fft1d( dat, dim, zf)


if nargin == 3
    a = ones(1,ndims(dat));
    a(1,dim) = zf;
    dat = matresize(dat,round(size(dat).*a./2)*2);
end

dat = fftshift( dat, dim ); 
dat = fft( dat, [], dim); 
dat = fftshift( dat, dim); 

%fftshift(fft(fftshift(dat,dim),[],dim),dim);

end