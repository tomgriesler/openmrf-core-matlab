function res = times(a,b)

if isa(a,'Wavelet') == 0
    error('In  A.*B only A can be Wavelet operator');
end

imSize_dyd = [max(2.^ceil(log2(a.imageSize(1:2)))), max(2.^ceil(log2(a.imageSize(1:2))))];
frames = size(b,3);

if a.adjoint % Wavelet domain to image
    res = zeros([a.imageSize(1:2) frames]);
    for frameidx = 1:frames
        tmp = b(:,:,frameidx);
        I = IWT2_PO(real(tmp),a.wavScale,a.qmf) + 1i*IWT2_PO(imag(tmp),a.wavScale,a.qmf);
        res(:,:,frameidx) = crop(I,a.imageSize(1:2));
    end
else % image to Wavelet domain
    res = zeros([imSize_dyd frames]);
    for frameidx = 1:frames
        tmp = b(:,:,frameidx);
        tmp = zpad(tmp,imSize_dyd);
        res(:,:,frameidx) = FWT2_PO(real(tmp),a.wavScale,a.qmf) + 1i* FWT2_PO(imag(tmp),a.wavScale,a.qmf);         
    end
end

