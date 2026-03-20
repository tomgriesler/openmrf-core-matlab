function dat = fft3d( dat )
%3D Fast Fourier Transform with MRI convention (image centre centred)
%
% image = fft3d( data );
%
% input-size:      [X1  X2  X3]
%                    or
%               [N  X1  X2  X3]
%
%           ==> spatial dimensions at the end

% reworked and accelerated for low-RAM-case
% Michael.Voelker@mr-bavaria.de, 2012

    sz      = size(dat);
    Ndim    = ndims( dat );

    shift(1, Ndim - (0:2)) = sz(Ndim - (0:2)) / 2;

    freeRAM  = getFreeMemory();
    byteFac  = 4 * (2 - isa(dat,'single')) * (2 - isreal(dat));     % one data point in dat needs byteFac Bytes
    datBytes = byteFac  * numel(dat);


    havePlentyRAM = ( freeRAM  >  2 * datBytes );


    if havePlentyRAM

        dat = cmshiftnd(  fft(fft(fft(  cmshiftnd( dat, shift ),  [],Ndim-2),[],Ndim-1),[],Ndim),  shift );

        % We're done. Go back home.
        return
    end


    % =====================================================================
    %
    % This Code is for the case that we have to be careful with RAM usage
    %
    % =====================================================================


    switch Ndim
        case 3
            % adjust size
            sz    = [ 1  sz   ];
            shift = [ 0 shift ];
            dat = reshape( dat, sz );
        case 4
            % fine
        otherwise
            error( 'fft3d:WrongNumDim', 'Input with 3 or 4 dimensions, please, with fft-dimensions at the end.' )
    end


    % Use a loop to save memory
    for step = 1:sz(1)          %#ok  <-- do not tell us to use parfor, here...

        tmp             = dat(step,:,:,:);          % no optimal memory access, but well...
        tmp             = cmshiftnd( tmp, shift );
        tmp             = fftn( tmp );
        tmp             = cmshiftnd( tmp, shift );
        dat(step,:,:,:) = tmp;

    end

    if Ndim == 3
       dat = reshape( dat, sz(2:end) );
    end

end % of fft3d()