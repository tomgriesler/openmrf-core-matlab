function [ im, cmap, wmap] = openadapt( im, doNorm, Rn, wmap, zf, dummy )
%Adaptive recon based on Walsh et al.
%
% [ recon, cmap, wmap] = openadapt( im, doNorm, Rn, wmap, zf )
%
% === This is the combination, if you want to do it manually ====
% recon = squeeze(sum( conj(wmap) .* im, 1) );
% ===============================================================
%
% === Input ===
%
%   im:         multi-channel images to be reconstructed, 2D/3D
%               (#coils x Ny x Nx x Nz)
%
%   doNorm:     TRUE/FALSE -- normalize image intensity?
%   (optional)  default: FALSE
%
%   Rn:         noise covariance matrix (#coils x #coils)
%   (optional)  *OR* pure noise data (#coils x Np   or   Np x #coils)
%
%   wmap:       coil combination weights, if known from anywhere
%   (optional)
%
%   zf:         factor describing how much the image was interpolated with [z]ero[f]illing
%   (optional)  May be set to increase speed.
%
%
%  The optional input parameters may be set as empty arrays to use defaults.
%
%  Example:
%
%    imCombined = openadapt( im, false, [], [], 2);
%
%      * Rn     is ignored (meaning Rn = eye(Ncoils))
%      * wmap   is calculated inside openadapt
%      * zf = 2 is used to increase interpolation step size
%
%
% === Output ===
%
%   recon:      Reconstructed (combined) image
%               (Ny x Nx x Nz)
%
%   cmap:       "Coil maps"
%
%   wmap:       coil combination weights; can be reused later as input
%
%
% === On noise covariance input Rn ===
%
%   Very ugly reconstructions are typical if Rn can hardly be inverted.
%   We try to detect that and help out the reconstruction by regularizing Rn:
%
%       L  = <some small value>;
%       Rn = (Rn  +  L*norm(Rn)*eye(Ncoil)) ./ (1+L);
%
%   The smallest possible L is found for which Rn can be inverted suffiently accurate.
%   If this happens, you will see a helpful message. No message, no cry! :-)
%   -- MiVö
%
%
% Reference:
%
%   Walsh DO, Gmitro AF, Marcellin MW.
%   Adaptive reconstruction of phased array MR imagery.
%   Magn Reson Med. 2000 May;43(5):682-690.
%   http://dx.doi.org/10.1002/(SICI)1522-2594(200005)43:5<682::AID-MRM10>3.0.CO;2-G
%
%    and
%
%   Mark Griswold, David Walsh, Robin Heidemann, Axel Haase, Peter Jakob.
%   The Use of an Adaptive Reconstruction for Array Coil Sensitivity Mapping and Intensity Normalization,
%   Proceedings of the Tenth Scientific Meeting of the International Society for Magnetic Resonance in Medicine pg 2410 (2002)
%
%
%   This function will calculate adaptively estimated coil maps
%   based on the Walsh algorithm for use in either optimal combination of
%   array images, or for parallel imaging applications. The doNorm flag can be
%   used to include the normalization described in the abstract above. This is
%   only optimal for birdcage type arrays, but will work reasonably for many other geometries.
%   This normalization is only applied to the recon'd image, not the coil maps.
%
%
%   The default block size is 4x4. One can also use interpolation to speed
%   up the calculation (see code), but this could cause some phase errors in practice.
%   Just pay attention to what you are doing.
%
%
%   1) This code is strictly for non-commercial applications
%   2) This code is strictly for research purposes, and should not be used in any
%      diagnostic setting.
%

%   10/1/2001  Mark Griswold
%   22/10/2004  MG - Updated help and changed interpolation to be more stable.
%   2010       Michael Völker - improved computation speed
%                             - added zerofill factor "zf"
%                               ---> smoothed images need about the same time as non-zerofilled ones
%                             - if no Rn is specified, skip dispensable calculations
%   2011       Michael Völker - speed increase on multicore due to parfor-loops
%                               (some changes were necessary in order to make it work nicely)
%                               call 'matlabpool open' before running openadapt()
%                             - further speed increase by using mtimesx() if available (http://www.mathworks.com/matlabcentral/fileexchange/25977)
%                             - imresize() accepts "multi-2D" data of size Ny x Nx x Nc,
%                               so we don't need loops there. This is faster *and* more readable.
%
%   2012       Weick / Völker - 2D and 3D possible
%              Michael Völker - updated help text and improved error checking
%                             - added adaptive regularization of Rn so that
%                               inversion of Rn is possible
%                             - Rn input may be the pure noise data, now.
%                               The covariance is then computed here.
%                             - removed license text
%                             - reduced z-blocksize (z-voxel-size is almost *always* way bigger than inplane)
%                             - added Info in help-text how to apply the weights manually
%                             - To match the notation in the original paper, we compute the "right" weights, that
%                               have to be applied as conj(weights) .* im
%
%   2014       Michael Völker - reduced overhead in 3D coilmap intepolation (no meshgrid needed)
%                             - Columning arrays with (:) is beautiful but
%                               inefficient. :-( This actually creates a temp
%                               copy of the entire array => Horrible with large 3D input.
%                               reshape() is fast.
%                             - modified im-variable prior to parfor loop to permit
%                               slicing (less memory overhead)
%                             - Removed mtimesx. Standard mtimes() seems to
%                               be faster than in older releases.
%                             - Removed a few esoteric optims and
%                               introduced new ones.
%                             - a bit of code deduplication with functions
%
%   2015       Michael Völker - Improved the resizing of complex arrays.
%                               Previously, the abs() was resized the "right" way, with bicubic splines,
%                               but the phase was not interpolated at all.
%                                 I found a good way to achieve an accurate bicubic interpolation for
%                               the phase as well now. Details are in complexResize().
%                               This can have a quite dramatic effect on the accuracy of the coil weights
%                               in regions with faster phase variations. The "blockiness" previously
%                               observed in some image regions is now gone. The phase in the combined image
%                               has always been blocky, this is also smooth now.
%                                 It is even possible now to choose ridiculously high step sizes (like 20)
%                               and still receive almost the same reconstruction as with the old value, 2.
%

% ----- basic input handling -------------------------------------------
    if ndims(im) > 4 || ismatrix(im) || ~isnumeric(im)
        error( 'openadapt:MainInput', 'Expected input image must be...\n\t...numeric\n\t...2D or 3D\n\t...matrix size: [coils] x [Ny] x [Nx] x [Nz]' )
    end
    if ~exist( 'doNorm', 'var') || isempty(doNorm)
        doNorm = false;
    end
    if ~exist( 'Rn', 'var')     || isempty(Rn)
        bNoise = false;     % No noise correlation data were passed to the function.
        Rn = [];            % in the parfor loop, Matlab wants an existing Rn, even if it is not used...
    else
        bNoise = true;
    end
    if exist( 'wmap' , 'var' ) && ~isempty(wmap)
        weights_present = true;
    else
        weights_present = false;
    end
    if ~exist( 'zf', 'var')  || isempty(zf)
        zf = 1;
    end
    % In a previous version, we had one more possible input.
    % Let us not break the interface and just ignore whatever we get here.
    if exist( 'dummy', 'var')
        dummy = dummy;  %#ok
        clear dummy
    end
% ----------------------------------------------------------------------

[nc, ny, nx, nz] = size(im);

% Handle nothing-to-do-case fast:
if nc == 1
    im = reshape( im, ny, nx, nz );
    return
end

st  = round( 3 * zf );     % interpolation step size
bsX = 2 * floor(11/2 * zf ) + 1;     % x-block size (odd)
bsY = 2 * floor(11/2 * zf ) + 1;     % y-block size (odd)
bsZ = 2 * floor( 7/2 * zf ) + 1;     % z-block size (odd)

% With actual 3D data we get good statistics for regional ROIs really
% easy. So we can reduce 3D block size a bit and still get lots of samples.
% => reduce x-blocksize, because it's the most expensive in the current
%    implementation: We need bsX copies of the low-res image.
if nz > 10
    bsX = 5;
end

imInput = im;       % keep a reference to orig input
                    % This cheap, no overhead, because that reference is
                    % kept in the calling workspace, anyway.
is3D = (nz > 1);
NumWorkers = inf;

% ----- more error checking --------------------------------------------
    if      numel( doNorm ) ~= 1                        ...     % no scalar
      || ~( islogical(doNorm) || isnumeric(doNorm) )    ...     % no usable data type
      ||    isempty(intersect( doNorm, [0 1] ))                 % something different than true/false or 1/0
        error('openadapt:doNorm',   '2nd argument ''doNorm'' must be a simple Boolean (true/false).' )
    end
    if weights_present && ~isequal( size(wmap), size(im) )
        error('openadapt:wmap',    '4th argument ''wmap'' must have the same size as the input data.')
    end
    if      numel( zf ) ~= 1        ...     % no scalar
      ||   ~isnumeric( zf )         ...     % no usable data type
      ||   ~isreal( zf )            ...     %
      ||    zf < 1
        error('openadapt:zf', '5th argument ''zf'' must be a scalar number >= 1' )
    end
% ----------------------------------------------------------------------

if bNoise
    Rn = noiseRegul( Rn, nc );
end

precis = class(im);                     % single or double precision?

% are we interested in coil maps?
wantCmaps = nargout > 1 || (doNorm && (~weights_present || bNoise));

NstepY = min( ceil( ny/st ), ny );
NstepX = min( ceil( nx/st ), nx );
NstepZ = min( ceil( nz/st ), nz );

% calculation of coil weights
if ~weights_present || wantCmaps

    NumWorkers = decideParPoolStuff( size(im), zf );

    % find coil with maximum intensity for correcting the phase of all
    % of the other coils.
    [~, maxcoil] = max(sum(abs(matrix(im)),2));

    % If necessary, enlarge im so we can index it without errors. Edges will be
    % wrong anyway.
    pad = [NstepY  NstepX  NstepZ];
    pad = mod( pad - rem([ny nx nz], pad), pad );
    if nnz(pad) > 0                                 % [n]umber of [n]on-[z]ero elements
        im = padarray( im, [0 pad], 0, 'post');
        [nc, ny, nx, nz] = size(im);        % update size
    end

    iy = myuint( ny );  % returns 'uint8', 'uint16', ...
    iz = myuint( nz );  %

    % Edges are cropped so the results near the edges of the image could
    % be in error. Not normally a problem. But watch out for aliased regions.
    ymin = cast( max( (st:st:ny) - floor(bsY/2),  1), iy );
    ymax = cast( min( (st:st:ny) + floor(bsY/2), ny), iy );

    st = nz / NstepZ;   % "fix" problems with small nz...
    zmin = cast( max( (st:st:nz) - floor(bsZ/2),  1), iz );
    zmax = cast( min( (st:st:nz) + floor(bsZ/2), nz), iz );

    if isempty(zmin) || isempty(zmax)
        zmin = 1;
        zmax = 1;
    end

    % Prepare im, such that matlab can slice it. Needs more memory.
    % edges in x-direction will be wrong, but *pfffft*...
    idx = bsxfun( @plus, (1:NstepX).' * nx/NstepX, -floor(bsX/2):floor(bsX/2) );
    idx = min(max(idx,1), nx);      % [ 1 1 1 2 3 4 5 ... (nx-2) (nx-1) nx nx nx nx ]
    im = reshape( im(:,:,idx,:), nc, ny, NstepX, [], nz );
    clear  pad  idx

    % different x-coordinates are treated simultaneously on several workers
%    for x = 1:NstepX      % for debugging/profiling
    parfor ( x = 1:NstepX, NumWorkers )

        % helper variables to make parfor-loop faster
        wTmp = complex(zeros( nc, NstepZ, NstepY, precis ));
        ClmTmp = [];        % work around a silly matlab warning
        if bNoise && wantCmaps
            ClmTmp = complex(zeros( nc, NstepZ, NstepY, precis ));
        end
        imBlockX = reshape( im(:,:,x,:,:), nc,ny,[],nz );   % im is nicely sliceable here

        for y = 1:NstepY
            yidx = ymin(y):ymax(y);     %#ok <-- suppress a warning in matlab editor

            for z = 1:NstepZ

                Rs = reshape( imBlockX(:,yidx,:,zmin(z):zmax(z)), nc, [] );     %#ok <-- suppress a warning in matlab editor
                Rs = Rs * Rs';

                % Eigenvector with max eigenvalue gives the correct combination coeffs.
                if bNoise
                    [eivec, eivals] = eig( Rn \ Rs);
                else
                    [eivec, eivals] = eig( Rs );
                end
                [~,ind] = max( diag(eivals) );

                normmf = eivec(:,ind);

                if bNoise
                    if ~weights_present
                        wTmp(:,z,y) = normmf / (normmf' * (Rn \ normmf));   % MiVö: Should the denomninator not be sqrt(normmf' * (Rn \ normmf)) ? --> Eq. (20) in the Walsh-Paper
                    end
                    if wantCmaps
                        ClmTmp(:,z,y) = normmf;
                    end
                else
                    wTmp(:,z,y) = normmf;
                end

            end % for z = 1:NstepZ
        end % for y = 1:NstepY

        if ~weights_present
            wmap(:,x) = reshape( wTmp, [], 1 );
        end

        if bNoise && wantCmaps
          cmap(:,x) = reshape( ClmTmp, [], 1 );
        end
    end % parfor x = 1:NstepX

    im = imInput;
    [nc, ny, nx, nz] = size(im);    % update size

    clear  imInput  Rn  wTmp  ClmTmp  imBlockX  yidx  normmf  eivec  eivals  Rs

    if ~weights_present
        wmap = permute( reshape(wmap,nc,NstepZ,NstepY,NstepX), [3 4 2 1] );  % NY x NX x NZ x NC
        wmap = bsxfun( @times, wmap, exp( -1i .* angle(wmap(:,:,:,maxcoil)) ) );
    end
    if bNoise && wantCmaps
        cmap = permute( reshape(cmap,nc,NstepZ,NstepY,NstepX), [3 4 2 1] );  % NY x NX x NZ x NC
        cmap = bsxfun( @times, cmap, exp( -1i .* angle(cmap(:,:,:,maxcoil)) ) );
    end

end % if ~weights_present || wantCmaps


if is3D && ( ~weights_present  ||  wantCmaps   )
    % When using matlab's "interp3", we need to manually specify *where* matlab
    % should calculate the interpolation.
    X = 1 + (NstepX-1)*(0:nx-1)./(nx-1);
    Y = 1 + (NstepY-1)*(0:ny-1)./(ny-1);
    Z = 1 + (NstepZ-1)*(0:nz-1)./(nz-1);
    Y = reshape( Y, [], 1);
else
    X = []; Y = []; Z = [];
end

% coil-by-coil interpolation of weights so the arrays have the size of the input images.
if ~weights_present
    wmap = complexResize( wmap, X, Y, Z, nc, nx, ny, nz, NumWorkers );
    wmap = permute( wmap, [4 1 2 3] );    % Coils back to first dimension
end

% coil-by-coil interpolation of coilmaps so the arrays have the size of the input images.
if bNoise && wantCmaps
    cmap = complexResize( cmap, X, Y, Z, nc, nx, ny, nz, NumWorkers );
    cmap = permute( cmap, [4 1 2 3] );  % Coils back to first dimension
end

clear X Y Z

% Combine coil signals:
im = sum(conj(wmap).*im, 1);
im = reshape( im, ny, nx, nz );

if ~bNoise          % we don't know Rn
    cmap = wmap;    % ---> then there is no difference between cmap and wmap
end
if nargout < 3
    clear wmap
end

if doNorm
    % This is the normalization proposed in the abstract
    % referenced in the header.
    im = im .* reshape( sum(abs(cmap),1).^2, ny, nx, nz );
end

end % of openadapt()



% =========== HELPERS =====================================================

function Rn = noiseRegul( Rn, nc )
    % Check if noise input is a covariance matrix
    % or the pure noise data.
    if isequal( size(Rn), [nc nc] )             % square matrix ncoil x ncoil
        % everything is fine
    elseif ismatrix(Rn) && size(Rn,1) == nc

        Rn = Rn * Rn' ./ size(Rn,2);            % correlation matrix (zero-mean noise)

    elseif ismatrix(Rn) && size(Rn,2) == nc

        Rn = conj(Rn' * Rn)  ./ size(Rn,1);     % correlation matrix (zero-mean noise)
    else
        error('openadapt:Rn', '3rd argument ''Rn'' must be a noise correlation matrix of size Ncoils x Ncoils or the direct noise data.' )
    end

    % Test the condition number and find an amount of regularization
    % so that the condition number does not exceed a threshold.
    cMax  = 5;                      % maximum tolerated condition number
    prec  = class(Rn);
    normR = norm(Rn);
    eiei  = eye( nc, prec );
    L     = 0;
    n     = log2(eps(normR));       % eps('double') = 2^-52
    maxN  = log2(realmax(prec));    % realmax('double') = 2^1024

    while  cond(Rn + L*normR*eiei) > cMax  &&  n < maxN
        L = 2^n;
        n = n + 1;
    end

    if L > 0
        L = max( L, max(eps(col(Rn))) );
        fprintf('Mixing %.3g%% * norm(Rn) * eye(N) to noise correlation matrix for regularization.\n', 100*L )
        Rn = ( Rn + L*normR*eiei );
        Rn = Rn .* (normR ./ norm(Rn) );    % don't change Rn's norm
    end
end % of noiseRegul()

function uIntType = myuint( n )
    mx = max(col(n));

    if mx <= intmax('uint8')
        uIntType = 'uint8';
    elseif mx <= intmax('uint16')
        uIntType = 'uint16';
    elseif mx <= intmax('uint32')
        uIntType = 'uint32';
    else
        uIntType = 'uint64';
    end
end % of myuint()

function x = col(x)
    x = reshape(x,[],1);   % *not* x(:)!
end

function x = matrix( x )
    x = reshape( x, size(x,1), [] );
end

function A = complexResize( A, X, Y, Z, nc, nx, ny, nz, NumWorkers )

    evalin( 'caller', ['clear ' inputname(1) ] )

    % How to treat complex data?
    %   a) resize real() and imag()  seperately, combine to full image later
    %   b) resize abs()  and angle() seperately, combine to full image later
    %   c) resize abs(). To treat the phase, do a) and only keep the resulting
    %                    phase.
    %
    % Use c).
    %
    % Why not a)?
    %   It's not energy-conserving. sqrt(R^2 + I^2) != R + I
    %   Areas with real(A) > 0 and those with real(A) < 0 are interpolated in a way
    %   that leads to signal voids. The same for imag(A).
    %
    % Why not b)? (the method used until 2015)
    %   angle(A) is a nonlinear transformation applied to A. The resizing algorithm
    %   respects energy conservation in magnitude data, i.e. each numerical value is interpreted
    %   as an intensity. This leads to silly interpolations of the phase, e.g.:
    %       (angle(1) + angle(-1))/2 = pi/2 -> 1i
    %   The workaround until 2015 was to use the dumb nearest-'interpolation',
    %   that is, each pixel is just resized.
    %
    absA = abs(A);
    R = real( A );
    I = imag( A );
    clear A

    if nz > 1
        parfor ( C=1:nc, NumWorkers )
            A(:,:,:,C) = interp3( absA(:,:,:,C),X,Y,Z,'cubic',0);

            phase(:,:,:,C) = atan2(                                         ...
                                    interp3( I(:,:,:,C),X,Y,Z,'cubic',0) , ...
                                    interp3( R(:,:,:,C),X,Y,Z,'cubic',0)   ...
                                  );
        end
    else
        % 2D (imresize does all channels at once)
        smY = size( absA, 1);   % "small y"
        smX = size( absA, 2);   % "small x"

        absA = reshape( absA, smY, smX, nc );   % \
        R    = reshape(    R, smY, smX, nc );   %    workaround for old matlab
        I    = reshape(    I, smY, smX, nc );   % /

        A = imresize( absA, [ny nx],'bicubic');

        phase = atan2(                                    ...
                        imresize( I, [ny nx],'bicubic'), ...
                        imresize( R, [ny nx],'bicubic')  ...
                     );

        A     = reshape(     A, ny, nx, 1, nc );
        phase = reshape( phase, ny, nx, 1, nc );
    end

    A = A .* exp( 1i .* phase );    % voilà :-)

end % of complexResize

% This is how it was done before 2015-02-13 for reference and easy comparison:
function A = complexResizeOld( A, X, Y, Z, nc, nx, ny, nz, NumWorkers )

    evalin( 'caller', ['clear ' inputname(1) ] )

    % call abs() and angle() once before the 3D-Loop for speed
    absA = abs(A);
    phsA = angle(A);
    clear A

    if nz > 1
        parfor ( C=1:nc, NumWorkers )
            A(:,:,:,C) = interp3( absA(:,:,:,C),X,Y,Z,'cubic',0) .* exp( 1i .* interp3( phsA(:,:,:,C), X,Y,Z, 'nearest',0));
        end
    else
        % 2D (imresize does all channels at once)
        smY = size( absA, 1);   % "small y"
        smX = size( absA, 2);   % "small x"
        absA = reshape( absA, smY, smX, nc );   % workaround for old matlab
        phsA = reshape( phsA, smY, smX, nc );   %
        A = imresize( absA, [ny nx],'bicubic') .* exp( 1i .* imresize( phsA, [ny nx], 'nearest'));
        A = reshape( A, ny, nx, 1, nc );
    end

end % of complexResizeOld

function NumWorkers = decideParPoolStuff( szIm, zf )
% returns 0 if no parfor is desired/possible

    NumWorkers = 0;

    % Define minimum number of pixels for which auto-opening a parpool
    % makes sense to counterbalance delay vs. speed.
    % (merely a wild guess, of course)
    NumPx4ParPool = 30 * 256 * 256;

    if exist('matlabpool','file') < 2  % we have no toolbox
        return
    end

    try
        foo = gcp('nocreate');
        if ~isempty(foo)    % newer matlab
            NumWorkers = foo.NumWorkers;
        end
    catch
        % old matlab
        NumWorkers = matlabpool('size');    %#ok <-- suppress a warning in matlab editor
    end

    if NumWorkers == 0  &&  prod(szIm(2:end))/zf >= NumPx4ParPool
        warning( 'off', 'parallel:cluster:DepfunError' )
        try
            foo = parpool('local');
            NumWorkers = foo.NumWorkers;
        catch
            % old matlab
            matlabpool open                     %#ok <-- suppress a warning in matlab editor
            NumWorkers = matlabpool('size');    %#ok <-- suppress a warning in matlab editor
        end
    end
end % of decideParPoolStuff()
