function rf = SIGPY_slr_dzrf(n, tb, ptype, ftype, d1, d2, cancel_alpha_phs)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%SIGPY_SLR_DZRF  MATLAB port of sigpy.mri.rf.slr.dzrf LIMITED TO ftype='ls'.
%
% This function behaves like sigpy's dzrf for the LS (firls) case and all
% pulse types (ptype): 'st','ex','se','inv','sat'.
%
% For other ftypes ('ms','pm','min','max'), this function errors out with a
% message explaining they were not robust after Python->MATLAB porting
% (due to solver/numerical differences, especially remez/minphase factorization).
%
% Inputs (same semantics as sigpy):
%   n (int)                  : number of time points (even recommended)
%   tb (float)               : time-bandwidth product
%   ptype (char/string)      : 'st','ex','se','inv','sat'
%   ftype (char/string)      : MUST be 'ls'
%   d1 (float)               : passband ripple
%   d2 (float)               : stopband ripple
%   cancel_alpha_phs (bool)  : only used for ptype='ex'
%
% Output:
%   rf (complex column vector) : RF waveform
%
% Requirements:
%   - Signal Processing Toolbox (firls)
%
% Notes:
%   - Keeps structure and naming close to sigpy source.
%   - Only removes functions that are exclusively needed for non-LS ftypes.
%   - Keeps b2rf/b2a/mag2mp/ab2rf unchanged (core SLR transform).

    % ---------------- defaults (match sigpy signature) ----------------
    if nargin < 1 || isempty(n), n = 64; end
    if nargin < 2 || isempty(tb), tb = 4; end
    if nargin < 3 || isempty(ptype), ptype = "st"; end
    if nargin < 4 || isempty(ftype), ftype = "ls"; end
    if nargin < 5 || isempty(d1), d1 = 0.01; end
    if nargin < 6 || isempty(d2), d2 = 0.01; end
    if nargin < 7 || isempty(cancel_alpha_phs), cancel_alpha_phs = false; end

    ptype = char(string(ptype));
    ftype = char(string(ftype));

    % robust parsing to allow 'True'/'False' strings
    if ischar(cancel_alpha_phs) || isstring(cancel_alpha_phs)
        cancel_alpha_phs = strcmpi(string(cancel_alpha_phs), "true");
    end
    cancel_alpha_phs = logical(cancel_alpha_phs);

    % ---------------- enforce LS-only ----------------
    if ~strcmpi(ftype, 'ls')
        error(['sigpy_slr_dzrf: Only ftype="ls" is supported in this MATLAB clone. ', ...
               'Other ftypes ("ms","pm","min","max") were not robust after Python->MATLAB porting ', ...
               'due to solver/numerical differences (especially remez/min-phase factorization) ', ...
               'which can lead to non-reproducible RF waveforms.']);
    end

    % ---------------- ripple adjustment (sigpy calc_ripples) ----------------
    [bsf, d1m, d2m] = calc_ripples(ptype, d1, d2);

    % ---------------- LS beta design (sigpy dzls) ----------------
    b = dzls(n, tb, d1m, d2m);

    % ---------------- ptype handling (sigpy dzrf) ----------------
    if strcmpi(ptype, 'st')
        rf = b;
    elseif strcmpi(ptype, 'ex')
        b = bsf .* b;
        rf = b2rf(b, cancel_alpha_phs);
    else
        b = bsf .* b;
        rf = b2rf(b, false);
    end

    rf = rf(:); % column

    % =====================================================================
    %                          support functions
    % =====================================================================

    function [bsf_, d1_, d2_] = calc_ripples(ptype_, d1in, d2in)
        % Direct port of sigpy calc_ripples()
        switch lower(ptype_)
            case 'st'
                bsf_ = 1;
                d1_  = d1in;
                d2_  = d2in;
            case 'ex'
                bsf_ = sqrt(1/2);
                d1_  = sqrt(d1in/2);
                d2_  = d2in/sqrt(2);
            case 'se'
                bsf_ = 1;
                d1_  = d1in/4;
                d2_  = sqrt(d2in);
            case 'inv'
                bsf_ = 1;
                d1_  = d1in/8;
                d2_  = sqrt(d2in/2);
            case 'sat'
                bsf_ = sqrt(1/2);
                d1_  = d1in/2;
                d2_  = sqrt(d2in);
            otherwise
                error('Pulse type ("%s") is not recognized.', ptype_);
        end
    end

    function di = dinf_local(d1__, d2__)
        % sigpy.mri.rf.util.dinf (same coefficients used in sigpy)
        a1 = 5.309e-3;
        a2 = 7.114e-2;
        a3 = -4.761e-1;
        a4 = -2.66e-3;
        a5 = -5.941e-1;
        a6 = -4.278e-1;

        l10d1 = log10(d1__);
        l10d2 = log10(d2__);
        di = (a1*l10d1.^2 + a2*l10d1 + a3).*l10d2 + ...
             (a4*l10d1.^2 + a5*l10d1 + a6);
    end

    function h = dzls(n_, tb_, d1__, d2__)
        % Direct port of sigpy dzls(), including the explicit half-sample shift.
        di = dinf_local(d1__, d2__);
        w  = di / tb_;

        f = [0, (1-w)*(tb_/2), (1+w)*(tb_/2), (n_/2)];
        f = f / (n_/2);
        m = [1, 1, 0, 0];
        wt = [1, d1__/d2__];

        % SciPy: signal.firls(n+1,...) returns length n+1
        % MATLAB: firls(order=n,...) returns length n+1 -> matches
        h = firls(n_, f, m, wt);  % length n+1

        % sigpy: shift filter half a sample to make it symmetric like MATLAB
        k = [0:(n_/2), (-n_/2):-1]; % length n+1
        c = exp(1j * 2*pi / (2*(n_+1)) * k);

        % sigpy uses sp.fft/ifft with center=False -> plain FFT/IFFT
        h = real(ifft(fft(h(:).') .* c));

        % lop off extra sample: h[:n]
        h = h(1:n_).';
    end

    function rf_ = b2rf(b_, cancel_alpha_phs_)
        % Direct port of sigpy b2rf()
        a_ = b2a(b_);
        if cancel_alpha_phs_
            b_a_phase = fft(b_(:)) .* exp(-1j * angle(fft(flipud(a_(:)))));
            b_ = ifft(b_a_phase);
        end
        rf_ = ab2rf(a_, b_);
    end

    function a_ = b2a(b_)
        % Direct port of sigpy b2a()
        n_ = numel(b_);

        npad = n_ * 16;
        bcp = complex(zeros(npad,1));
        bcp(1:n_) = b_(:);

        bf = fft(bcp);
        bfmax = max(abs(bf));
        if bfmax >= 1
            bf = bf / (1e-7 + bfmax);
        end

        afa = mag2mp(sqrt(1 - abs(bf).^2));
        a_ = fft(afa) / npad;
        a_ = a_(1:n_);
        a_ = flipud(a_);
    end

    function a_ = mag2mp(x_)
        % Direct port of sigpy mag2mp()
        n_ = numel(x_);
        xl  = log(abs(x_(:)));
        xlf = fft(xl);

        xlfp = xlf;
        xlfp(1) = xlf(1);                 % DC
        xlfp(2:(n_/2)) = 2*xlf(2:(n_/2));  % double positive freqs
        xlfp(n_/2+1) = xlf(n_/2+1);        % keep Nyquist
        xlfp(n_/2+2:end) = 0;              % zero negative freqs

        xlaf = ifft(xlfp);
        a_ = exp(xlaf);
    end

    function rf_ = ab2rf(a_, b_)
        % Direct port of sigpy ab2rf()
        n_ = numel(a_);
        rf_ = complex(zeros(n_,1));

        a_ = complex(a_(:));
        b_ = complex(b_(:));

        for ii = n_:-1:1
            cj = sqrt(1 / (1 + abs(b_(ii)/a_(ii))^2));
            sj = conj(cj * b_(ii) / a_(ii));
            theta = atan2(abs(sj), cj);
            psi = angle(sj);
            rf_(ii) = 2 * theta * exp(1j * psi);

            if ii > 1
                at = cj*a_ + sj*b_;
                bt = -conj(sj)*a_ + cj*b_;
                a_ = at(2:ii);
                b_ = bt(1:ii-1);
            end
        end
    end

end
