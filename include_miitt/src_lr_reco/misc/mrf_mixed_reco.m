function reco = mrf_mixed_reco(DATA, ktraj, Nxy, fov)

    Ncoils     = size(DATA,1);
    N          = [Nxy, Nxy];
    kxy        = ktraj(:,:)';
    nufft_args = {N, [6 6], 2*N, N/2, 'table', 2^12, 'minmax:kb'};
    nufft_st   = nufft_init( kxy / max(abs(kxy(:)))*pi, N, [6 6], N*2, N/2, 'minmax:kb');
    G          = Gmri(kxy, true(N), 'fov', fov, 'basis', {'rect'}, 'nufft', nufft_args);
    dcf        = abs(mri_density_comp( kxy, 'pipe', 'G', G.arg.Gnufft)) *2^2 *Nxy^2;
    reco       = nufft_adj( permute(DATA(:,:),[2,1]) .* repmat(dcf,1,Ncoils), nufft_st );
    reco       = permute(reco, [3,1,2]);
    reco(end+1,:,:) = openadapt(reco);
    reco = reco / max(abs(reco(:)));
    xtv(reco)

end