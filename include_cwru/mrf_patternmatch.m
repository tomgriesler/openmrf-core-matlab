function match = mrf_patternmatch(images, dict, look_up, look_up_names, roimask, numBlocks)

    % ----- input: -----
    % images:        Ny x Nx x TRs MRF images
    % dict:          TRs x tissuePropertyValues the MRF dictionary
    % look_up:       tissuePropertyValues the list of tissue property combinations in the dictionary
    % look_up_names: names of the tissue parameters of the look_up table
    % roimask:       Ny x Nx a binary mask
    % numBlocks:     if greater than 1, then we divide the images into smaller blocks rather than matchign the entire image at once.
     
    % ----- output: -----
    % match: structured object
    %  .IND: global index map for best matches
    % .IP:   inner product indicating the goodness of fit
    % .M0:   M0 map
    % .'X':  parameter maps depending on look_up table
    
    [ny,nx,~] = size(images);
    
    if isempty(roimask)
        roimask = true(ny,nx);
    else
        roimask = (roimask==1);
    end
    if numBlocks < 1
        numBlocks = 1;
    end
    
    images   = permute(images,[3 1 2]); %
    idxROI   = find(roimask==1);
    n        = numel(idxROI);
    images   = images(:,idxROI);
    pixels   = ceil(n/numBlocks);
    dictnorm = sqrt(sum(dict.*conj(dict)));
    
    IND = zeros(ny,nx, 'single');
    IP  = zeros(ny,nx, 'single');
    M0  = zeros(ny,nx, 'single');
    
    for iblock=1:numBlocks
        index = (iblock-1)*pixels+1:iblock*pixels;
        index(index > n)=[];
        Nindex = length(index);
        if isempty(index)
            continue;
        end
        
        xx = squeeze(images(:,index)).';
        xxNorm = sqrt(sum(xx.*conj(xx),2));
        normAll = xxNorm*dictnorm;
        innerProduct = conj(xx)*dict./normAll;
        
        [values, indexm] = max(abs(innerProduct),[],2);
        dictCol = dict(:,indexm);
        
        coefSave = complex(zeros(1,Nindex),zeros(1,Nindex));
        for iindex = 1:Nindex
            coefSave(iindex) = pinv(dictCol(:,iindex))*xx(iindex,:).';
        end
        
        IND(idxROI(index)) = indexm;
        IP(idxROI(index))  = values(:);
        M0(idxROI(index))  = coefSave;
    end
    
    match.IND = IND;
    match.IP  = IP;
    match.M0  = M0;
    IND(IND==0) = 1;
    for j=1:numel(look_up_names)
        temp_map = look_up(:,j);
        temp_map = temp_map(IND) .* roimask;
        eval(['match.' look_up_names{j} '=temp_map;']);
    end

end

