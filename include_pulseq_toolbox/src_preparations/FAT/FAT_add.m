% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% add fat suppression
if strcmp(FAT.mode, 'on')
    seq.addBlock(FAT.rf);
    seq.addBlock(FAT.gz_crush);
end