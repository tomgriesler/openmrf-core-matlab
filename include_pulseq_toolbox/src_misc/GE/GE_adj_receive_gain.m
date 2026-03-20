function seq = GE_adj_receive_gain(system, NR, Trec, adc, alpha, dz, external_path, wip_id)

    % Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 14.03.2026

    % calculate sequence objects
    [rf, gz, gz_reph] = mr.makeSincPulse( alpha, ...
                                          system, ...
                                          'Duration', 2.5*1e-3,...
                                          'SliceThickness', dz * 1.25, ...
                                          'timeBwProduct', 4, ...
                                          'apodization', 0.5, ...
                                          'PhaseOffset', -pi/2, ...
                                          'maxSlew', system.maxSlew/sqrt(3), ...
			                              'use', 'excitation' );
    adc.phaseOffset = 0;
    [gx_crush, gy_crush, gz_crush] = CRUSH_x_y_z(2, 2, 4, 1e-3, 1e-3, dz, 1/sqrt(3), 1/sqrt(3), system);
    Trec = mr.makeDelay(round(Trec/system.blockDurationRaster)*system.blockDurationRaster);

    % built sequence for pre-scans
    seq  = mr.Sequence(system);
    for loop_rep = 1 : NR
        seq.addTRID('adj_receive');
        seq.addBlock(rf, gz);
        seq.addBlock(gz_reph);
        seq.addBlock(adc, mr.makeDelay(ceil(mr.calcDuration(adc)*1.05/system.blockDurationRaster)*system.blockDurationRaster));
        seq.addBlock(gx_crush, gy_crush, gz_crush);
        seq.addBlock(Trec);
    end

    % set definitions & export pre-scan .seq file
    idx      = strfind(external_path, '/');
    seq_name = [external_path(idx(end-1)+1:idx(end)-1) '_' external_path(idx(end-3)+1:idx(end-2)-1) '_adj_receive_gain.seq'];
    seq.setDefinition('Name',       seq_name);
    seq.setDefinition('Scan_ID',    wip_id);
    seq.setDefinition('FOV',        [0.3 0.3 dz]);
    seq.setDefinition('Rot_Matrix', eye(3));
    seq.write([external_path(1:idx(end)) seq_name]);

end
