function wave_hash = pulseq_get_wave_hash(waveform)
    % convert the waveform to a string of bytes
    data = typecast(waveform(:), 'uint8');
    
    % use Java's MD5 hashing
    md = java.security.MessageDigest.getInstance('MD5');
    md.update(data);
    hashBytes = typecast(md.digest(), 'uint8');
    
    % convert the hash bytes to a hexadecimal string
    wave_hash = dec2hex(hashBytes)';
    wave_hash = lower(wave_hash(:)');
end

