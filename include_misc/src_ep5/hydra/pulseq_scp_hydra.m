function [out] = pulseq_scp_hydra(pulseq_user, external_path)

    temp_file = strrep(which('pulseq_scp_hydra'),'\','/');
    temp_file = fopen([temp_file(1:end-19), '/hydra.txt']);
    temp_cmd  = ['pscp -pw "' char(str2num(fgetl(temp_file))-ceil(exp(pi))) '" -sftp ' external_path ' ' pulseq_user '@hydra:/skyra/ep5/pulseq_' pulseq_user '/'];
    system(temp_cmd);
    fclose(temp_file);
    clear temp_file temp_line1 temp_line2 temp_cmd;
    out = 0;

end

