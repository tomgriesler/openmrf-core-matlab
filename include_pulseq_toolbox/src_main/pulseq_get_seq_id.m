function [seq_id, scan_id, wip_id] = pulseq_get_seq_id()
% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026
    temp1 = year(datetime('now'));
    temp2 = month(datetime('now'));
    temp3 = day(datetime('now'));
    temp4 = hour(datetime('now'));
    temp5 = minute(datetime('now'));
    
    temp1 = num2str(temp1);
    temp1(1:2) = [];
    
    if temp2<10
        temp2 = ['0' num2str(temp2)];
    else
        temp2 = num2str(temp2);
    end
    
    if temp3<10
        temp3 = ['0' num2str(temp3)];
    else
        temp3 = num2str(temp3);
    end
    
    if temp4<10
        temp4 = ['0' num2str(temp4)];
    else
        temp4 = num2str(temp4);
    end
    
    if temp5<10
        temp5 = ['0' num2str(temp5)];
    else
        temp5 = num2str(temp5);
    end
    
    seq_id  = [ temp1 temp2 temp3 '_' temp4 temp5 ];
    scan_id = int64(str2num([ temp1 temp2 temp3 temp4 temp5 ]));
    wip_id  = int64(str2num([ temp1(2) temp2 temp3 temp4 temp5 ]));
    % wip_id is used in seq.setDefinition() and saved in the sWipMemBlock of meas.dat
    % here only 9 decimal places can be used, so we remove the 2 of the year 2X
    % works fine until 2030 :D

end

