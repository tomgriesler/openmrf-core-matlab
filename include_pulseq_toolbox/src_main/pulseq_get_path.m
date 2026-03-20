function fun_path = pulseq_get_path(fun_name)
% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026
    temp = strsplit(which(fun_name), filesep);
    for j=1:numel(temp)-1
        if j==1
            fun_path = [ temp{j} '/'];
        else
            fun_path = [ fun_path temp{j} '/'];
        end
    end
end