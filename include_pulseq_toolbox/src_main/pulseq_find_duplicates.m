function pulseq_find_duplicates(mypath)
% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026
    mypath = fileparts(mypath);
    P = path;
    P = strsplit(P, pathsep());
    P = P(strncmpi(mypath,P,length(mypath)));
    P = cellfun(@(x) what(x),P,'UniformOutput',false);
    P = vertcat(P{:});
    Q = arrayfun(@(x) x.m,P,'UniformOutput',false);
    Q = vertcat(Q{:});
    R = arrayfun(@(x) repmat({x.path},size(x.m)),P,'UniformOutput',false);
    R = vertcat(R{:});
    [C,~,~] =unique(Q);
    for c=1:numel(C)
       ind=strcmpi(C{c},Q);
       if sum(ind)>1
           fprintf('duplicate %s at paths\n\t',C{c});
           fprintf('%s\n\t',R{ind});
           fprintf('\n');
       end
    end
end