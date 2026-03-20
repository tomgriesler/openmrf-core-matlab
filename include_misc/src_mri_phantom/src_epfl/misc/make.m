cd misc;

insidepoly_install

archstr = computer('arch');
switch archstr
    case 'glnxa64'
        archstr = 'a64 ';
end
include = ' ';% -I/opt/local/var/macports/software/gcc46/4.6-20110312_0/opt/local/lib/gcc46/gcc/x86_64-apple-darwin10/4.6.0/include ';
opts = {['arch=' archstr],'CXX=g++-4.2','LDCXX=g++-4.2','CXXOPTIMFLAGS="-O3 -DNDEBUG"','LDOPTIMFLAGS="-O3"','-v'};%' -v CC=gcc-4.4 CXX=g++-4.4 -fopenmp -lpthread -O3 -march=x86_64 -msse3 -mssse3 ';%'-O3 -march=core2 ';% '-O3  ';
ext = {'.c','.cpp'};
list = {'myerfzparts_c','biot_savart_map_c'}; %,'myerfzparts_dev_c','erf_c'
for i = 1:length(list)
    for j = 1:length(ext)
        if exist([list{i} ext{j}],'file')%&&~exist([list{i} '.mex' archstr],'file')
            fprintf('Compiling %s\n',[list{i} ext{j}])
            mex(opts{:},which([list{i} ext{j}])); %['arch=' archstr include opts]
        end
    end
end

cd ..;
clear archstr ext i j list opts;