function [T1, T2] = NIST_references(ref_name)

% ----- input: -----
% ref_name: e.g. '1.5T_NiCl2' '1.5T_MnCl2' '3.0T_NiCl2' '3.0T_MnCl2'
% or a result .mat file name created via reco_spin_echo_T1_T2

% ----- output: -----
% T1: 14x1 array [s]
% T2: 14x1 array [s]

% Stupic KF, et al.
% A standard system phantom for magnetic resonance imaging.
% Magn Reson Med. 2021 Sep;86(3):1194-1211.
% doi: 10.1002/mrm.28779.
% onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fmrm.28779&file=mrm28779-sup-0001-Supinfo.pdf
switch ref_name
    case '1.5T_NiCl2'
        T1 = [22.859  32.45  46.32  65.40  92.70  132.7  184.6  260.1  367.9  514.1  730.8  1012  1489  2033]' *1e-3;
        T2 = [20.2  28.7  41.0  57.7  81.9  117.1  161.9  227.7  321.2  446.3  628.5  859.3  1244  1669]' *1e-3;
    case '1.5T_MnCl2'
        T1 = [83.33  106.4  160.2  194.9  292.9  413.4  550.2  752.2  1030  1237  1539  1870  2183  2736]' *1e-3;
        T2 = [8.15  10.47  15.99  19.76  30.62  45.28  64.84  91.76  140.6  184.9  267.0  416.5  594.3  939.4]' *1e-3;
    case '3.0T_NiCl2'
        T1 = [21.719  30.84  44.53  62.70  89.00  125.9  175.3  247.13  351.5  496.7  706.0  984.1  1454  1989]' *1e-3;
        T2 = [15.83  22.38  31.86  45.70  64.30  90.30  127.3  180.8  255.5  359.6  510.1  717.9  1077  1465]' *1e-3;
    case '3.0T_MnCl2'
        T1 = [90.90  126.9  176.9  244.2  336.5  458.4  608.6  801.7  1044  1332  1604  1907  2173  2480]' *1e-3;
        T2 = [5.59  7.91  11.24  15.81  22.56  31.97  46.42  64.07  96.89  133.27  190.94  278.1  403.5  581.3]' *1e-3;
    otherwise
        try
            load([ref_name '.mat'], 'T1', 'T2'); % created via reco_spin_echo_T1_T2
        catch
            error(['no reference data found for: ' ref_name]);
        end
end

end
