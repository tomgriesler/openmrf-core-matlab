% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% temporary MLEV objects
d1           = MLEV.MLEV_objs.d1;
t_inter_1    = MLEV.MLEV_objs.t_inter_1;
t_inter_2    = MLEV.MLEV_objs.t_inter_2;
rf_90_td     = MLEV.MLEV_objs.rf_90_td;
rf_90_tu     = MLEV.MLEV_objs.rf_90_tu;
rf_comp_pos  = MLEV.MLEV_objs.rf_comp_pos;
rf_comp_neg  = MLEV.MLEV_objs.rf_comp_neg;
gx_crush     = MLEV.MLEV_objs.gx_crush;
gy_crush     = MLEV.MLEV_objs.gy_crush;
gz_crush     = MLEV.MLEV_objs.gz_crush;
cycling_list = MLEV.MLEV_objs.cycling_list';
cycling_list = cycling_list(:);

%% add MLEV T2 and T2p preparation
seq.addBlock(rf_90_td);
seq.addBlock(t_inter_1);

for temp_loop = 1 : MLEV.n_composite(loop_MLEV)

    if cycling_list(temp_loop) == 1
        seq.addBlock(rf_comp_pos);
    end
    if cycling_list(temp_loop) == -1
        seq.addBlock(rf_comp_neg);
    end
    if temp_loop < MLEV.n_composite(loop_MLEV)
        seq.addBlock(t_inter_2);
    end

end
clear temp_loop;

seq.addBlock(t_inter_1);
seq.addBlock(rf_90_tu);
seq.addBlock(d1);

%% add crusher gradients
seq.addBlock(gx_crush, gy_crush, gz_crush);
seq.addBlock(d1);

%% clear temp objects
clear d1 t_inter_1 t_inter_2 rf_90_td rf_90_tu rf_comp_pos rf_comp_neg gx_crush gy_crush gz_crush cycling_list;