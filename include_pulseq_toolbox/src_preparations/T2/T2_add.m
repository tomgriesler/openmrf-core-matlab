% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% temporary T2 objects
d1           = T2.T2_objs.d1;
t_inter      = T2.T2_objs.t_inter;
rf_90_td     = T2.T2_objs.rf_90_td;
rf_90_tu     = T2.T2_objs.rf_90_tu;
rf_comp_pos  = T2.T2_objs.rf_comp_pos;
rf_comp_neg  = T2.T2_objs.rf_comp_neg;
gx_crush     = T2.T2_objs.gx_crush;
gy_crush     = T2.T2_objs.gy_crush;
gz_crush     = T2.T2_objs.gz_crush;

%% add T2 T2 and T2p preparation
seq.addBlock(rf_90_td);
seq.addBlock(t_inter(loop_T2));
seq.addBlock(rf_comp_pos);
seq.addBlock(t_inter(loop_T2));
seq.addBlock(t_inter(loop_T2));
seq.addBlock(rf_comp_neg);
seq.addBlock(t_inter(loop_T2));
seq.addBlock(rf_90_tu);
seq.addBlock(d1);

%% add crusher gradients
seq.addBlock(gx_crush, gy_crush, gz_crush);
seq.addBlock(d1);

%% clear temp objects
clear d1 t_inter rf_90_td rf_90_tu rf_comp_pos rf_comp_neg gx_crush gy_crush gz_crush cycling_list;