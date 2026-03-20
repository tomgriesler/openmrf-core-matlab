%%
clear
load('temp.mat')

t1lims  = [0 2];
t2lims  = [0 0.12];
t1plims = [0 0.18];

flag_log = 1;

[t1cmp,  T1_Map_clip]  = get_cmp('T1',      1000, flag_log, t1lims,  T1_Map);
[t2cmp,  T2_Map_clip]  = get_cmp('T2',      1000, flag_log, t2lims,  T2_Map);
[t1pcmp, T1p_Map_clip] = get_cmp('inferno', 1000, flag_log, t1plims, T1p_Map);

figure()
ax1 = subplot(2,3,1);
imagesc(T1_Map, t1lims); axis image; axis off; colormap(ax1, t1cmp); colorbar;

ax2 = subplot(2,3,2);
imagesc(T2_Map, t2lims); axis image; axis off; colormap(ax2, t2cmp); colorbar;

ax3 = subplot(2,3,3);
imagesc(T1p_Map, t1plims); axis image; axis off; colormap(ax3, t1pcmp); colorbar;

ax4 = subplot(2,3,4);
imagesc(T1_Map_clip, t1lims); axis image; axis off; colormap(ax4, t1cmp); colorbar;

ax5 = subplot(2,3,5);
imagesc(T2_Map_clip, t2lims); axis image; axis off; colormap(ax5, t2cmp); colorbar;

ax6 = subplot(2,3,6);
imagesc(T1p_Map_clip, t1plims); axis image; axis off; colormap(ax6, t1pcmp); colorbar;
