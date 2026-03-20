function [gx_crush, gy_crush, gz_crush] = CRUSH_x_y_z(crush_nTwists_x, crush_nTwists_y, crush_nTwists_z, dx, dy, dz, lim_grad, lim_slew, system)

    % Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

    crush_area_x   = crush_nTwists_x / dx;
    gx_crush       = mr.makeTrapezoid('x', 'Area', crush_area_x, 'maxGrad', system.maxGrad*lim_grad, 'maxSlew', system.maxSlew*lim_slew, 'system', system);

    crush_area_y   = crush_nTwists_y / dy;
    gy_crush       = mr.makeTrapezoid('y', 'Area', crush_area_y, 'maxGrad', system.maxGrad*lim_grad, 'maxSlew', system.maxSlew*lim_slew, 'system', system);
    gy_crush.delay = ceil(gx_crush.riseTime*1.05 / system.gradRasterTime) * system.gradRasterTime; 

    crush_area_z   = crush_nTwists_z / dz;
    gz_crush       = mr.makeTrapezoid('z', 'Area', crush_area_z, 'maxGrad', system.maxGrad*lim_grad, 'maxSlew', system.maxSlew*lim_slew, 'system', system);
    gz_crush.delay = ceil(gx_crush.riseTime*1.05 / system.gradRasterTime) * system.gradRasterTime + ...
                     ceil(gy_crush.riseTime*1.05 / system.gradRasterTime) * system.gradRasterTime;  

end
