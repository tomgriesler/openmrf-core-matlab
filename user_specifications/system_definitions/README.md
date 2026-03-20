# Contents: /system_definitions

This folder contains `.csv` files that define hardware specifications of different MRI scanners. These configuration files are used to ensure compliance with scanner-specific constraints during sequence compilation.

## Example Format

```
B0;1.5
gamma;4.257747851783256e+07
maxGrad;43
gradUnit;mT/m
maxSlew;180
slewUnit;T/m/s
maxB1;20
rfDeadTime;1e-4
rfRingdownTime;1e-4
adcDeadTime;1e-5
adcRasterTime;1e-7
rfRasterTime;1e-6
gradRasterTime;1e-5
blockDurationRaster;1e-5
ascfile;Sola_MP_GPA_K2368_2250V_950A_GC04XQ.asc
```

### Notes

- The **ascfile** field is optional. It refers to a gradient system calibration file used for PNS (Peripheral Nerve Stimulation) simulation.
- Redistribution of `.asc` files is **not permitted** due to vendor restrictions.

_Maximilian Gram: 11.07.2025_