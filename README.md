# FDTD-in-CEM
  The matlab codes show Finite Difference Time Domain (FDTD) method applied in electromagnetic problem

  FDTD method can calculate the time response in a region when some stimulation is introduced by applying the Maxwell equations. In this situation, Gaussian impulse is introduced. Frequency response can be determined by a Fourier transformation to time response. The disadvantage of this method is large time-comsuming and memory comsuming.
  A 2-D problem would be raised and the absorb boundary would be treated as MOR absorb boundary which can only absorb the electromagnetic wave vertically insert.
  The other file named "FDTD_obstacle" has introduced a metal block in a transmission line. The transmission line can be treated as a filter when the block is introduced.
