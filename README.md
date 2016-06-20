# Photon-efficient imaging with a single-photon camera

Repo includes data and MATLAB code used in <br />
"Photon-efficient imaging with a single-photon camera" <br />
by D.Shin&#42;, F.Xu&#42;, D.Venkatraman, R.Lussana, F.Villa, F.Zappa, V.K.Goyal, F.N.C.Wong, and J.H.Shapiro <br />
in Nature Communications

The given MATLAB implementation uses photon arrival data obtained by a SPAD camera
to compute accurate 3D and reflectivity of room-scale objects as proof-of-concept. The results are compared
to those from the baseline filtered histogram approach.

The computational imaging framework uses variants of 
SPIRAL-TAP software by Harmany et al.
(http://drz.ac/code/spiraltap/)
and 
OMP implementation by Becker
(http://www.mathworks.com/matlabcentral/fileexchange/32402-cosamp-and-omp-for-sparse-recovery).
