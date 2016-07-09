# Photon-efficient imaging with a single-photon camera

Repo includes data and MATLAB code used in <br />
"Photon-efficient imaging with a single-photon camera" <br />
by D.Shin&#42;, F.Xu&#42;, D.Venkatraman, R.Lussana, F.Villa, F.Zappa, V.K.Goyal, F.N.C.Wong, and J.H.Shapiro <br />
published in Nature Communications

[Link to paper](http://www.nature.com/ncomms/2016/160624/ncomms12046/full/ncomms12046.html)
<br>
[Link to project page](https://photon-efficient-imaging.github.io/single-photon-camera-project/)
<br>

How to cite (BibTeX): <br />
&nbsp; @article{shin2016photon,<br />
&nbsp;   title={Photon-efficient imaging with a single-photon camera},<br />
&nbsp;   author={Shin, Dongeek and Xu, Feihu and Venkatraman, Dheera and Lussana, Rudi and Villa, Federica and <br />
&nbsp;   Zappa, Franco and Goyal, Vivek K and Wong, Franco NC and Shapiro, Jeffrey H},<br />
&nbsp;   journal={Nature Communications},<br />
&nbsp;   volume={7},<br />
&nbsp;   year={2016},<br />
&nbsp;   publisher={Nature Publishing Group}<br />
&nbsp; }

The given MATLAB implementation uses photon arrival data obtained by a SPAD camera
to compute accurate 3D and reflectivity of room-scale objects as proof-of-concept. The results are compared
to those from the baseline filtered histogram approach.

The computational imaging framework uses variants of 
SPIRAL-TAP software by Harmany et al.
(http://drz.ac/code/spiraltap/)
and 
OMP implementation by Becker
(http://www.mathworks.com/matlabcentral/fileexchange/32402-cosamp-and-omp-for-sparse-recovery).
