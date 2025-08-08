#!/bin/bash
cat <<EOF > parameters_L_$1_tau_$2_$3.dat

&parametersfile
inCluster = T
L = $1,$1,$1
N_thermalization = 1000
N_measurements = 5000,
N_skip = 10,
isbeta = F
beta_i = 2.0,
beta_f = 0.1,
n_beta = 41,
algorithm = "$3"
Nhmc=20
Thmc = 1.0
start = "cold"
equilibrium=F
tau_Q = $2
/

EOF
