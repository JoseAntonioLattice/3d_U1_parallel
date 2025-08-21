#!/bin/bash
cat <<EOF > parameters_L_$1_tau_$2_$3.dat

&parametersfile
inCluster = T
L = $1,$1,$1
save_thermalized_conf = T
N_thermalization = 10000
N_measurements = 5000,
N_skip = 10,
isbeta = F
readbeta = T
beta_i = 2.0,
beta_f = 0.1,
n_beta = 41,
algorithm = "$3"
Nhmc=20
Thmc = 1.0
start = "cold"
equilibrium=T
tau_Q = $2
savelastconf=F
readconfig=F
/

EOF
