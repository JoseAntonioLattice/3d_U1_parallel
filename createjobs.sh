#!/bin/bash


filename=3d_U1_parallel_2_$1x$2x$3
jobname=3dU1_2_$1x$2x$3-$7
#touch $filename$1.slurm

nodes=$8

cat <<EOF > jobs/$filename-$(($1*$2*$3)).slurm
#!/bin/bash
#SBATCH --job-name=logs/$jobname 			# Nombre del trabajo
#SBATCH --output=logs/$jobname.log   		# Archivo de registro de salida
#SBATCH --error=logs/$jobname.err    		# Archivo de registro de errores
#SBATCH --partition=QuantPhysMC   	# Nombre de la partición o cola de trabajos
#SBATCH --nodes=$nodes			# Número de nodos a utilizar (puedes cambiarlo)
#SBATCH --ntasks-per-node=$(($1*$2*$3/$nodes))		# Número de tareas por nodo (1 para ejecución serial)
#SBATCH --cpus-per-task=1 		# Número de CPUs por tarea (puedes cambiarlo)
#SBATCH --mem=4G      			# Memoria RAM necesaria (puedes cambiarlo)

module load lamod/coarrays/2.10 
cd ~/3d_U1_parallel
# Comando para ejecutar tu programa
for i in 8 16 32 64 128 256 512; do { echo $1 $2 $3; echo input/parameters_L_$4_tau_\$i_$6.dat; echo measurements_$7.dat; } | LD_LIBRARY_PATH=$HOME/Fortran/lib cafrun -n $(( $1*$2*$3 )) build/3d_U1_parallel; done

EOF
