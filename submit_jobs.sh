#! /bin/bash

queue=${queue:-test}

epochs_max=50000
epochs_step=10000

learning_rate_lambda_step=1000

radius_lambda_step=1000

num_procs=${num_procs:-4}
num_procs_step=4
num_procs_max=${num_procs_max:-256}


while [ $num_procs -le $num_procs_max ]
do

    epochs=10000
    learning_rate_lambda=1000
    radius_lambda=1000

    while [ $epochs -le $epochs_max ]
    do

        echo sbatch -n $num_procs -o stdout_${epochs}_epochs_${num_procs}_processes.txt -p $queue ompi kohonen_learn -i test_data/experiment_data.txt -x 40 -y 40 -m 3 -e $epochs -a0 0.1 -al $learning_rate_lambda -r0 40.0 -rl $radius_lambda -o neuron_weights_${epochs}_epochs_${num_procs}_processes.txt

        epochs=$[$epochs+$epochs_step]
        learning_rate_lambda=$[$learning_rate_lambda+$learning_rate_lambda_step]
        radius_lambda=$[$radius_lambda+$radius_lambda_step]
    done

    num_procs=$[$num_procs+$num_procs_step]
done
