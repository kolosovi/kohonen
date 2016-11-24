#! /bin/bash

queue=${queue:-test}

epochs_max=50000
epochs_step=10000

learning_rate_lambda_step=1000

radius_lambda_step=1000

num_procs=${num_procs:-8}
num_procs_step=${num_procs_step:-1}
num_procs_max=${num_procs_max:-256}

used_slots=0
max_slots=3


QUEUE_CHECK_INTERVAL=60


while [ $num_procs -le $num_procs_max ]
do

    epochs=10000
    learning_rate_lambda=1000
    radius_lambda=1000

    while [ $epochs -le $epochs_max ]
    do
        used_slots=$(squeue | grep 'zackw' | wc -l)

        # Stop until there are slots for another job
        while [ $used_slots -ge $max_slots ]
        do
            sleep $QUEUE_CHECK_INTERVAL
            used_slots=$(squeue | grep 'zackw' | wc -l)
        done

        sbatch -t 15 -n $num_procs -o stdout_${epochs}_epochs_${num_procs}_processes.txt -p $queue ompi -n $num_procs kohonen_learn -i test_data/experiment_data.txt -x 40 -y 40 -m 3 -e $epochs -a0 0.1 -al $learning_rate_lambda -r0 40.0 -rl $radius_lambda -o neuron_weights_${epochs}_epochs_${num_procs}_processes.txt

        epochs=$[$epochs+$epochs_step]
        learning_rate_lambda=$[$learning_rate_lambda+$learning_rate_lambda_step]
        radius_lambda=$[$radius_lambda+$radius_lambda_step]
    done

    num_procs=$[$num_procs+$num_procs_step]
done
