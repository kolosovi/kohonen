#! /bin/bash

queue=${queue:-test}

neurons_start=${neurons_start:-1000}
neurons_max=${neurons_max:-10000}
neurons_step=${neurons_step:-1000}

learning_rate_lambda_step=1000

radius_lambda_step=1000

num_procs=${num_procs:-1}
num_procs_step=${num_procs_step:-1}
num_procs_max=${num_procs_max:-8}

used_slots=0
max_slots=3


QUEUE_CHECK_INTERVAL=60

while [ $num_procs -le $num_procs_max ]
do

    neurons=${neurons_start}
    learning_rate_lambda=1000
    radius_lambda=1000

    while [ $neurons -le $neurons_max ]
    do
        used_slots=$(squeue -h -u `whoami` | wc -l)

        # Stop until there are slots for another job
        while [ $used_slots -ge $max_slots ]
        do
            sleep $QUEUE_CHECK_INTERVAL
            used_slots=$(squeue -h -u `whoami` | wc -l)
        done

        sbatch -t 15 -n $num_procs -o stdout_${neurons}_neurons_${num_procs}_processes.txt -p $queue ompi -n $num_procs kohonen_learn -i test_data/experiment_data.txt -x 1 -y $neurons -m 3 -e 10000 -a0 0.1 -al $learning_rate_lambda -r0 $neurons -rl $radius_lambda -o neuron_weights_${neurons}_neurons_${num_procs}_processes.txt

        neurons=$[$neurons+$neurons_step]
    done

    num_procs=$[$num_procs+$num_procs_step]
done
