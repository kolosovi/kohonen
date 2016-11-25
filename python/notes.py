import re
name_re = re.compile(r"stdout_(?P<num_epochs>\d+)_(epochs|neurons)_(?P<num_procs>\d+)_processes.txt")
import os
filenames = [name for name in os.listdir(".") if "stdout" in name]

all_neurons = range(1000, 11000, 1000)
all_procs = range(9, 49)
all_combinations = [(neurons, procs) for neurons in all_neurons for procs in all_procs]
retrieved_combinations = []

for filename in filenames:
    match = re.search(name_re, filename)
    if not match:
            continue
    num_n = int(match.groupdict()["num_epochs"])
    num_p = int(match.groupdict()["num_procs"])
    retrieved_combinations.append((num_n, num_p))

combinations_to_do = set(all_combinations) - set(retrieved_combinations)

def num_with_procs(seq, procs):
    return len([i for i in seq if i[1] == procs])

[i for i in xrange(9, 49) if num_with_procs(combinations_to_do, i) == 10]

invalid_filenames = []

for name in filenames:
    with open(name, "r") as stdout_file:
        try:
            execution_time = float(stdout_file.read().strip())
        except:
            invalid_filenames.append(name)
