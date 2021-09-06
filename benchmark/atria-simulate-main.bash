#! bash

for i in 16 20 24 28 33
do
    echo "Start: adapter length $i"
    echo "Start: adapter length $i"
    echo "Start: adapter length $i"
    bash $atria/benchmark/atria-simulate.bash $i
done

working_dir=~/analysis/atria-benchmark/atria_simulate

cd $working_dir

atria statplot -i auto -l DIR2
