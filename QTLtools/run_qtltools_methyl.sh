for CHUNK in $(seq 1 1000);
do
  sbatch --export=ALL,CHUNK=$CHUNK qtltools_methyl.sh

  sleep 5
done

