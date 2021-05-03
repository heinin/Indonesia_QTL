for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ;
do
  sbatch --export=ALL,CHR=$CHR qtltools_exp.sh

  sleep 5
done

