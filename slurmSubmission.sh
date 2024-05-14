module load python

if [ ! -d slurmOut/ ];
then
				mkdir slurmOut
fi

if [ -d .snakemake/lock ];
then
        echo "Pipeline is locked. Unlocking..."
        snakemake --unlock --use-envmodules --snakefile ./src/Snakefile.smk --cluster-config ./src/slurmConfig.json -R all --latency-wait 60 --cluster "sbatch -J {rule} -o slurmOut/slurm-%j.out -e slurmOut/slurm-%j.err -N1 -n {cluster.threads} --time {cluster.time} --mem={cluster.mem}" --jobs 100
fi

snakemake --use-envmodules --snakefile ./src/Snakefile.smk --cluster-config ./src/slurmConfig.json -R all --latency-wait 60 --cluster "sbatch -J {rule} -o slurmOut/slurm-%j.out -e slurmOut/slurm-%j.err -N1 -n {cluster.threads} --time {cluster.time} --mem={cluster.mem}" --jobs 100
