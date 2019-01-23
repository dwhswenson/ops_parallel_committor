# Example: Alanine Dipeptide

The `ops_parallel_committor` script takes an OPS file as input. That file
must include named states, as well as snapshots that are outside of any
named state (in the transition region).

To generate that file, we'll run transition path sampling. The scripts in
this directory are based on the OPS examples for TPS of alanine dipeptide.
Run them will the following commands:

```bash
python AD_tps_1_trajectory.py
python AD_tps_2a_run_flex.py
```

The second script, in particular, will take some time. Once both scripts
have been run, you should have the file `alanine_dipeptide_tps.nc` in this
directory. That will be the file that you draw initial points from.

For a quick test, run

`python ../ops_parallel_committor.py --n-snapshots 2 --n-per-snapshot 2 -S alpha_R -S C_7eq alanine_dipeptide_tps.nc -o committor.nc`

For a longer test, you can increase the `--n-snapshots`, which will increase
the number of tasks, and the `--n-per-snapshot`, which will increase time
taken for each task.

For large numbers of snapshots, you may need to increase the number of steps
taken in the `AD_tps_2a_run_flex.py` script. Currently it is set to
`n_steps=500`. I would recommend using `n_snapshots <= n_steps * 5`.
