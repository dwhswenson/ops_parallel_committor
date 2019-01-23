# Parallel Committor for OpenPathSampling

This is a hacky way to get a parallel committor in OPS. It creates many
tasks, each represented by a specific input file. Then the individual input
files can be run as an OPS task, each in a separate process. Finally, the
results are returns and made into the tools necessary for further analysis
with OPS.
