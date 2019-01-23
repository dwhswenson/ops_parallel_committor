
# coding: utf-8

# This is file runs the main calculation for the flexible length TPS simulation. It requires the file `alanine_dipeptide_tps_equil.nc`, which is written in the notebook `alanine_dipeptide_tps_first_traj.ipynb`.
# 
# In this file, you will learn:
# * how to set up and run a flexible length TPS simulation
# 
# NB: This is a long calculation. In practice, it would be best to export the Python from this notebook, remove the `live_visualizer`, and run non-interactively on a computing node.

# In[1]:
n_steps = 500


import openpathsampling as paths


# ## Load engine, trajectory, and states from file

# In[2]:


old_storage = paths.Storage("alanine_dipeptide_tps_equil.nc", "r")


# In[3]:


engine = old_storage.engines['300K']
C_7eq = old_storage.volumes.find('C_7eq')
alpha_R = old_storage.volumes.find('alpha_R')
traj = old_storage.samplesets[len(old_storage.samplesets)-1][0].trajectory
phi = old_storage.cvs.find('phi')
psi = old_storage.cvs.find('psi')
template = old_storage.snapshots[0]


# In[4]:


print(engine.name)


# ## TPS
# 
# As always, the process for setting up a simulation is:
# 
# 1. Create a `network`
# 2. Create a `move_scheme`
# 3. Set up `initial_conditions`
# 4. Create the `PathSampling` object and run it.
# 
# Since we already created all the input to these when we set up the first trajectory, we can load use the versions we loaded above.

# In[5]:


network = paths.TPSNetwork(C_7eq, alpha_R)


# In[6]:


scheme = paths.OneWayShootingMoveScheme(network, selector=paths.UniformSelector(), engine=engine)


# In[7]:


initial_conditions = scheme.initial_conditions_from_trajectories(traj)


# In[8]:


storage = paths.Storage("alanine_dipeptide_tps.nc", "w", template)
sampler = paths.PathSampling(storage=storage,
                             move_scheme=scheme,
                             sample_set=initial_conditions)


# Note: 10000 steps will take a long time. If you just want to run a little bit, reduce this number.

# In[9]:


#sampler.live_visualizer = paths.StepVisualization2D(network, phi, psi, [-3.14, 3.14], [-3.14, 3.14])
sampler.run(n_steps)


# With this done, you can go on to do the flexible-length parts of the analysis in `alanine_dipeptide_tps_analysis.ipynb`.
