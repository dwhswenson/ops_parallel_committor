from __future__ import print_function
import os
import argparse
import collections
import random
import warnings

import openpathsampling as paths
from tqdm import tqdm


def create_randomizer(engine):
    from simtk import unit as u
    # only support OpenMM so far
    # this might raise an AttributeError if using an engine or
    # integrator that isn't supported
    # integrator = engine.integrator
    # temperature = integrator.temperature
    temperature = 300 * u.kelvin  # TODO: hard coding this for now
    beta = 1.0 / (temperature * u.BOLTZMANN_CONSTANT_kB)

    randomizer = paths.RandomVelocities(beta=beta, engine=engine)
    return randomizer

def file_snapshots_committor(task_file_name, n_per_snapshot):
    """Run the committor for all snapshots in an input file.

    Returns
    -------
    list
        tuples of ``(uuid, state_name)``, where ``uuid`` is the UUID of the
        original snapshot and ``state_name`` is the name of the final volume
        the shot landed in
    """
    tmp_filename = "tmp_out_" + task_file_name
    task_file = paths.Storage(task_file_name, mode='r')
    states = task_file.tags['states']
    initial_snapshots=task_file.tags['initial_snapshots']
    storage = paths.Storage(tmp_filename, mode='w')
    sim = paths.CommittorSimulation(
        storage=storage,
        engine=task_file.tags['engine'],
        states=states,
        randomizer=task_file.tags['randomizer'],
        initial_snapshots=initial_snapshots
    )
    sim.output_stream.write("Running file: " + task_file_name + "\n")
    sim.run(n_per_snapshot)

    # now we convert the results to a format that can be returned over dask
    analyzer = paths.ShootingPointAnalysis(steps=None, states=states)

    key_to_init_uuid = {analyzer.hash_function(snap): snap.__uuid__
                   for snap in initial_snapshots}

    def final_state_name(step):
        state_list = analyzer.analyze_single_step(step)
        if len(state_list) > 1:
            raise RuntimeError("State definitions are overlapping!")
        return state_list[0].name

    def initial_snapshot_uuid(step):
        key = analyzer.step_key(step)
        return key_to_init_uuid[analyzer.hash_function(key)]

    # return pairs (initial_snapshot_uuid, final_state_name)
    # this is and (int, str) tuple for each step; should be trivial to
    # communicate back
    return_value = [(initial_snapshot_uuid(step), final_state_name(step))
                    for step in storage.steps]

    # cleanup
    task_file.close()
    storage.close()
    os.remove(tmp_filename)

    return return_value


class ParallelCommittorSimulation(object):
    """
    Parameters
    ----------
    tps_file : str
        filename for file from which initial snapshots will come
    state_names : list
        name of initial state (must be saved in ``tps_file``)
    n_snapshots : int
        number of snapshots to draw from the ``tps_file``
    n_per_snapshot : int
        number of shots to try per snapshot
    n_per_task : int
        currently unused (will define number of shots in a task)
    """
    def __init__(self, tps_file, state_names, n_snapshots, n_per_snapshot,
                 n_per_task=None, output=None):
        self.n_snapshots = n_snapshots
        self.n_per_snapshot = n_per_snapshot
        self.n_per_task = n_per_snapshot if None else n_per_task
        self.state_names = state_names
        self.tps_file = tps_file
        self.output = output
        self.states = []
        if self.tps_file is not None:
            storage = paths.Storage(self.tps_file, mode='r')
            self.states = [storage.volumes[name]
                           for name in self.state_names]
            storage.close()

    def sequential_main(self):
        tps_storage = paths.Storage(self.tps_file, mode='r')
        # get initial conditions, and stick them in files for each task
        initial_snapshots = self.select_snapshots(tps_storage)
        engine = self.most_probable_engine(initial_snapshots)
        randomizer = create_randomizer(engine)
        if self.output:
            out_storage = paths.Storage(self.output, mode='w')
            out_storage.save(engine)
            out_storage.tags['randomizer'] = randomizer
            out_storage.save(paths.Trajectory(initial_snapshots))
        task_files = self.write_task_files(initial_snapshots, engine,
                                           randomizer)
        # tps_storage.close()
        # this is the part that will be executed remotely
        result_list = [
            file_snapshots_committor(task_file, self.n_per_snapshot)
            for task_file in task_files
        ]

        runs = self.create_individual_runs(tps_storage, result_list)
        tps_storage.close()
        if self.output:
            out_storage.tags['individual_runs'] = runs
            out_storage.close()
        return initial_snapshots, runs

    def select_snapshots(self, tps_storage):
        snapshots = []
        snap_uuids = set([])
        print("Preparing initial snapshots")
        pbar = tqdm(total=self.n_snapshots)
        while len(snapshots) < self.n_snapshots:
            # OPS stores snapshots in pairs (one with reversed velocities),
            # which is what the divide/multiply by 2 handles; snapshot 2n and
            # 2n*1 have the same coordinates
            n_stored_snapshots = len(tps_storage.snapshots)
            random_choice = random.randint(0, (n_stored_snapshots // 2) - 1)
            snap = tps_storage.snapshots[random_choice * 2]
            in_state = any(state(snap) for state in self.states)
            if not in_state and snap.__uuid__ not in snap_uuids:
                snapshots.append(snap)
                snap_uuids.add(snap.__uuid__)
                pbar.update(1)
        pbar.close()
        return snapshots

    @staticmethod
    def most_probable_engine(snapshots):
        engines = {snap.engine.__uuid__: snap.engine for snap in snapshots}
        engine_counter = collections.Counter(snap.engine.__uuid__
                                             for snap in snapshots)
        if len(engine_counter) > 1:
            warnings.warn("More than one engine found.")
        return engines[engine_counter.most_common()[0][0]]

    def write_task_files(self, snapshots, engine, randomizer):
        print("Writing task files")
        filename_template = "committor_task_{:06d}.nc"
        pbar = tqdm(total=self.n_snapshots)
        task_files = []
        file_num = 0
        for (snap_num, snapshot) in enumerate(snapshots):
            filename = filename_template.format(file_num)
            storage = paths.Storage(filename, mode='w')
            storage.tag['initial_snapshots'] = [snapshot]
            storage.tag['states'] = self.states
            storage.tag['engine'] = engine
            storage.tag['randomizer'] = randomizer
            storage.close()
            task_files.append(filename)
            file_num += 1
            pbar.update(1)
        pbar.close()
        return task_files

    def create_individual_runs(self, tps_storage, task_results):
        name_to_volume = {state.name: state for state in self.states}
        def load_uuid(store, uuid):
            return store[store.index[uuid]]

        # the next part is where we turn futures into real results
        run_results = []
        for result in task_results:
            local_result = [
                (load_uuid(tps_storage.snapshots, uuid), name_to_volume[name])
                for (uuid, name) in result
            ]
            run_results += local_result

        # this will use paths.ShootingPointAnalysis.from_individual_runs
        return run_results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('tps_file', type=str,
                        help="file from which initial snapshots come")
    parser.add_argument('--n-snapshots', type=int,
                        help="number of snapshots to draw from file")
    parser.add_argument('--n-per-snapshot', type=int,
                        help="number of shots per snapshots")
    parser.add_argument('--state', '-S', type=str, action='append',
                        help="name of state in tps_file; specify multiple")
    parser.add_argument('--output', '-o', type=str, default=None,
                        help="output file")
    opts = parser.parse_args()

    kwargs = dict(vars(opts))
    # del kwargs['output']
    kwargs['state_names'] = kwargs.pop('state')
    committor_sim = ParallelCommittorSimulation(**kwargs)

    initial_snapshots, runs = committor_sim.sequential_main()

    if opts.output is None:
        spa = paths.ShootingPointAnalysis.from_individual_runs(runs)
        print(spa.to_pandas())

