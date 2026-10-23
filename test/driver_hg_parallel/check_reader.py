#!/usr/bin/env python3
"""Small, inspectable reader check; saves inputs, assignments and logs.

Run from any directory. --build names an already built MPI Autotools Zoltan
build. Uses only Python's standard library, an MPI C compiler, and mpirun.
"""
import argparse
import json
from pathlib import Path
import re
import shlex
import shutil
import subprocess

import run as checks


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def run_logged(command, log_path, cwd=None):
    try:
        output = checks.run(command, cwd=cwd)
    except Exception as error:
        log_path.write_text(shlex.join(command)+'\n'+str(error)+'\n')
        raise
    log_path.write_text(shlex.join(command)+'\n'+output)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--driver', type=Path, help='Default: BUILD/src/driver/zdrive.exe')
    parser.add_argument('--output', type=Path, required=True, help='New directory; all results are retained')
    parser.add_argument('--ranks', default='1,2', help='Comma-separated MPI rank counts')
    parser.add_argument('--mpicc', default='mpicc', help='Use the compiler matching the driver MPI')
    parser.add_argument('--mpiexec', default='mpirun')
    parser.add_argument('--mpi-args', default='', help='Extra launcher options, e.g. --mpi-args="--oversubscribe"')
    args = parser.parse_args()
    build = args.build.resolve()
    driver = (args.driver or build/'src/driver/zdrive.exe').resolve()
    require(driver.is_file(), f'Driver does not exist: {driver}; build zdrive first')
    ranks = [int(n) for n in args.ranks.split(',')]
    require(all(n > 0 for n in ranks), 'Rank counts must be positive')
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    here = Path(__file__).resolve().parent
    repo = here.parents[1]
    probe = output/'reader_compare'
    command = [*shlex.split(args.mpicc), '-O2', '-g', '-ffunction-sections', '-fdata-sections',
               f'-I{build}/src/include', f'-I{repo}/src/include',
               f'-I{repo}/src/driver', f'-I{repo}/src/ch',
               str(here/'reader_compare.c'), str(repo/'src/driver/dr_hg_readfile.c'),
               str(repo/'src/driver/dr_err.c'), str(repo/'src/driver/dr_util.c'),
               str(repo/'src/ch/ch_dist_graph.c'), str(repo/'src/ch/ch_init_dist.c'),
               '-Wl,--gc-sections', f'-L{build}/src', f'-Wl,-rpath,{build}/src',
               '-lzoltan', '-lm', '-o', str(probe)]
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (output/'build.log').write_text(shlex.join(command)+'\n'+result.stdout)
    require(result.returncode == 0, f'Probe compilation failed; see {output}/build.log')

    expected = json.loads((here/'fixtures/expected.json').read_text())
    report = {'driver': str(driver), 'input_checks': [], 'partition_checks': []}
    for name, examples in expected.items():
        source = here/'fixtures'/(name+'.hg')
        # A hand-checkable oracle validates the independent evaluator itself.
        for example in examples.values():
            parts = {int(gid): part for gid, part in example['assignments'].items()}
            actual = checks.evaluate_file(source, parts, len(example['part_weights']))
            require(actual == {k: v for k, v in example.items() if k != 'assignments'},
                    f'Objective evaluator disagrees with hand-worked example: {name}')
        print(f'PASS {name}: hand-worked cut, km1 and balance examples', flush=True)
        for nr in ranks:
            mpi = [*shlex.split(args.mpiexec), '-np', str(nr), *shlex.split(args.mpi_args)]
            for distribution, code in [('linear', 1), ('cyclic', 2), ('owner', 3)]:
                run_logged([*mpi, str(probe), str(source), str(code), '-1'],
                           output/f'{name}-{nr}-{distribution}-arrays.log')
                report['input_checks'].append({'fixture': name, 'ranks': nr,
                                                'distribution': distribution, 'equal': True})
            for objective in ('connectivity', 'hyperedges'):
                results = []
                for mode in ('old', 'parallel'):
                    directory = output/f'{name}-{nr}-{objective}-{mode}'
                    directory.mkdir()
                    shutil.copyfile(source, directory/'graph.hg')
                    config = checks.config(mode, objective).replace('num_global_parts=4', 'num_global_parts=2')
                    if name == 'unweighted-one':
                        config = config.replace('obj_weight_dim=1', 'obj_weight_dim=0')
                        config = config.replace('edge_weight_dim=1', 'edge_weight_dim=0')
                    (directory/'zdrive.inp').write_text(config)
                    log = run_logged([*mpi, str(driver), 'zdrive.inp'], directory/'run.log', cwd=directory)
                    require('DRIVER:  Zoltan_LB_Partition time =' in log, f'Partitioning failed: {directory}')
                    if mode == 'parallel':
                        require('Read hypergraph in parallel using MPI-IO' in log,
                                'The rebuilt driver did not activate the parallel reader')
                    parts = checks.assignments(directory, nr)
                    metrics = checks.evaluate_file(source, parts, 2)
                    for label, key in [('CUTN', 'cut_net'), ('CUTL', 'connectivity')]:
                        values = re.findall(r'\b'+label+r' \(.*?\):\s+([\d.eE+-]+)', log)
                        require(values and abs(float(values[-1])-metrics[key]) <= 0.00051,
                                f'Independent {label} disagrees with driver evaluation: {directory}')
                    require(metrics['imbalance'] <= 1.5 + 1e-6, f'Balance constraint violated: {directory}')
                    (directory/'assignments.json').write_text(json.dumps(parts, indent=2)+'\n')
                    (directory/'metrics.json').write_text(json.dumps(metrics, indent=2)+'\n')
                    results.append((parts, metrics))
                require(results[0] == results[1], f'Reader results differ: {name}, {nr} ranks, {objective}')
                item = {'fixture': name, 'ranks': nr, 'objective': objective,
                        'different_assignments': 0, **results[0][1]}
                report['partition_checks'].append(item)
                print(f'PASS {name}, {nr} ranks, {objective}: identical assignments; '
                      f'CUTN={item["cut_net"]:g}, km1={item["connectivity"]:g}, '
                      f'imbalance={item["imbalance"]:g}', flush=True)
    (output/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    print(f'PASS: {len(report["input_checks"])} exact input comparisons and '
          f'{len(report["partition_checks"])} final partition comparisons. Results: {output}')


if __name__ == '__main__':
    main()
