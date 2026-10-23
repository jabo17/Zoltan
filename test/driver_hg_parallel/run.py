#!/usr/bin/env python3
"""Compare reader arrays and independently evaluate final PHG assignments.

Requires an MPI zdrive and the reader_compare probe (see README.md).
All generated inputs and outputs live in a temporary directory.
"""
import argparse
import math
import json
import statistics
from pathlib import Path
import re
import subprocess
import tempfile


def run(command, cwd=None, timeout=60, expect_error=False):
    result = subprocess.run(command, cwd=cwd, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=timeout)
    if result.returncode and not expect_error:
        raise RuntimeError(f"{' '.join(map(str, command))}\n{result.stdout}")
    return result.stdout


def fixture(path, nv=17, base=1, vd=0, ed=0, patoh=False, empty=False, long=False):
    edges = [[(e*3+j*5) % nv + base for j in range(3+e % 4)] for e in range(19)]
    if empty:
        edges[2] = []
        edges[-1] = []
    if long:
        edges[5] = [j % nv + base for j in range(80000)]
    ew = [[float(1+(e+j) % 5)/4 for j in range(ed)] for e in range(len(edges))]
    vw = [[float(1+(v+j) % 3)/2 for j in range(vd)] for v in range(nv)]
    pins = sum(map(len, edges))
    if patoh:
        code = int(bool(vd)) + 2*int(bool(ed))
        header = f'{base} {nv} {len(edges)} {pins} {code}' + (f' {vd}' if vd else '')
    else:
        header = f'{nv} {len(edges)} {pins} {vd}{ed}'
    data = '% equality fixture\n' + header + '\n'
    for edge, weights in zip(edges, ew):
        data += ' '.join([*(str(w) for w in weights), *map(str, edge)]) + '\n'
    if vd:
        data += ''.join(' '.join(map(str, row))+'\n' for row in vw)
    path.write_text(data)
    return edges, ew, vw


def assignments(directory, nranks):
    parts = {}
    for rank in range(nranks):
        path = directory / f'graph.out.{nranks}.{rank}'
        # Driver naming zero-pads ranks only at >=10 processes.
        text = path.read_text().split('GID\tPart\tPerm\tIPerm\n', 1)[1]
        for line in text.splitlines():
            fields = list(map(int, line.split()))
            if len(fields) != 4:
                continue
            gid, part, _, _ = fields
            assert gid not in parts
            parts[gid] = part
    return parts


def evaluate(parts, edges, ew, vw, base, k):
    cut = connectivity = 0.0
    for e, edge in enumerate(edges):
        occupied = {parts[v] for v in edge}
        weight = ew[e][0] if ew[e] else 1.0
        cut += weight * (len(occupied) > 1)
        connectivity += weight * max(0, len(occupied)-1)
    weights = [0.0]*k
    for gid, part in parts.items():
        weights[part] += vw[gid-base][0] if vw[gid-base] else 1.0
    return cut, connectivity, weights


def config(mode, objective='connectivity', distribution='linear', subset=-1):
    return f'''Decomposition Method = hypergraph
Zoltan Parameters = lb_approach=partition
Zoltan Parameters = hypergraph_package=phg
Zoltan Parameters = num_global_parts=4
Zoltan Parameters = seed=12345
Zoltan Parameters = imbalance_tol=1.5
Zoltan Parameters = phg_edge_size_threshold=1.0
Zoltan Parameters = phg_cut_objective={objective}
Zoltan Parameters = obj_weight_dim=1
Zoltan Parameters = edge_weight_dim=1
File Type = hypergraph, initial distribution={distribution}, initial procs={subset}{', initial read=parallel' if mode == 'parallel' else ''}
File Name = graph
Parallel Disk Info = number=0
Text Output = 1
'''


def evaluate_file(path, parts, k):
    """Compute both objectives from the original .hg, without reader output."""
    with path.open() as source:
        header = next(source)
        while header.startswith('%'):
            header = next(source)
        numbers = list(map(int, header.split()))
        if numbers[0] <= 1:
            base, nv, ne, _, *rest = numbers
            code = rest[0] if rest else 0
            vd = (rest[1] if len(rest) > 1 else 1) if code in (1, 3) else 0
            ed = int(code in (2, 3))
        else:
            nv, ne, _, code = numbers
            base, vd, ed = 1, code//10, code % 10
        assert len(parts) == nv
        cut = connectivity = 0.0
        weights = [0.0]*k
        for _ in range(ne):
            tokens = next(source).split()
            weight = float(tokens[0]) if ed else 1.0
            occupied = {parts[int(pin)] for pin in tokens[ed:]}
            cut += weight*(len(occupied) > 1)
            connectivity += weight*max(0, len(occupied)-1)
        for v in range(nv):
            weight = float(next(source).split()[0]) if vd else 1.0
            weights[parts[v+base]] += weight
    return {'cut_net': cut, 'connectivity': connectivity, 'part_weights': weights,
            'imbalance': max(weights)/(sum(weights)/k)}


def benchmark(args, ranks):
    source = args.benchmark.resolve()
    report = {'input': str(source), 'bytes': source.stat().st_size,
              'repetitions': args.repetitions, 'results': []}
    with tempfile.TemporaryDirectory(prefix='zoltan-hg-bench-') as temp:
        root = Path(temp)
        (root/'graph.hg').symlink_to(source)
        for nr in ranks:
            mpi = [args.mpiexec, '-np', str(nr), '--bind-to', 'core']
            # Check the actual benchmark's complete rank-local arrays too.
            run([*mpi, str(args.probe), str(source), '1', '-1'], timeout=180)
            times = {'old': [], 'parallel': []}
            for repetition in range(args.repetitions+1):
                # Alternate order to reduce systematic cache/order bias.
                modes = ('old', 'parallel') if repetition % 2 == 0 else ('parallel', 'old')
                for mode in modes:
                    inp = config(mode) + '\nZdrive action = 0\nZdrive debug level = 0\nText output = 0\n'
                    (root/'zdrive.inp').write_text(inp)
                    log = run([*mpi, str(args.driver), 'zdrive.inp'], cwd=root, timeout=180)
                    matches = re.findall(r'DRIVER:  Input loading time = ([\d.eE+-]+)', log)
                    assert matches, log
                    if repetition:
                        times[mode].append(float(matches[-1]))
            medians = {m: statistics.median(t) for m, t in times.items()}
            item = {'ranks': nr, 'seconds': times, 'median_seconds': medians,
                    'speedup': medians['old']/medians['parallel']}
            if args.check_partitions:
                solutions = []
                metrics = []
                for mode in ('old', 'parallel'):
                    inp = config(mode).replace('num_global_parts=4', 'num_global_parts=8')
                    inp += '\nZdrive debug level = 0\n'
                    (root/'zdrive.inp').write_text(inp)
                    log = run([*mpi, str(args.driver), 'zdrive.inp'], cwd=root, timeout=300)
                    assert 'DRIVER:  Zoltan_LB_Partition time =' in log, log
                    parts = assignments(root, nr)
                    metrics.append(evaluate_file(source, parts, 8))
                    solutions.append(parts)
                differences = sum(solutions[0][gid] != solutions[1][gid] for gid in solutions[0])
                assert differences == 0 and metrics[0] == metrics[1]
                assert metrics[0]['imbalance'] <= 1.5 + 1e-6
                item['partition_check'] = {'different_assignments': differences, **metrics[0]}
            report['results'].append(item)
            print(json.dumps(item), flush=True)
    if args.report:
        args.report.write_text(json.dumps(report, indent=2)+'\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--driver', required=True, type=Path)
    ap.add_argument('--probe', required=True, type=Path)
    ap.add_argument('--mpiexec', default='mpirun')
    ap.add_argument('--ranks', default='1,2,4')
    ap.add_argument('--benchmark', type=Path)
    ap.add_argument('--repetitions', type=int, default=5)
    ap.add_argument('--report', type=Path)
    ap.add_argument('--check-partitions', action='store_true')
    args = ap.parse_args()
    args.driver = args.driver.resolve()
    args.probe = args.probe.resolve()
    ranks = list(map(int, args.ranks.split(',')))
    if args.benchmark:
        benchmark(args, ranks)
        return
    probes = partitions = invalid = 0
    with tempfile.TemporaryDirectory(prefix='zoltan-hg-') as temp:
        root = Path(temp)
        inputs = []
        for name, kw in [
            ('ibm', {}), ('ibm-weighted', dict(vd=2, ed=2)),
            ('patoh-zero', dict(base=0, patoh=True)),
            ('patoh-one', dict(base=1, patoh=True, vd=2, ed=1)),
            ('patoh-edge', dict(base=0, patoh=True, ed=1)),
            ('patoh-vertex', dict(base=1, patoh=True, vd=1)),
            ('empty', dict(empty=True)), ('long', dict(long=True)),
        ]:
            path = root / (name+'.hg')
            fixture(path, **kw)
            inputs.append((name, path))
        # More ranks than edges; CRLF and a leading comment.
        path = root / 'tiny.hg'
        path.write_bytes(b'% tiny\r\n2 1 1 00\r\n1\r\n')
        inputs.append(('tiny', path))
        path = root / 'trailing-blanks.hg'
        fixture(path, patoh=True, base=0, vd=1, ed=1)
        with path.open('a') as output:
            output.write('\n \t\n')
        inputs.append(('trailing-blanks', path))
        path = root / 'no-edges.hg'
        path.write_text('2 0 0 00\n')
        inputs.append(('no-edges', path))
        for nr in ranks:
            mpi = [args.mpiexec, '-np', str(nr)]
            for name, path in inputs:
                for distribution in (1, 2):  # linear, cyclic
                    for subset in sorted({-1, max(1, nr-1)}):
                        run([*mpi, str(args.probe), str(path), str(distribution), str(subset)])
                        probes += 1
                if name not in ('empty', 'no-edges'):
                    run([*mpi, str(args.probe), str(path), '3', '-1'])
                    probes += 1
                if name not in ('tiny', 'no-edges'):
                    run([*mpi, str(args.probe), str(path), '0', '-1'])
                    probes += 1
            # Final assignments: run both configured cut objectives and initial
            # distributions, then evaluate the original complete input in Python.
            for base in (0, 1):
                for objective in ('connectivity', 'hyperedges'):
                    for distribution in ('linear', 'cyclic', 'owner'):
                        results = []
                        for mode in ('old', 'parallel'):
                            directory = root / f'p-{nr}-{base}-{objective}-{distribution}-{mode}'
                            directory.mkdir()
                            edges, ew, vw = fixture(directory/'graph.hg', base=base,
                                                     patoh=True, vd=1, ed=1)
                            (directory/'zdrive.inp').write_text(config(mode, objective, distribution))
                            log = run([*mpi, str(args.driver), 'zdrive.inp'], cwd=directory)
                            assert 'DRIVER:  Input loading time =' in log, log
                            parts = assignments(directory, nr)
                            assert len(parts) == 17
                            metrics = evaluate(parts, edges, ew, vw, base, 4)
                            # Ensure the independent cut computations also agree
                            # with Zoltan's reported final evaluation.
                            for label, value in zip(('CUTN', 'CUTL'), metrics[:2]):
                                reported = re.findall(r'\b'+label+r' \(.*?\):\s+([\d.eE+-]+)', log)
                                assert reported and math.isclose(float(reported[-1]), value,
                                                                  rel_tol=0, abs_tol=0.00051), (label, value, reported)
                            assert max(metrics[2]) / (sum(metrics[2])/4) <= 1.5 + 1e-6
                            results.append((parts, metrics))
                        assert results[0] == results[1], results
                        partitions += 1
            # Errors must terminate collectively, including failures on ranks
            # whose byte range falls near the end of the file.
            for body in ('2 1 1 00\n9\n', '2 1 2 00\n1\n',
                         '2 1 1 00\n1', '2 1 1 00\n1 nope\n',
                         '2 1 1 11\n1 1\n0.5\nbad\n',
                         '2 1 1 00\n', '2 -1 0 00\n'):
                directory = root / f'bad-{nr}-{invalid}'
                directory.mkdir()
                (directory/'graph.hg').write_text(body)
                (directory/'zdrive.inp').write_text(config('parallel'))
                log = run([*mpi, str(args.driver), 'zdrive.inp'], cwd=directory, expect_error=True)
                assert 'Parallel HG reader:' in log, log
                assert 'DRIVER:  Input loading time =' not in log
                invalid += 1
    print(f'PASS: {probes} exact per-rank input comparisons; {partitions} final '
          f'assignment/objective/balance comparisons; {invalid} collective rejection checks')


if __name__ == '__main__':
    main()
