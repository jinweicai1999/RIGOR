"""Recheck published means, SDs, BH adjustments, and standardized exact tests."""
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
import csv
import statistics

ROOT = Path(__file__).resolve().parent


def read(name):
    with (ROOT / name).open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f))


def signed_rank_p(differences):
    groups = defaultdict(list)
    for value in differences:
        if value:
            groups[abs(value)].append(value)
    ranked = []
    rank = 1
    for value in sorted(groups):
        group = groups[value]
        rank2 = 2 * rank + len(group) - 1
        ranked.extend((rank2, x > 0) for x in group)
        rank += len(group)
    distribution = {0: 1}
    for rank2, _ in ranked:
        new = defaultdict(int)
        for total, n in distribution.items():
            new[total] += n
            new[total + rank2] += n
        distribution = new
    observed = sum(rank2 for rank2, positive in ranked if positive)
    ways = 2 ** len(ranked)
    low = sum(n for s, n in distribution.items() if s <= observed) / ways
    high = sum(n for s, n in distribution.items() if s >= observed) / ways
    return min(1., 2 * min(low, high))


scores = read('per_gene_condition_scores.csv')
summaries = read('condition_summary.csv')
tests = read('paired_tests_and_BH.csv')
assert len(scores) == 26 and len({r['benchmark_id'] for r in scores}) == 26
for row in summaries:
    values = [float(r[row['condition']]) for r in scores]
    assert len(values) == int(row['n'])
    assert abs(statistics.mean(values) - float(row['mean'])) < 1e-9
    assert abs(statistics.stdev(values) - float(row['sample_sd'])) < 1e-9
for row in tests:
    if row['analysis'] != 'standardized_ablation':
        continue
    # Recover the exact rational scores underlying the stored decimal values.
    differences = [Fraction(r['standardized_standardized_full_evidence']).limit_denominator(1000000)
                   - Fraction(r['standardized_' + row['condition_b']]).limit_denominator(1000000)
                   for r in scores]
    assert abs(signed_rank_p(differences) - float(row['raw_p'])) < 1e-12
for group in ('baseline', 'standardized_ablation'):
    subset = sorted((r for r in tests if r['analysis'] == group), key=lambda r: float(r['raw_p']))
    running = 1.
    for i in range(len(subset) - 1, -1, -1):
        running = min(running, float(subset[i]['raw_p']) * len(subset) / (i + 1))
        assert abs(running - float(subset[i]['BH_adjusted_p'])) < 1e-12
print('Verified 26-gene summaries, standardized exact tests, and both BH families.')
