#!/usr/bin/env python3
import re
from pathlib import Path

files = [
    'examples/NiO/exchanges.distance.out.reference',
    'examples/NiO/exchanges.list.out.reference',
    'examples/KCuF3/exchanges.out.reference'
]

num_re = re.compile(r'([+-]?\d+\.\d+(?:[eE][+-]?\d+)?|[+-]?\d+)')

for f in files:
    p = Path(f)
    if not p.exists():
        print('missing', f)
        continue
    lines = p.read_text().splitlines()
    out = []
    i = 0
    while i < len(lines):
        line = lines[i]
        out.append(line)
        # detect summary line containing 'meV' and 'K ('
        if 'meV' in line and 'K' in line and 'distance' in line:
            # negate first numeric before 'meV' and numeric before 'K'
            # find part before '(distance'
            try:
                before_dist, after = line.split('(distance',1)
            except ValueError:
                before_dist = line
                after = ''
            # find all numbers in before_dist
            nums = list(num_re.finditer(before_dist))
            if nums:
                # negate last two numbers that correspond to meV and K
                # usually pattern: <spaces> {mev} meV = {K} K
                # find indices of numbers: take the last two numeric matches
                if len(nums) >= 2:
                    toneg = nums[-2:]
                else:
                    toneg = nums
                new_before = before_dist
                # process in reverse to not break indices
                for m in reversed(toneg):
                    s = m.group(0)
                    try:
                        val = float(s)
                        neg = -val
                        # preserve formatting: if integer-like, keep as integer
                        if re.fullmatch(r'[+-]?\d+', s):
                            new = str(int(neg))
                        else:
                            # preserve decimals count
                            if 'e' in s or 'E' in s:
                                new = ('%.'+str(6)+'e')%neg
                            else:
                                dec = len(s.split('.')[-1])
                                fmt = '%.'+str(dec)+'f'
                                new = fmt%neg
                    except:
                        new = s
                    start, end = m.span()
                    new_before = new_before[:start] + new + new_before[end:]
                new_line = new_before + '(distance' + after
                out[-1] = new_line
        # detect matrix lines: heuristics - lines that contain many numbers and not 'Orbital' words
        elif re.search(r'\d', line) and ('Orbital exchange' not in line) and ('For atom' not in line):
            # if line contains at least 2 numeric tokens separated by spaces we treat as matrix row
            nums = list(num_re.finditer(line))
            if len(nums) >= 2:
                new_line = line
                for m in reversed(nums):
                    s = m.group(0)
                    try:
                        val = float(s)
                        neg = -val
                        if re.fullmatch(r'[+-]?\d+', s):
                            new = str(int(neg))
                        else:
                            if 'e' in s or 'E' in s:
                                new = ('%.'+str(6)+'e')%neg
                            else:
                                dec = len(s.split('.')[-1])
                                fmt = '%.'+str(dec)+'f'
                                new = fmt%neg
                    except:
                        new = s
                    start, end = m.span()
                    new_line = new_line[:start] + new + new_line[end:]
                out[-1] = new_line
        i += 1
    p.write_text('\n'.join(out)+"\n")
    print('processed', f)
