#!/usr/bin/env python3
import re
from pathlib import Path

files = [
    'examples/NiO/exchanges.distance.out.reference',
    'examples/NiO/exchanges.list.out.reference',
    'examples/KCuF3/exchanges.out.reference'
]
num_re = re.compile(r'[+-]?\d+\.\d+(?:[eE][+-]?\d+)?|[+-]?\d+')

def negate_token(token, width):
    # preserve format: decimals or integer, then right-align to field width
    if re.fullmatch(r'[+-]?\d+', token):
        val = int(token)
        neg = -val
        s = str(neg)
        return s.rjust(width)
    else:
        # float
        if 'e' in token or 'E' in token:
            val = float(token)
            neg = -val
            # format with same exponent style, use 6 decimals
            s = ('%.6e') % neg
            return s.rjust(width)
        else:
            parts = token.split('.')
            dec = len(parts[1])
            val = float(token)
            neg = -val
            fmt = ('%.' + str(dec) + 'f') % neg
            return fmt.rjust(width)

for f in files:
    p = Path(f)
    if not p.exists():
        print('missing', f)
        continue
    lines = p.read_text().splitlines()
    out = []
    in_matrix = False
    for i,line in enumerate(lines):
        if 'Orbital exchange interaction matrix' in line:
            out.append(line)
            in_matrix = True
            continue
        if in_matrix:
            # matrix lines until blank line or line without digits
            if line.strip()=='' or not re.search(r'\d', line):
                in_matrix = False
                out.append(line)
                continue
            # replace each numeric token with its negation, preserving spacing via spans
            new_line = line
            matches = list(num_re.finditer(line))
            if len(matches) >= 1:
                for m in reversed(matches):
                    token = m.group(0)
                    a,b = m.span()
                    field_width = b - a
                    newtok = negate_token(token, field_width)
                    new_line = new_line[:a] + newtok + new_line[b:]
            out.append(new_line)
            continue
        # summary lines: look for pattern with 'meV' and 'K' and 'distance'
        if 'meV' in line and 'K' in line and 'distance' in line:
            # limit to this line only: replace numeric immediately before 'meV' and before 'K'
            # find 'meV' index
            idx_mev = line.find('meV')
            # find last number before idx_mev
            before_mev = line[:idx_mev]
            mnums = list(num_re.finditer(before_mev))
            new_line = line
            if mnums:
                last = mnums[-1]
                token = last.group(0)
                a,b = last.span()
                field_width = b - a
                newtok = negate_token(token, field_width)
                new_line = new_line[:a] + newtok + new_line[b:]
            # now find 'K' (the K before '(distance') - find ' K (' pattern
            kpos = new_line.find('K')
            # find number before that K (search back)
            before_k = new_line[:kpos]
            mnums2 = list(num_re.finditer(before_k))
            if mnums2:
                last2 = mnums2[-1]
                token2 = last2.group(0)
                a2,b2 = last2.span()
                field_width2 = b2 - a2
                newtok2 = negate_token(token2, field_width2)
                new_line = new_line[:a2] + newtok2 + new_line[b2:]
            out.append(new_line)
            continue
        out.append(line)
    p.write_text('\n'.join(out)+"\n")
    print('processed', f)
