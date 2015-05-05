#!/usr/bin/env python3


def translate(pairlist, output, mapoutput):
    """Translate a pairlist into mfinder numbers and a map file."""
    mapping = {}
    start = 1

    with open(output, 'w') as f:
        for src, dst in pairlist:
            if src not in mapping:
                mapping[src] = start
                start += 1

            if dst not in mapping:
                mapping[dst] = start
                start += 1

            if src != dst:  # disallow self loops for mfinder
                print('%d %d 1' % (mapping[src], mapping[dst]), file=f)

    with open(mapoutput, 'w') as f:
        for gene, number in mapping.items():
            print('%s %d' % (gene, number), file=f)
