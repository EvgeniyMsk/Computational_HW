from sys import stdout

def f_range(a, b, d):
    while a <= b:
        yield a
        a += d

xs   = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

f_xs = [1614419,
        1656832,
        1694888,
        1728506,
        1758030,
        1783225,
        1804279,
        1821299,
        1834414,
        1843768]

n = len(f_xs)

diffs = []
diffs.append(f_xs)

# calculate differences
for i in range(1, n):
    i_diff      = []
    prev_diff   = diffs[i - 1]
    for j in range(1, len(prev_diff)):
        i_diff.append(prev_diff[j] - prev_diff[j - 1])
    diffs.append(i_diff)

lens = list()
for i in range(0, len(diffs)):
    max = 0
    for j in range(0, len(diffs[i])):
        if max < len(str(diffs[i][j])):
            max = len(str(diffs[i][j]))
    lens.append(max + 3)

# print differences
for i in range(0, len(diffs[0])):
    for j in range(0, len(diffs)):
        if i < len(diffs[j]):
            stdout.write(str(diffs[j][i]) + ' '*(lens[j] - len(str(diffs[j][i]))))
    stdout.write('\n')
