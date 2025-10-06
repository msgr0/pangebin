import pandas as pd

import itertools as it

thr = [200, 350, 500, 1000, 2000, 5000]
pctid = [95, 98, 93, 90]
cutlen = [0, 100, 500, 1000]
pty = [0, 1]

tests = list(it.product(thr, pctid, cutlen, pty))
df = pd.DataFrame(tests, columns=['thr', 'pctid', 'cutlen', 'pty'])
## add a unique test id column
df.insert(0, 'test_id', range(1, len(df) + 1))
## for test ID use alfanumeric with leading zeros
df['test_id'] = df['test_id'].apply(lambda x: f'T{x:03d}')
df.to_csv('test_params.csv', index=False)
# print(len(tests))