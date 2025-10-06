import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.legend_handler import HandlerTuple

import re 
import os
import sys
import argparse as ap

sns.set_theme(style="whitegrid")

# data_folder = "../dataset"
# test_id = 'ismb2025'
# test_id = 'sept25-multiple'

data_folder = "../sept25"
types = ['binning']#, 'labeling']
type_tags = ['bin']#, 'lab']

df_sample = pd.read_csv("samples.tsv", sep="\t", header=None)
df_sample.columns = ['sample_id']
samples = df_sample['sample_id'].tolist()
res_columns = ['sample_id', 'model','bintype', 'tool', 'pty', 'thr', 'pctid', 'precision', 'recall', 'f1']
df_sample = pd.DataFrame(columns=res_columns)


tests = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', # with penalty score
         'n', 'o', 'p', 'q']
failed_tests = {'a':[], 'b':[], 'c':[], 'd':[], 'e':[], 'f':[], 'g':[], 'h':[], 'i':[], 'j':[], 'k':[],  'l':[], 'm':[],
                'n':[], 'o':[], 'p':[], 'q':[]} 

## MAP test to thr, pctid, pty ###################
tests_mapping_tsv = pd.read_csv("tests.csv", sep=",", header=0)
tests_mapping_tsv = tests_mapping_tsv.set_index('testid')
# print(tests_mapping_tsv)

test_map = {}
for test in tests:
    if test in tests_mapping_tsv.index:
        row = tests_mapping_tsv.loc[test]
        test_map[test] = {'thr': row['thr'], 'pctid': row['pctid'], 'pty': row['pty']}
    else:
        print(f"Warning: test {test} not found in tests.csv")
        test_map[test] = {'thr': "N/A", 'pctid': "N/A", 'pty': "N/A"}
# print(test_map)
##################################################

def read_score_file(samplefolder, score):
    cols = ["sample", "species", "model", "bintype", "tool", "prec", "rec", "f1"]
    out_tuple = pd.DataFrame(columns=cols)
    _sample = "null"
    for filename in os.listdir(samplefolder):
        fields = filename.split(".")
        if (not filename.endswith(".txt")): # salto i file non txt, non sono file di score
            continue
        _species = "null"
        _sample = fields[0] #nome del sample
        _type = fields[1] # tipo in input {ske, uni, pan}
        _thr = fields[2] # threshold (default 1)
        _model = fields[3] # modello {pbf, ml}, default pbf
        assert(_model == "pbf")
        _score = fields[-2] # score {lab, bin} labeling o binning
        _ref = fields[-3] # reference {uni, ske}
        _bins = fields[-4] # tipo di binning {pred, nve, ovl} pred= predizione del modello, nve= naive, ovl= graph-overlap

        if _type == 'pan':
            if _score == "lab" and score == "lab": # se il tipo di score del file è del tipo richiesto
                            assert(filename.endswith(".txt"))
                            with open(samplefolder+ "/" + filename) as file:
                                lines = file.readlines()
                                prec, rec, f1 = float(lines[-4].split("\t")[-1].strip()), float(lines[-3].split("\t")[-1].strip()), float(lines[-2].split("\t")[-1].strip())

            elif _score == "bin" and score == "bin":
                assert(filename.endswith(".txt"))
                with open(samplefolder+ "/" + filename) as file:
                    lines = file.readlines()
                    # print(lines)
                    prec, rec, f1 = float(lines[-3].split("\t")[-1].strip()), float(lines[-2].split("\t")[-1].strip()), float(lines[-1].split("\t")[-1].strip())
            else:
                continue
            out_tuple = pd.concat([pd.DataFrame([[_sample, _species, _model, _bins, f"{_thr}.{_type}", prec, rec, f1]], columns=cols), out_tuple if not out_tuple.empty else None], ignore_index=True)
            # print(out_tuple)
        else:
            continue
    if len(out_tuple) > 1:
        out_tuple = out_tuple.groupby(["sample", "model", "species", "bintype", "tool"]).mean(numeric_only=True).reset_index()
        # out_tuple.insert(loc=3, column="tool", value = f"{_thr}.{_type}")

    return out_tuple


def load_results(data_folder, testid, df_sample, samples):
    ### i'm in a test folder like "pange-data/a/", iterate over all samples in df_sample
    for sample_id in samples:
        # result_file = os.path.join(data_folder, testid, f"{sample_id}_results.tsv")
        # if os.path.exists(result_file):
        #     df_result = pd.read_csv(result_file, sep="\t")
        #     # Process the result dataframe as needed
        # else:
        #     print(f"Warning: result file {result_file} not found")
        df_result = read_score_file(os.path.join(data_folder, testid, sample_id), "bin")
        # print(df_result)
        # print(f"Processing sample {sample_id} in test {testid}")
    
        if df_result.empty:
            # print(f"Warning: No results for sample {sample_id} in test {testid}")
            failed_tests[testid].append(sample_id)
            continue
        # print(df_result)
        # print(df_result)
        # print(df_result)
        for t, tag in zip(types, type_tags):
            # print(f"Processing type {t} with tag {tag}")

            # Assuming only one row per type
            res_row = df_result.iloc[0]
            # df_sample['model'].dtype = 'string']
            # from already existing df_result, add pty, thr, pctid columns
            df_result['pty'] = test_map[testid]['pty']
            df_result['thr'] = test_map[testid]['thr']
            df_result['pctid'] = test_map[testid]['pctid']
            df_result['sample_id'] = sample_id
            df_result = df_result.rename(columns={'prec': 'precision', 'rec': 'recall', 'f1': 'f1'})
            df_result = df_result[['sample_id', 'model', 'bintype', 'tool', 'pty', 'thr', 'pctid', 'precision', 'recall', 'f1']]
            # take only df_result with bintype == 'pred'
            df_result = df_result[df_result['bintype'] == 'pred']
            # print(df_result)
            df_sample = pd.concat([df_sample, df_result], ignore_index=True)

            # df_sample.loc[df_sample['sample_id'] == sample_id, 'type'] = tag  
        
    
    return df_sample

df_sample = load_results(data_folder, tests[0], df_sample, samples)
# print(df_sample)

for test in tests[1:]:
    df_sample = load_results(data_folder, test, df_sample, samples)

print(df_sample)

print(failed_tests)
# count failed tests
total_failed = sum(len(v) for v in failed_tests.values())
print(f"Total failed tests: {total_failed}")
# count fail tests per test
for test, failed in failed_tests.items():
    print(f"Test {test}: {len(failed)} failed samples")


test_colors = {
    'a': 'tab:blue', 'b': 'tab:orange', 'c': 'tab:green', 'd': 'tab:red', 'e': 'tab:purple',
    'f': 'tab:brown', 'g': 'tab:pink', 'h': 'tab:gray', 'i': 'tab:olive', 'l': 'tab:cyan',
    'm': 'lightblue', 'n': 'lightgreen', 'o': 'lightcoral', 'p': 'lightpink', 'q': 'lightgray',
    'r': 'lightyellow', 'j': 'lightcyan', 'k': 'lavender'}
def boxplot(df, output_file, tests, metric):
    plt.figure(figsize=(12, 8))
    ax = sns.boxplot(x='test', y=metric, data=df, palette=test_colors, order=tests, hue='test', dodge=False)
    ax.set_title(f'Comparison of {metric} across tests')
    ax.set_xlabel('Test')
    ax.set_ylabel(metric.capitalize())
    # plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

df_sample['test'] = df_sample.apply(lambda row: next((t for t in tests if row['pty'] == test_map[t]['pty'] and str(row['thr']) == str(test_map[t]['thr']) and str(row['pctid']) == str(test_map[t]['pctid'])), 'unknown'), axis=1)
# print(df_sample)
# remove column model, bintype, tool, pty, thr, pctid
df_sample = df_sample.drop(columns=['model', 'bintype', 'tool', 'pty', 'thr', 'pctid'])
# print only rows with test != 'unknown'


def plot(testset, title, df):
    if len(testset) == 0:
        testset = tests
    boxplot(df_sample, f"{title}_precision.png", testset, "precision")
    boxplot(df_sample, f"{title}_recall.png", testset, "recall")
    boxplot(df_sample, f"{title}_f1.png", testset, "f1")

# tests with less than 4 failed tests are inserted in variable passed_tests


def subsample(df, tests):
    df_sample_subset = df[df['test'].isin(tests)]
    samples_passed_all = set(df_sample_subset['sample_id'])
    for test in tests:
        samples_passed_all = samples_passed_all.intersection(set(df[df['test'] == test]['sample_id']))
    df_sample_subset = df_sample_subset[df_sample_subset['sample_id'].isin(samples_passed_all)]
    return df_sample_subset

subtests = ['a', 'b', 'c', 'j', 'k', 'l']
df_sample_subset = subsample(df_sample, subtests)
print(f"Subsample with tests {subtests}, number of samples: {len(df_sample_subset['sample_id'].unique())}")
plot(subtests, "best_subset", df_sample_subset)

# everysub = tests
# df_every = subsample(df_sample, everysub)
# print(f"Subsample with all tests {everysub}, number of samples: {len(df_every['sample_id'].unique())}")
# plot(everysub, "every_subset", df_every)

passed_tests = [test for test in tests if len(failed_tests[test]) < 4]
df_passed_subset = subsample(df_sample, passed_tests)
print(f"Subsample with passed tests {passed_tests}, number of samples: {len(df_passed_subset['sample_id'].unique())}")
plot(passed_tests, "passed_sub", df_passed_subset)
# data_folder = "../dataset"

plot(tests, "all_tests", df_sample)

def print_mean_metric_std_dev(df, tests, metric):
    means = {}
    std_devs = {}
    for test in tests:
        df_test = df[df['test'] == test]
        mean_value = df_test[metric].mean()
        std_dev_value = df_test[metric].std()
        means[test] = mean_value
        std_devs[test] = std_dev_value

    #print std_dev next to mean
    for test in tests:
        if test in means and test in std_devs:
            print(f"Test {test}: Mean {metric} = {means[test]:.4f}, Std Dev = {std_devs[test]:.4f}")
        else:
            print(f"Test {test}: No data available")
def print_metrics(df, tests):
    print_mean_metric_std_dev(df, tests, "precision")
    print_mean_metric_std_dev(df, tests, "recall")
    print_mean_metric_std_dev(df, tests, "f1")

print("Metrics for best subset:")
print_metrics(df_sample_subset, subtests)
# print_metrics(df_every, everysub)
print("Metrics for passed tests subset:")
print_metrics(df_passed_subset, passed_tests)
print("Metrics for all tests:")
print_metrics(df_sample, tests)

