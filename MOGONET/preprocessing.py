import argparse
import pandas as pd
import numpy as np


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--label")
    parser.add_argument("--t")
    parser.add_argument("--m")
    parser.add_argument("--p")
    parser.add_argument("--clin")
    parser.add_argument("--train_samples")
    parser.add_argument("--test_samples")
    args = parser.parse_args()

    label=args.label
    transcriptome = pd.read_csv(args.t,sep='\t',index_col=0)
    transcriptome.columns = transcriptome.columns.astype(int)
    
    methylome= pd.read_csv(args.m,sep='\t',index_col=0)
    methylome.columns = methylome.columns.astype(int)
    
    proteome = pd.read_csv(args.p,sep='\t',index_col=0)
    proteome.columns = proteome.columns.astype(int)
    
    
    clinical_raw = pd.read_csv(args.clin,sep='\t',index_col=0)
    dict_clinical = clinical_raw.reset_index().groupby(label)['SUBJNO'].apply(list).to_dict()
    all_samples_clinical = {v for i in dict_clinical.values() for v in i}
    samples_common_omics = set.intersection(set(transcriptome.columns), set(methylome.columns), set(proteome.columns))
  

    train_samples = set([int(l.strip()) for l in  open(args.train_samples).readlines()]).intersection(samples_common_omics)
    test_samples = set([int(l.strip()) for l in  open(args.test_samples).readlines()]).intersection(samples_common_omics)
    dict_clinical_r = clinical_raw.loc[:,label].to_dict()

    all_samples_omics_clinical = samples_common_omics.intersection(all_samples_clinical)

    transcriptome_filt = transcriptome.loc[:,list(all_samples_omics_clinical)].T
    methylome_filt = methylome.loc[:,list(all_samples_omics_clinical)].T
    proteome_filt = proteome.loc[:,list(all_samples_omics_clinical)].T

    label_train = list(map(lambda x:1 if dict_clinical_r[x]==1 else 0,train_samples))
    np.savetxt("labels_tr.csv",np.array(label_train),delimiter=',')
    label_test = list(map(lambda x:1 if dict_clinical_r[x]==1 else 0,test_samples))
    np.savetxt("labels_te.csv",np.array(label_test),delimiter=',')
    
    i=1
    dfs = [transcriptome_filt, methylome_filt, proteome_filt]
    for df in dfs:
        np.savetxt("{}_tr.csv".format(i),df.loc[list(train_samples),:].to_numpy(),delimiter=',')
        np.savetxt("{}_te.csv".format(i),df.loc[list(test_samples),:].to_numpy(),delimiter=',')
        i += 1
