import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import pandas as pd
from msda import batch_normalization as bn
from msda import preprocessing as pr

# Preprocess data from mass spec core first
# -----------------------------------------
#dfs = pr.preprocess_dataset(data_filename.xlsx)
#dfs.to_csv(data_filename.csv)

def read_norm(data_filename, metadata_filename, batch_name):
    """Read in preproceesed data for a given batch
       and normalize samples with respect to corresponding bridge
    
    Parameters
    ----------
    data_filename : str
        Path to  preprocessed data file (.csv)
    metadata_filename : str
        Path to metadata file
    batch_name : str
       Batch name

    Returns
    -------
    dfs : pandas.DataFrame
       Bridge normalized mass spec data
    """

    # NOTE: Ensure that the file has already being preprocessed once using
    # the commented out function belo
    #dfs = pr.preprocess_dataset(data_filename)

    dfs = pd.read_csv(data_filename)
    dfm = pd.read_csv(metadata_filename)
    dfm = dfm[dfm.batch == batch_name].copy()
    mdict =  {m:c for m, c in zip(dfm.tmt_label, dfm.cell_line)}
    dfs = dfs.rename(columns=mdict)
    samples = dfm['sample'].tolist()
    
    # Normalize summed intensities of each protein
    # by Number of peptides reported
    peptide_col = 'Set%s Peptides' % batch_name[-1]
    if 'Number_of_peptides' in dfs.columns.tolist():
        dfs = dfs.rename(columns={'Number_of_peptides' : peptide_col})
    dfs[samples] = dfs[samples].div(dfs[peptide_col], axis=0)

    # Get sample name of bridge sample. If 2 exist, get old bridge
    bridge_samples = [s for s in samples if 'Bridge' in s]
    if len(bridge_samples) > 1:
        bridge_sample = [s for s in bridge_samples if 'Old' in s][0]
    else:
        bridge_sample = bridge_samples[0]

    # Normalize samples relative to Bridge Sample
    dfn = bn.normalize_within_batch(dfs, samples, control=bridge_sample)
    return dfn



def normalize_between_batches(df, df_ref, samples, control='Bridge'):
    """Function normalizes samples across batches by ratio of Bridge samples
    Parameters
    ----------
    df : pandas dataframe
        dataframe corresponding to mass spec results of one batch
        that is already been normalzied wrt its bridge
    df_ref : pandas dataframe
        dataframe whose Bridge serves as primary normalizer
    samples : list[str]
        list of sample names
    Returns
    -------
    df_norm_inter : pandas dataframe
       dataframe corresponding to mass spec results of one batch
       normalized by its ratio of bridge samples across batches
    """

    true_samples = samples#[s for s in samples if not "Bridge" in s]
    #true_samples.remove(control)
    refset_bridge_mean = float(df_ref.loc[:, 'Bridge4'].sum())
    set_bridge = [s for s in samples if s.startswith("Bridge")
                  or "New Bridge" in s][0]
    set_bridge_mean = float(df.loc[:, set_bridge].sum())
    nrm = refset_bridge_mean / set_bridge_mean
    df_norm_inter = df.copy()
    df_norm_inter.loc[:, true_samples] = nrm * df_norm_inter.loc[:,
                                                                 true_samples]
    if 'Uniprot_Id' not in df_norm_inter.columns.tolist():
        df_norm_inter['Uniprot_Id'] = [s.split('|')[1] for s in
                                       df_norm_inter['Protein Id'].tolist()]
    df_norm_inter.index = df_norm_inter.Uniprot_Id
    # Take mean of duplicate protein entries
    df_norm_inter = df_norm_inter.groupby(df_norm_inter.index).mean()
    return df_norm_inter




