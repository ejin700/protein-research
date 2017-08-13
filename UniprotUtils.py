import pandas as pd

def get_sites(seqrec, stype='glycosylation'):
    """
    pull out sequence meta data
    """
    features = [f for f in seqrec.features if stype in f.type]
    return features

def features_to_frame(seq,features):
    """
    features list to features frame
    Note: evidence is actually key number.  need to use evidence series
    """
    L = []
    for f in features:
        d = f.qualifiers.copy()
        d['start'] = f.location.start.position
        d['end'] = f.location.end.position
        d['val'] = f.extract(seq)
        L.append(d)
    df = pd.DataFrame(L).rename(columns={'evidence':'evidence_key'})
    return df

def get_evidence_map(evidence):
    """
    create mapping between evidence key and evidence type
    """
    df = pd.DataFrame(evidence).rename(columns={'type': 'evidence_type'})
    df['key'] = df['key'].astype(int)
    return df.set_index('key')['evidence_type']


def split_string_mult_rows(s, name=None, split_on=' '):
    """
    take a column of df, split, and add rows
    """
    return (s
            .str
            .split(expand=True)
            .stack()
            .reset_index(drop=True,level=1)
            .rename(name)
            )


def evidence_keys_to_evidence_type(ekeys, emap, split_on=' ', recombine=True):
    """
    create series with same indexed like ekeys, with etype mapping
    """
    # mapping
    s = emap[ekeys]
    
    # reindex
    s.index = ekeys.index
    
    if recombine:
        s = s.groupby(level=0).apply(lambda x: split_on.join(x.tolist())).rename(emap.name)
        
    return s


def get_feature_frame(seqrec, stype='glycosylation', filter_val='N'):
    """
    create a dataframe with features
    """
    
    # features
    features = get_sites(seqrec,stype)
    df = features_to_frame(str(seqrec.seq), features)

    if filter_val is not None:
        df = df.query('val=="{}"'.format(filter_val))
    # evidence
    emap = get_evidence_map(seqrec.annotations['evidence'])
    
    # evidence_keys
    ekeys = split_string_mult_rows(df['evidence_key'], 'key').astype(int)
    #return df, emap, ekeys
    
    df = df.join(evidence_keys_to_evidence_type(ekeys, emap))
    
    return df
