import urllib.request
import shutil
import pandas as pd
from Bio import GenBank
from collections import defaultdict
import sys
from pathlib import Path



def download_mibig(version: str = "4.0", base_dir: Path = Path('.')):
    DATADIR = base_dir / f"mibig_{version}"
    DATADIR.mkdir(exist_ok=True)

    if not (DATADIR / "mibig_json.tar.gz").exists():
        urllib.request.urlretrieve(
            f"https://dl.secondarymetabolites.org/mibig/mibig_json_{version}.tar.gz", 
            (DATADIR / "mibig_json.tar.gz").as_posix()
        )
        urllib.request.urlretrieve(
            f"https://dl.secondarymetabolites.org/mibig/mibig_prot_seqs_{version}.fasta", 
            (DATADIR / "mibig_prot_seqs.fasta").as_posix()
        )
    if not (DATADIR / "mibig_gbk.tar.gz").exists():
        urllib.request.urlretrieve(
            f"https://dl.secondarymetabolites.org/mibig/mibig_gbk_{version}.tar.gz", 
            (DATADIR / "mibig_gbk.tar.gz").as_posix()
        )
    
    if not (DATADIR / "mibig_json_{version}").exists():
        shutil.unpack_archive(
            (DATADIR / "mibig_json.tar.gz").as_posix(),
            DATADIR.as_posix(),
            "tar"
        )
    if not (DATADIR / "mibig_gbk_{version}").exists():
        shutil.unpack_archive(
            (DATADIR / "mibig_gbk.tar.gz").as_posix(),
            DATADIR.as_posix(),
            "tar"
        )
    return DATADIR



def read_genbank_file(filepath):
    accession = filepath.stem
    with open(filepath.as_posix()) as handle:
        for record in GenBank.parse(handle):
            assert "_".join(record.accession) == accession, f"{accession}, {record.accession}"
            # print(record.id + "_" + seguid(record.seq))
            yield record


def record2list(record):
    qualifiers_list = []
    for feat in record.features:
        if feat.key !='aSDomain':
            continue
        qualifiers = dict()
        for q in feat.qualifiers:
            key = q.key[1:-1]
            value = eval(q.value)
            if key in qualifiers:
                if isinstance(qualifiers[key], list):
                    qualifiers[key].append(value)
                else:
                    qualifiers[key] = [qualifiers[key], value]
            else:
                qualifiers[key] = value
        qualifiers_list.append(qualifiers)

    return qualifiers_list    


def extract_domain_information(mibig_dir: Path, savedir="preprocessed_features", force=False) -> Path:
    if not mibig_dir.exists():
        return
    version = mibig_dir.name.split("_")[1]
    gbk_dir = mibig_dir / f"mibig_gbk_{version}"
    gbk_files = list(gbk_dir.glob("*.gbk"))
    
    features_dir = mibig_dir / savedir
    if features_dir.exists() and not force:
        return
    features_dir.mkdir(exist_ok=True)
    for sample_file in gbk_files:
        for record in read_genbank_file(sample_file):
            accession = "_".join(record.accession)
            savepath = features_dir / f"{accession}.csv"
            if savepath.exists() and not force:
                continue
            qualifiers_list = record2list(record)
            df = pd.DataFrame(qualifiers_list)
            df.to_csv(savepath, index=None)
            # break
    return features_dir


def collect_domain_information(domain_info_dir: Path, savepath: Path):
    # selected_columns = ['aSDomain', 'locus_tag', 'translation']
    domain_data = []
    for feature_path in domain_info_dir.glob('*.csv'):
        try:
            accession = feature_path.stem
            df = pd.read_csv(feature_path)
            # domain_names[accession] = df.aSDomain.values
            # selected_columns = df.columns
            # df = df.loc[:, selected_columns]
            df['accession'] = accession
            domain_data.extend(df.to_dict('records'))
        except:
            # print(accession)
            continue
    collected_df = pd.DataFrame(domain_data)  # [['accession'] + selected_columns]
    collected_df.to_csv(savepath, index=False)


if __name__ == "__main__":
    data_dir = Path("../data")
    data_dir.mkdir(exist_ok=True)
    print("Download data...")
    mibig_dir = download_mibig(version="4.0", base_dir=data_dir)
    print('Start domain info extraction...')
    features_dir = extract_domain_information(mibig_dir)
    if features_dir is not None:
        collect_domain_information(features_dir, mibig_dir / "gbk_domain_info.csv")
    print("Everything finished")