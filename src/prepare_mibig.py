import urllib.request
import shutil
import pandas as pd
from collections import defaultdict
from Bio import GenBank
from pathlib import Path
import click
import json


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


def read_json(json_path):
    with open(json_path) as f:
        data = json.load(f)
    return data


def collect_json_info(datadir, version="4.0"):
    all_data = []
    json_files = sorted(datadir.glob(f"mibig_json_{version}/*.json"))
    for json_file in json_files:
        bgc_data = read_json(json_file)
        all_data.append(bgc_data)
    return all_data


def extract_compound_activity_info(all_data):
    all_info = []
    for bgc_data in all_data:
        accession = bgc_data['accession']
        # print(list(bgc_data))
        compounds_info = bgc_data['compounds']
        for info in compounds_info:
            # print(list(info))
            bioactivities = info.get('bioactivities', [])
            for activity in bioactivities:
                # print(activity)
                if isinstance(activity['name'], dict):
                    # print(activity['name'])
                    activity['name'] = activity['name']['activity']
                sel_columns = "name	observed	references".split()
                activity = {col: activity.get(col) for col in sel_columns if col in activity}


                for column in ['name', 'evidence', 'classes', 'structure', 'cyclic', 'synonyms', 'moieties']:
                    activity[f"compound_{column}"] = info.get(column)
                activity['accession'] = accession
                # activity['compound_name'] = compound_name
                all_info.append(activity)
            # info['accession'] = accession
        # all_info.append(info)
        # print(compound_info)
    return all_info


def extract_taxonomy_info(all_data):
    all_info = []
    for bgc_data in all_data:
        accession = bgc_data['accession']
        tax_info = bgc_data['taxonomy']
        tax_info['accession'] = accession
        # todo: add taxonomy extention
        all_info.append(tax_info)
    return all_info


def extract_biosyn_class_info(all_data):
    all_info = []
    for bgc_data in all_data:
        accession = bgc_data['accession']
        biosyn_info = bgc_data['biosynthesis']
        biosyn_classes = biosyn_info['classes']
        for info in biosyn_classes:
            info = {k: info[k] for k in "class	subclass".split() if k in info}
            info['accession'] = accession
            all_info.append(info)
    return all_info


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


def read_fasta_iter(path):
    desc, sequences = '', []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if len(sequences) > 0:
                    yield desc, "".join(sequences)
                    sequences = []
                desc = line[1:]
            else:
                sequences.append(line)
        if len(sequences) > 0:
            yield desc, "".join(sequences)


def extract_protein_seqs(path):
    all_data = []
    for desc, seq in read_fasta_iter(path):
        data = desc.split('|')
        accession = data[0].split('.')[0]
        position = data[2]
        direction = data[3]
        name = data[5]
        all_data.append({
            'accession': accession, 
            'position': position,
            'direction': direction,
            'name': name,
            'sequence': seq
        })
    return all_data


@click.command()
@click.option('--datadir', default='data', help='Path to the directory where the data should be saved')
@click.option('--version', default='4.0', help='Version of mibig dataset to be saved')
def main(datadir, version):
    data_dir = Path(datadir)
    data_dir.mkdir(exist_ok=True)
    print(f"Download data to {repr(str(data_dir))}...")
    mibig_dir = download_mibig(version=version, base_dir=data_dir)
    all_data = collect_json_info(mibig_dir, version=version)
    df = pd.DataFrame(all_data)
    # print(list(df['compounds']))
    df.to_csv(mibig_dir / "annotation.csv", index=None)
    print('Extract taxonomy info...')
    taxonomy_data = extract_taxonomy_info(all_data)
    pd.DataFrame(taxonomy_data).to_csv(mibig_dir / 'taxonomy_info.csv', index=None)

    print('Extract biosynthetic classes info...')
    biosyn_class_info = extract_biosyn_class_info(all_data)
    pd.DataFrame(biosyn_class_info).to_csv(mibig_dir / 'biosyn_classes.csv', index=None)

    print('Start bgc info extraction...')
    activity_data = extract_compound_activity_info(all_data)
    pd.DataFrame(activity_data).to_csv(mibig_dir / "compound_activity_info.csv", index=None)
    print('Start domain info extraction...')
    features_dir = extract_domain_information(mibig_dir)
    if features_dir is not None:
        collect_domain_information(features_dir, mibig_dir / "gbk_domain_info.csv")
    print('Start protein sequences extraction...')
    all_data = extract_protein_seqs(mibig_dir / "mibig_prot_seqs.fasta")
    pd.DataFrame(all_data).to_csv(mibig_dir / "protein_sequences.csv", index=None)
    print("Everything finished")



if __name__ == "__main__":
    main()
