# What can be found in this folder?

All my previous experiments were run with Colab Pro +, this folder contains some of the notebooks which were used during project development phase. All the notebooks in this folder are reproducible, each of them is standalone, although they use/store data at Google Drive (so if you want to rerun them, change the locations according to your personal setup). 

## Data preparation

Currently I use MIBiG dataset 3.1 in this project, the dataset itself is downloaded in notebooks.

I additionally processed the protein sequences: annotated them with HMMer with Pfam-A database and extracted domain sequences. 

The resulting collection of domains after redundant annotations removal is stored in the file [domain_data.csv](https://drive.google.com/file/d/1XFebRlF-kyw7JZyvD4VOAiPtvhj27ksO/view).

## Training notebooks

Currently there are two notebooks available:

- `train_loco_split_single_domains.ipynb`

- `train_loco_split_pair_domains.ipynb`

In the first case Siamese network gets one domain per BGC as an input, in the second case it takes a pair of consecutive domains, each of them is processed with transformer head, after that the extracted embeddings for each BGC are concatenated and cosine similarity loss between the pair of BGCs is computed.

## Post-processing

(TBA)