#!/usr/bin/env python
import argparse

import pandas as pd


def get_neighbor(fimo, seq, start, stop, inner, outer):
    fimo[fimo.sequence_name == seq]
    fimo_neighbor = fimo[fimo.sequence_name == seq].query('start > ' +
                                                          str(stop - outer) +
                                                          ' & start <= ' +
                                                          str(stop + outer) +
                                                          ' & stop <= ' +
                                                          str(start + inner) +
                                                          ' & stop >' +
                                                          str(start - outer))
    return fimo_neighbor


def main(arg):
    annot = pd.read_csv(arg.annot, sep='\t', header=0, index_col=0)
    get_tf = annot.loc[:, ['Transcription factor']].drop_duplicates().to_dict(
        orient='dict')['Transcription factor']
    annot.index = annot.index.str.replace('_MOUSE.*', '')
    get_annot = annot.loc[:, ['Transcription factor']].drop_duplicates(
    ).reset_index().set_index('Transcription factor').to_dict(
        orient='dict')['Model']

    fimo = pd.read_csv(arg.fimo, sep='\t', header=0, comment='#')
    fimo_target = fimo[fimo.motif_id.str.startswith(get_annot[arg.target])]
    fimo_neighbors = pd.concat([
        get_neighbor(fimo, row.sequence_name, row.start, row.stop, arg.inner,
                     arg.outer) for i, row in fimo_target.iterrows()
    ],
                               axis=0)
    # fimo_neighbors = fimo_neighbors.iloc[:, ]
    fimo_neighbors = fimo_neighbors.join(
        fimo_neighbors['sequence_name'].str.split(
            '::', 1, expand=True).rename(columns={
                0: 'target_gene',
                1: 'location'
            }))
    fimo_neighbors.drop(['sequence_name', 'motif_alt_id'], axis=1, inplace=True)
    fimo_neighbors['tf'] = [
        get_tf[f] for f in fimo_neighbors['motif_id'].values
    ]
    fimo_neighbors.reset_index(drop=True).to_csv(arg.output)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='find coregulators based on FIMO output')

    parser.add_argument(
        '-f',
        '--fimo',
        type=str,
        default='fimo.tsv',
        help='Input file. (A FIMO output TSV) (default: fimo.tsv)')

    parser.add_argument(
        '-t',
        '--target',
        type=str,
        default=None,
        help=
        'The Target TF which we want to find its coregulators. (default: None)'
    )

    parser.add_argument(
        '-a',
        '--annot',
        type=str,
        default=None,
        help='Annotation file for HOCOMOCO motifs (default: None)')

    parser.add_argument('-i',
                        '--inner',
                        type=int,
                        default=10,
                        help='Inner distance to target motif (default: 10)')

    parser.add_argument('-o',
                        '--outer',
                        type=int,
                        default=40,
                        help='Outer distance to target motif (default: 40)')

    parser.add_argument('-O',
                        '--output',
                        type=str,
                        default='output.csv',
                        help='Output CSV file (default: output.csv)')

    args = parser.parse_args()

    main(args)
