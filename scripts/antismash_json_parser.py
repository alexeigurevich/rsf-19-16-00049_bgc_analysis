#!/usr/bin/env python

import argparse
import os
import sys
import json
from collections import OrderedDict
# from matplotlib import pyplot as plt
# import numpy as np


MIN_FULL_REGION_LENGTH = 5000
# The list below is originally from antiSMASH v.5, but newer names (e.g., from v.6)
# will be automatically added to the report if at least one BGC with such product is found.
# The most up-to-date list is here: https://docs.antismash.secondarymetabolites.org/glossary/#clustertypes
ALL_PRODUCTS = ['T1PKS', 'T2PKS', 'T3PKS', 'transAT-PKS', 'transAT-PKS-like', 'PpyS-KS', 'hglE-KS', 'CDPS', 'PKS-like',
                'arylpolyene', 'resorcinol', 'ladderane', 'PUFA', 'NRPS', 'NRPS-like', 'thioamide-NRP', 'terpene',
                'lanthipeptide', 'lipolanthine', 'bacteriocin', 'betalactone', 'thiopeptide', 'linaridin',
                'cyanobactin', 'glycocin', 'LAP', 'lassopeptide', 'sactipeptide', 'bottromycin', 'head_to_tail',
                'microviridin', 'proteusin', 'blactam', 'amglyccycl', 'aminocoumarin', 'siderophore', 'ectoine',
                'butyrolactone', 'indole', 'nucleoside', 'phosphoglycolipid', 'melanin', 'oligosaccharide',
                'furan', 'hserlactone', 'phenazine', 'phosphonate', 'fused', 'PBDE', 'acyl_amino_acids',
                'tropodithietic-acid', 'NAGGN', 'RaS-RiPP', 'fungal-RiPP', 'TfuA-related', 'saccharide',
                'fatty_acid', 'halogenated', 'other']
# NRPS_PRODUCTS = ['NRPS', 'NRPS-like', 'thioamide-NRP']
# PKS_PRODUCTS = ['T1PKS', 'T2PKS', 'T3PKS', 'transAT-PKS', 'transAT-PKS-like', 'PKS-like', 'hglE-KS']
# RiPP_PRODUCTS = ['']
# MIXED_PRODUCTS = ['NRPS-PKS', 'NRPS-xor-PKS', 'not-NRPS-PKS-mix']


def _is_NRPS(product):
    return 'NRP' in product


def _is_PKS(product):
    return 'PKS' in product or 'hglE-KS' in product  # hglE-KS is a special case: "heterocyst glycolipid synthase-like PKS"


def _is_terpene(product):
    return 'terpene' in product


def _is_RiPP(product):
    for ripp_kind in ['RiPP', 'lanthipeptide', 'lipolanthine', 'bacteriocin', 'thiopeptide', 'linaridin', 'cyanobactin',
                      'glycocin', 'LAP', 'lassopeptide', 'sactipeptide', 'bottromycin', 'microviridin', 'proteusin']:
        if ripp_kind in product:
            return True
    return False


MIX_PREFIX = 'mix:'
MIX_SEPARATOR = ','


def info(msg, verbose=True):
    if verbose:
        print(msg)


def error(msg, exit=True):
    print('ERROR! ' + msg)
    if exit:
        sys.exit(1)


def __parse_location(location):
    # e.g. 'location' = '[351:486](+)'
    coords = location.split(']')[0][1:]
    start, end = map(int, coords.split(':'))
    if '(' in location:
        strand = location.split('](')[1][0]
    else:
        strand = ''
    return start, end, strand


def handle_single_input(path, cluster_type_to_print='none',
                        output_dir='', save_seq_mode='', save_gene_kind='',
                        skip_partial_regions=False):
    info('\nProcessing ' + path)
    path = os.path.abspath(path)
    main_json_path = path
    if os.path.isdir(path):
        # TODO: if there is a single JSON inside the dir, just take it independently of its name
        main_json_path = os.path.join(path, os.path.basename(os.path.normpath(path)) + '.json')
    if not os.path.isfile(main_json_path):
        error('Main antiSMASH v.5 JSON file not found: %s. Skipping this input..' % main_json_path)
        return None, None

    run_name = os.path.splitext(os.path.basename(main_json_path))[0]
    info('Processing JSON %s (run: %s)' % (main_json_path, run_name))
    need_to_print_seqs = True if output_dir else False
    if need_to_print_seqs:
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        nucleotide_sequences = [] if 'n' in save_seq_mode else None
        protein_sequences = [] if 'p' in save_seq_mode else None

    with open(main_json_path, 'r') as f:
        data = json.load(f)

    regions_per_product = dict()  # 'product_name': (num_full, num_partial)
    region_ids_per_product = []
    region_lengths = []
    full_region_lengths = []
    suitable_gene_kinds = []
    if 'add' in save_gene_kind:
        suitable_gene_kinds.append('biosynthetic-additional')
    if 'core' in save_gene_kind:
        suitable_gene_kinds.append('biosynthetic')
    for contig_data in data["records"]:
        product_name = 'undefined product'
        current_region_is_partial = False
        for feature in contig_data['features']:
            if feature['type'] == 'region':
                if len(feature['qualifiers']['candidate_cluster_numbers']) > 1:
                    pass  # ignore them for know -- let's consider such cases as a simple mixed cluster (>1 product)
                if len(feature['qualifiers']['product']) == 1:
                    product_name = feature['qualifiers']['product'][0]
                else:
                    product_name = MIX_PREFIX + MIX_SEPARATOR.join(feature['qualifiers']['product'])
                if product_name not in regions_per_product:
                    regions_per_product[product_name] = [0, 0]

                is_on_contig_edge = feature['qualifiers']['contig_edge'][0] == 'True'
                start, end, _ = __parse_location(feature['location'])
                region_length = end - start  # "+ 1" is not needed since the coordinates here are 0-based
                region_lengths.append(region_length)
                if is_on_contig_edge or region_length < MIN_FULL_REGION_LENGTH:
                    regions_per_product[product_name][1] += 1   # partial cluster
                    current_region_is_partial = True
                else:
                    regions_per_product[product_name][0] += 1   # full cluster
                    full_region_lengths.append(region_length)
                    current_region_is_partial = False

                    # we want to print only full clusters
                    if product_name == cluster_type_to_print:  # if cluster_type_to_print == 'none' this is simply False
                        region_ids_per_product.append(contig_data['id'])
            elif need_to_print_seqs and feature['type'] == 'CDS':
                if skip_partial_regions and current_region_is_partial:
                    continue
                if 'ripp' in save_gene_kind and not _is_RiPP(product_name):
                    continue

                if 'gene_kind' in feature['qualifiers'] and feature['qualifiers']['gene_kind'][0] in suitable_gene_kinds:
                    assert len(feature['qualifiers']['locus_tag']) == 1, "cannot work with multiple locus tags, yet"
                    locus_tag = feature['qualifiers']['locus_tag'][0]
                    gene_kind = feature['qualifiers']['gene_kind'][0]
                    location = feature['location']
                    start, end = location[1:-4].split(':')
                    strand = location[-2]
                    record_id = '_'.join([locus_tag, start, end, 'pos' if strand == '+' else 'neg'])
                    record_description = product_name + ' ' + gene_kind + ' gene'
                    if nucleotide_sequences is not None:
                        nucl_seq = Seq(contig_data['seq']['data'][int(start):int(end)])
                        if strand == '-':
                            nucl_seq = nucl_seq.reverse_complement()
                        nucl_record = SeqRecord(nucl_seq, id=record_id, description=record_description)
                        nucleotide_sequences.append(nucl_record)
                    if protein_sequences is not None:
                        prot_record = SeqRecord(Seq(feature['qualifiers']['translation'][0]),
                                                id=record_id, description=record_description)
                        protein_sequences.append(prot_record)

    print(regions_per_product)
    # print(len(region_lengths))
    # region_lengths.sort(reverse=True)
    # a = np.array(region_lengths)
    # plt.hist(a, bins='auto')
    # plt.title("histogram")
    # plt.show()
    # a = np.array(full_region_lengths)
    # plt.hist(a, bins='auto')
    # plt.title("histogram full")
    # plt.show()
    # hist, bins = np.histogram(a, bins=[0, 10000, 50000, 100000, 500000])
    # print(hist)
    # print(bins)

    if region_ids_per_product:
        print('Contigs with full clusters of interest (%s):' % cluster_type_to_print)
        print(region_ids_per_product)

    if need_to_print_seqs:
        output_dir = os.path.abspath(output_dir)
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        from Bio import SeqIO
        if 'm' in save_seq_mode:  # 'merged mode': single multi-FASTA file
            if nucleotide_sequences is not None:
                nucleotide_out = os.path.join(output_dir, run_name + '_nucl.fasta')
                with open(nucleotide_out, 'w') as output_handle:
                    SeqIO.write(nucleotide_sequences, output_handle, 'fasta')
                info('Saved nucl sequences of %s genes to %s under %s' % (save_gene_kind, os.path.basename(nucleotide_out), output_dir))

            if protein_sequences is not None:
                protein_out = os.path.join(output_dir, run_name + '_prot.fasta')
                with open(protein_out, 'w') as output_handle:
                    SeqIO.write(protein_sequences, output_handle, 'fasta')
                info('Saved prot sequences of %s genes to %s under %s' % (save_gene_kind, os.path.basename(protein_out), output_dir))

        if 's' in save_seq_mode:  # 'split mode': many single-FASTA files
            if nucleotide_sequences is not None:
                base_output_dir = os.path.join(output_dir, run_name + '_nucl')
                if not os.path.isdir(base_output_dir):
                    os.makedirs(base_output_dir)
                for seq in nucleotide_sequences:
                    nucleotide_out = os.path.join(base_output_dir, seq.id + '.fasta')
                    with open(nucleotide_out, 'w') as output_handle:
                        SeqIO.write(seq, output_handle, 'fasta')
                info('Saved nucl sequences of %s genes under %s' % (save_gene_kind, base_output_dir))

            if protein_sequences is not None:
                base_output_dir = os.path.join(output_dir, run_name + '_prot')
                if not os.path.isdir(base_output_dir):
                    os.makedirs(base_output_dir)
                for seq in protein_sequences:
                    protein_out = os.path.join(base_output_dir, seq.id + '.fasta')
                    with open(protein_out, 'w') as output_handle:
                        SeqIO.write(seq, output_handle, 'fasta')
                info('Saved prot sequences of %s genes under %s' % (save_gene_kind, base_output_dir))

    return run_name, regions_per_product


def print_summary_tables(summary, product_names, skip_empty_lines=True, skip_partial_regions=False):
    separator = ','
    full_and_partial_separator = ' + '
    empty_value = [0, 0]
    header = ['RUN'] + list(summary.keys())

    full_values = []
    included_product_names = []
    mixed_values = OrderedDict([('NRPS-PKS', [[0, 0] for _ in summary.keys()]),
                                ('NRPS-xor-PKS', [[0, 0] for _ in summary.keys()]),
                                ('not-NRPS-PKS-mix', [[0, 0] for _ in summary.keys()])])

    clustered_values = OrderedDict([('NRPS', [[0, 0] for _ in summary.keys()]),
                                    ('PKS', [[0, 0] for _ in summary.keys()]),
                                    ('NRPS or/and PKS hybrids', [[0, 0] for _ in summary.keys()]),
                                    ('terpene', [[0, 0] for _ in summary.keys()]),
                                    ('RiPP', [[0, 0] for _ in summary.keys()]),
                                    ('OTHER', [[0, 0] for _ in summary.keys()])])

    def __increment(d, key, addition, idx=None):
        if idx is None:  # special case: update the full line
            pass
            for idx in range(len(d[key])):
                for i, v in enumerate(addition[idx]):
                    d[key][idx][i] += v
        else:
            for i, v in enumerate(addition):
                d[key][idx][i] += v

    for product in product_names:
        if product.startswith(MIX_PREFIX):  # special case: hybrids. They are processed below in a separate loop
            continue
        values = []
        for run_products in summary.values():
            if product in run_products:
                values.append(run_products[product])
            else:
                values.append(empty_value)
        if skip_empty_lines and values.count(empty_value) == len(values):
            continue
        if _is_NRPS(product):
            __increment(clustered_values, 'NRPS', values)
        elif _is_PKS(product):
            __increment(clustered_values, 'PKS', values)
        elif _is_terpene(product):
            __increment(clustered_values, 'terpene', values)
        elif _is_RiPP(product):
            __increment(clustered_values, 'RiPP', values)
        else:
            __increment(clustered_values, 'OTHER', values)
        full_values.append(values)
        included_product_names.append(product)

    # special case: mixes
    # format in the dict: product_name = 'mix:' + ','.join(feature['qualifiers']['product'])
    # mixed_values = OrderedDict([('NRPS-PKS', [[0, 0] for _ in summary.keys()]),
    #                             ('NRPS-xor-PKS', [[0, 0] for _ in summary.keys()]),
    #                             ('not-NRPS-PKS-mix', [[0, 0] for _ in summary.keys()])])
    for run_name, run_products in summary.items():
        run_index = list(summary.keys()).index(run_name)
        for product, value in run_products.items():
            if product.startswith(MIX_PREFIX):
                has_PKS = False
                has_NRPS = False
                for sub_product in product[len(MIX_PREFIX):].split(MIX_SEPARATOR):
                    if _is_NRPS(sub_product):
                        has_NRPS = True
                    elif _is_PKS(sub_product):
                        has_PKS = True
                if has_PKS and has_NRPS:
                    __increment(mixed_values, 'NRPS-PKS', value, idx=run_index)
                elif not has_PKS and not has_NRPS:
                    __increment(mixed_values, 'not-NRPS-PKS-mix', value, idx=run_index)
                else:
                    __increment(mixed_values, 'NRPS-xor-PKS', value, idx=run_index)

    __increment(clustered_values, 'NRPS or/and PKS hybrids', mixed_values['NRPS-PKS'])
    __increment(clustered_values, 'NRPS or/and PKS hybrids', mixed_values['NRPS-xor-PKS'])
    __increment(clustered_values, 'OTHER', mixed_values['not-NRPS-PKS-mix'])

    def __print_line(stream, row_name, values, is_header=False):
        def __format_cell_value(value):
            return full_and_partial_separator.join(map(str, value[0:1] if skip_partial_regions else value))

        if is_header:
            stream.write(separator.join([row_name] + values) + '\n')
        else:
            stream.write(separator.join([row_name] + list(map(__format_cell_value, values))) + '\n')

    # with open(fpath, 'w') as f:
    with sys.stdout as f:
        # table full
        __print_line(f, 'RUN', list(summary.keys()), is_header=True)
        for i, product in enumerate(included_product_names):
            __print_line(f, product, full_values[i])
        for product, values in mixed_values.items():
            __print_line(f, product, values)

    # with sys.stdout as f:
        # table clustered
        f.write("== LARGE CLASS GROUPED TABLE ==\n")
        __print_line(f, 'RUN', list(summary.keys()), is_header=True)
        for product, values in clustered_values.items():
            __print_line(f, product, values)


def main():
    parser = argparse.ArgumentParser(description='Parser of antiSMASH v5 output (main JSON files) and '
                                                 'generator of a summary report on found BGCs')
    parser.add_argument(
        'inputs',
        metavar='FILE/DIR',
        type=str,
        nargs='+',
        help='paths to antiSMASH v.5 JSONs or its output directories containing JSONs'
    )
    parser.add_argument(
        '-p', '--print-full-clusters-info',
        dest='cluster_type_to_print',
        default='none',
        help='Print IDs of full clusters for the specified product type ("all" for all types, "none" for not printing)'
        # TODO: 'all' is not implemented yet
    )
    parser.add_argument(
        '-f', '--full-clusters-only',
        dest='skip_partial_regions',
        action='store_true',
        default=False,
        help='Skip partial (short or on contig edge) clusters, a.k.a. "regions", in summary count tables '
             'and in core gene sequence saving (if requested)'
    )
    parser.add_argument(
        '-s', '--save-gene-sequences',
        dest='output_dir',
        default='',
        help='Save nucleotide and/or protein sequences of BGC core and/or additional biosynthetic genes in '
             'the specified directory (don\'t specify to disable, you can choose specific save mode with --save-mode '
             'and kind of BGC genes to save with --save-gene-kind) '
             '[requires biopython!]'
    )
    parser.add_argument(
        '--save-mode',
        dest='save_seq_mode',
        default='npm',
        choices=['ns', 'ps', 'nps', 'nm', 'pm', 'npm'],
        help='Save nucleotide (n_), or protein (p_), or both types (np_) of sequences of core genes in the '
             'specified directory either split into many single-FASTA files (_s), one seq per file, or merged together '
             'into a large multi-FASTA file (_m). '
             '(used only if -s/--save-gene-sequences is specified)'
    )
    parser.add_argument(
        '--save-gene-kind',
        dest='save_gene_kind',
        default='core',
        choices=['core', 'add', 'core-and-add', 'ripp-add'],
        help='Save core (core), or additional (add), or both types (core-and-add) of biosynthetic genes '
             'for all types of BGCs or only additional biosynthetic genes of RiPP and RiPP-like BGCs (ripp-add)'
             '(used only if -s/--save-gene-sequences is specified)'
    )


    args = parser.parse_args()
    summary = OrderedDict()
    identified_products = set()
    for input_path in args.inputs:
        run_name, regions_per_product = handle_single_input(input_path, args.cluster_type_to_print,
                                                            args.output_dir, args.save_seq_mode, args.save_gene_kind,
                                                            skip_partial_regions=args.skip_partial_regions)
        if run_name is None:
            continue
        run_name_for_report = run_name
        i = 2
        while run_name_for_report in summary:
            run_name_for_report = run_name + '_' + str(i)
            i += 1
        summary[run_name_for_report] = regions_per_product
        identified_products |= set(regions_per_product.keys())

    product_names_for_report = sorted(list(identified_products | set(ALL_PRODUCTS)))
    print_summary_tables(summary, product_names_for_report, skip_partial_regions=args.skip_partial_regions)


if __name__ == "__main__":
    main()
