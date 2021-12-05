#!/usr/bin/env python

import argparse
import os
import sys
from Bio.Blast import NCBIXML
from Bio import SeqIO


MAX_QUERY_DESCRIPTION_WORDS = 3
MAX_MATCH_DESCRIPTION_WORDS = 5
QUERY_EXT = '.fasta'
REPORT_SUFFIX = '_blastp_hits.tsv'


def info(msg):
    sys.stdout.write(str(msg) + '\n')
    sys.stdout.flush()


def warn(msg):
    sys.stderr.write('WARNING! ' + str(msg) + '\n')
    sys.stderr.flush()


def _truncate_words(phrase, max_words):
    return ' '.join(phrase.split()[:max_words])


class BlastHit(object):
    def __init__(self, blast_record, alignment, hsp, query_description=''):
        self.query_id = blast_record.query_id
        self.query_description = query_description if query_description else blast_record.query
        self.query_len = blast_record.query_length
        self.title = alignment.title
        self.description = alignment.hit_def
        self.accession = alignment.accession
        self.accession_link = 'https://www.ncbi.nlm.nih.gov/protein/{}?report=genbank'.format(self.accession)
        self.acc_length = alignment.length
        self.bit_score = hsp.bits
        self.query_cover = 100.0 * (hsp.align_length - hsp.gaps) / self.query_len
        self.e_value = hsp.expect
        self.per_idy = 100.0 * hsp.identities / hsp.align_length

    def __str__(self):
        '''
        standard BLAST output is:
            Description
            Scientific Name
            Max Score
            Total Score
            Query Cover
            E value
            Per. Ident
            Acc. Len
            Accession
        '''
        # return '\t'.join(map(str, [_truncate_words(self.query_description, MAX_QUERY_DESCRIPTION_WORDS),
        #                            _truncate_words(self.description, MAX_MATCH_DESCRIPTION_WORDS),
        #                            self.accession,
        #                            int(self.bit_score), '%.2f%%' % self.query_cover, self.e_value,
        #                            '%.2f%%' % self.per_idy]))
        #
        return '\t'.join(map(str, [_truncate_words(self.query_description, MAX_QUERY_DESCRIPTION_WORDS),
                                   '%.2f%%' % self.query_cover, '%.2f%%' % self.per_idy,
                                   int(self.bit_score), self.e_value,
                                   _truncate_words(self.description, MAX_MATCH_DESCRIPTION_WORDS),
                                   self.accession_link]))

    def is_good(self, filter_options):
        return self.e_value <= filter_options.e_value and self.query_cover >= filter_options.query_cover \
               and self.per_idy >= filter_options.per_idy

    def get_sort_key(self):
        return self.e_value, -self.bit_score

    @staticmethod
    def get_header_columns():
        # return '\t'.join(['Query Description',
        #                   'Match Description', 'Accession', 'Score (bits)', 'Query Cover', 'E-value', 'Per. Ident'])
        return ['Query Description', 'Query Cover', 'Per. Ident', 'Score (bits)', 'E-value',
                'Match Description', 'Accession Link']


def handle_single_xml(fpath, filter_options, query_fpath=None):
    query_description = ''
    if query_fpath is not None and os.path.isfile(query_fpath):
        ### Here we expect a single-entry FASTA as input
        for record in SeqIO.parse(query_fpath, "fasta"):
            query_description = record.description
            break

    with open(fpath) as result_handle:
        ### Here we expect a single query per XML, if there would be multi-entry queries, the following loop would be needed:
        # for blast_record in NCBIXML.parse(result_handle):
        blast_record = NCBIXML.read(result_handle)
        filtered_hits = []
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                hit = BlastHit(blast_record, alignment, hsp, query_description)
                if hit.is_good(filter_options):
                    filtered_hits.append(hit)
        filtered_hits.sort(key=lambda x:x.get_sort_key())
        return filtered_hits[:filter_options.num_top_hits]


def dir_contains_xml_files(dirpath):
    return any(map(lambda x: os.path.splitext(x)[1] == '.xml', os.listdir(dirpath)))


def summarise_dir(args, dirpath, report_fpath, query_dir=None):
    with open(report_fpath, 'w') as report:
        report.write('\t'.join(['QueryIdx', 'HitIdx'] + BlastHit.get_header_columns()) + '\n')

        query_idx = 0
        for fpath in os.listdir(dirpath):
            if os.path.splitext(fpath)[1] != '.xml':
                continue

            query_idx += 1
            query_path = None
            if query_dir is not None and os.path.isdir(query_dir):
                query_path = os.path.join(query_dir, os.path.splitext(fpath)[0] + QUERY_EXT)

            blast_xml_fullpath = os.path.join(dirpath, fpath)
            hits = handle_single_xml(blast_xml_fullpath, args, query_path)
            for idx, hit in enumerate(hits):
                hit_idx = idx + 1
                report.write('{}\t{}\t{}\n'.format(query_idx, hit_idx, str(hit)))

    info('Report created: {}'.format(report_fpath))
    return


def main():
    parser = argparse.ArgumentParser(description='Parser of BLAST (blastp) results in the XML format. '
                                                 'Summarizes all hits in a single TSV file.')
    parser.add_argument(
        'inputs',
        metavar='DIR',
        type=str,
        nargs='+',
        help='paths to dir with blastp XML files or dir containing several subdirs with such files'
    )
    parser.add_argument(
        '-o', '--output-dir',
        default='.',
        help='Output dir for summary reports (TSV files)'
    )
    parser.add_argument(
        '--query-root-dir',
        default=None,
        help='Path to root dir with original queries in the FASTA format (to extract descriptions)'
    )
    parser.add_argument(
        '-t', '--num-top-hits',
        default=3,
        type=int,
        help='Number of top BLAST hits per entry to report'
    )
    parser.add_argument(
        '-e', '--e-value',
        default=0.001,
        type=float,
        help='Max E-value to report (0.0 - 1.0)'
    )
    parser.add_argument(
        '-c', '--query-cover',
        default=50.0,
        type=float,
        help='Min Query Coverage to report, in %%'
    )
    parser.add_argument(
        '-i', '--per-idy',
        default=50.0,
        type=float,
        help='Min Percent Identity to report, in %%'
    )

    args = parser.parse_args()
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)

    for input_path in args.inputs:
        if not os.path.isdir(input_path):
            warn('{} is not a directory! Skipping'.format(input_path))
            continue
        if dir_contains_xml_files(input_path):
            info('{} contains XML files, summarising their hits!'.format(input_path))
            report_fpath = os.path.join(args.output_dir, os.path.basename(input_path) + REPORT_SUFFIX)
            query_dir = None if args.query_root_dir is None else args.query_root_dir
            summarise_dir(args, input_path, report_fpath, query_dir)
        else:
            info('{} does not contain XML files, processing its subdirectories independently!'.format(input_path))
            has_valid_subdir = False
            for subdir_basename in os.listdir(input_path):
                subdir_fullpath = os.path.join(input_path, subdir_basename)
                if not os.path.isdir(subdir_fullpath):
                    continue
                if not dir_contains_xml_files(subdir_fullpath):
                    continue
                has_valid_subdir = True
                report_fpath = os.path.join(args.output_dir, subdir_basename + REPORT_SUFFIX)
                query_dir = None if args.query_root_dir is None else os.path.join(args.query_root_dir, subdir_basename)
                summarise_dir(args, subdir_fullpath, report_fpath, query_dir)
            if not has_valid_subdir:
                warn('{} does not contain valid subdirs as well! Skipping'.format(input_path))

    info("Done! All results are saved in {}".format(args.output_dir))


if __name__ == "__main__":
    main()
