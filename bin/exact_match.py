#!/usr/bin/env python3
"""
exact_match.py - Find all exact substring matches between two sequences.
Drop-in replacement for blastn with similar interface and outfmt 6 output.

Usage:
    exact_match.py -query <query.fa> -db <target.fa> [-word_size <n>] [-outfmt 6] [-strand plus|minus|both]
"""

import argparse
import sys


def parse_fasta(filepath):
    """Parse a simple FASTA file, return list of (name, sequence) tuples."""
    sequences = []
    current_name = None
    current_seq = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name is not None:
                    sequences.append((current_name, ''.join(current_seq).upper()))
                current_name = line[1:].split()[0]  # Take first word after >
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_name is not None:
            sequences.append((current_name, ''.join(current_seq).upper()))
    
    return sequences


def reverse_complement(seq):
    """Return reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                  'N': 'N', 'n': 'n'}
    return ''.join(complement.get(base, base) for base in reversed(seq))


def find_exact_matches(query_seq, target_seq, min_length=4):
    """
    Find all exact substring matches between query and target.
    Returns list of (qstart, qend, sstart, send, length) tuples.
    Positions are 1-based, inclusive (like BLAST).
    
    Reports all matches with no positional redundancy:
    - Different matches CAN overlap (same query pos can match different target pos)
    - But we don't report a 7bp submatch contained within an 8bp match 
      at the same (qstart, tstart) position
    """
    matches = []
    qlen = len(query_seq)
    tlen = len(target_seq)
    
    # For each starting position in query
    for qstart in range(qlen):
        # For each starting position in target
        for tstart in range(tlen):
            # Check if sequences match at this position
            if query_seq[qstart] != target_seq[tstart]:
                continue
            
            # Check if this is a continuation of a match that started earlier
            # (i.e., both qstart-1 and tstart-1 also matched)
            if (qstart > 0 and tstart > 0 and 
                query_seq[qstart - 1] == target_seq[tstart - 1]):
                # This is a submatch, skip it
                continue
            
            # Find how long the match extends
            match_len = 0
            while (qstart + match_len < qlen and 
                   tstart + match_len < tlen and
                   query_seq[qstart + match_len] == target_seq[tstart + match_len]):
                match_len += 1
            
            # Record the maximal match if >= min_length
            if match_len >= min_length:
                # Convert to 1-based inclusive coordinates
                qend = qstart + match_len
                send = tstart + match_len
                matches.append((qstart + 1, qend, tstart + 1, send, match_len))
    
    return matches


def format_outfmt6(qname, sname, qstart, qend, sstart, send, length, strand='+'):
    """
    Format output like BLAST outfmt 6:
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    """
    pident = 100.0  # Always 100% for exact matches
    mismatch = 0
    gapopen = 0
    # For exact matches, evalue and bitscore are placeholders
    # Rough approximation: longer matches = lower evalue, higher bitscore
    evalue = 10 ** (-(length - 3))  # Very rough approximation
    bitscore = length * 2  # Rough approximation
    
    if strand == '-':
        # For minus strand, BLAST reports target coords in reverse
        sstart, send = send, sstart
    
    return f"{qname}\t{sname}\t{pident:.3f}\t{length}\t{mismatch}\t{gapopen}\t{qstart}\t{qend}\t{sstart}\t{send}\t{evalue:.2e}\t{bitscore:.1f}"


def main():
    parser = argparse.ArgumentParser(
        description='Find exact substring matches (blastn drop-in replacement)')
    parser.add_argument('-query', required=True, help='Query FASTA file')
    parser.add_argument('-db', required=True, help='Target/database FASTA file')
    parser.add_argument('-word_size', type=int, default=4, 
                        help='Minimum match length (default: 4)')
    parser.add_argument('-outfmt', type=str, default='6',
                        help='Output format (only 6 supported)')
    parser.add_argument('-strand', choices=['plus', 'minus', 'both'], 
                        default='plus', help='Strand to search')
    parser.add_argument('-out', type=str, default=None,
                        help='Output file (default: stdout)')
    
    args = parser.parse_args()
    
    # Parse input files
    query_seqs = parse_fasta(args.query)
    target_seqs = parse_fasta(args.db)
    
    # Prepare output
    output_lines = []
    
    # For each query-target pair
    for qname, qseq in query_seqs:
        for tname, tseq in target_seqs:
            # Plus strand
            if args.strand in ['plus', 'both']:
                matches = find_exact_matches(qseq, tseq, args.word_size)
                for qstart, qend, sstart, send, length in matches:
                    output_lines.append(
                        format_outfmt6(qname, tname, qstart, qend, sstart, send, length, '+')
                    )
            
            # Minus strand (reverse complement of target)
            if args.strand in ['minus', 'both']:
                tseq_rc = reverse_complement(tseq)
                matches = find_exact_matches(qseq, tseq_rc, args.word_size)
                for qstart, qend, sstart, send, length in matches:
                    # Convert coordinates back to original target orientation
                    # When we RC the target, position i becomes len-i-1
                    orig_sstart = len(tseq) - send + 1
                    orig_send = len(tseq) - sstart + 1
                    output_lines.append(
                        format_outfmt6(qname, tname, qstart, qend, orig_sstart, orig_send, length, '-')
                    )
    
    # Sort by length (descending), then by qstart
    def sort_key(line):
        parts = line.split('\t')
        return (-int(parts[3]), int(parts[6]))  # -length, qstart
    
    output_lines.sort(key=sort_key)
    
    # Output
    output_text = '\n'.join(output_lines)
    if args.out:
        with open(args.out, 'w') as f:
            f.write(output_text + '\n')
    else:
        print(output_text)


if __name__ == '__main__':
    main()
