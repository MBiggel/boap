#!/usr/bin/env python3

import argparse
import sys
import os
import multiprocessing
from collections import defaultdict
from scipy.stats import binomtest

VERSION = "0.1.0"

# ==============================================================================
# CONFIGURATION PARAMETERS
# ==============================================================================
LQB_THRESHOLD = 5           # Phred Q-score threshold (inclusive)
MIN_CONTIG_LEN = 50000      # Minimum contig length to include in analysis
LQB_WINDOW_SIZE = 5000      # Window size (bp) for scanning LQB dense regions
DENSITY_MULTIPLIER = 5.0    # Multiplier for baseline density to trigger masking
MASK_FLOOR = 0.001          # Absolute minimum LQB density (0.001 = 0.1%)
# ==============================================================================

def get_significance_stars(p_value, num_tests):
    if num_tests == 0: return ""
    alpha_1, alpha_2, alpha_3 = 0.05/num_tests, 0.01/num_tests, 0.001/num_tests
    if p_value < alpha_3: return "***"
    if p_value < alpha_2: return "**"
    if p_value < alpha_1: return "*"
    return ""

def identify_lqb_dense_regions(scores, lqb_threshold, window_size, dynamic_density_threshold):
    seq_len = len(scores)
    indices_to_mask = set()
    
    is_lqb = [1 if 1 <= s <= lqb_threshold else 0 for s in scores]
    
    if seq_len < window_size:
        if (sum(is_lqb) / seq_len) > dynamic_density_threshold:
            return set(range(seq_len))
        return set()

    step = int(window_size * 0.1) if window_size > 10 else 1
    
    for start in range(0, seq_len - window_size + 1, step):
        end = start + window_size
        window_lqb_count = sum(is_lqb[start:end])
        if (window_lqb_count / window_size) > dynamic_density_threshold:
            for i in range(start, end):
                indices_to_mask.add(i)
    
    tail_start = max(0, seq_len - window_size)
    tail_end = seq_len
    tail_lqb_count = sum(is_lqb[tail_start:tail_end])
    tail_len = tail_end - tail_start
    
    if tail_len > 0 and (tail_lqb_count / tail_len) > dynamic_density_threshold:
        for i in range(tail_start, tail_end):
            indices_to_mask.add(i)
                
    return indices_to_mask

def analyze_single_file(file_path, args):
    file_total_bases = 0 
    file_total_lqb = 0   
    analyzed_base_count = 0 
    lqb_filtered_count = 0 
    absolute_total_q_sum = 0 
    scope_total_bases = 0
    total_masked_bases = 0
    
    k_list = [4, 5, 6]
    kmer_stats = {k: defaultdict(lambda: [0, 0]) for k in k_list}
    
    records = []
    
    try:
        with open(file_path, 'r') as f:
            while True:
                header = f.readline()
                if not header: break
                sequence = f.readline().strip().upper()
                f.readline() 
                quality_line = f.readline().strip()
                
                if len(sequence) == len(quality_line):
                    slen = len(sequence)
                    q_scores = [ord(char) - 33 for char in quality_line]
                    records.append({'head': header, 'seq': sequence, 'qual': quality_line, 'len': slen, 'q_scores': q_scores})
                    file_total_bases += slen
                    
                    file_total_lqb += sum(1 for s in q_scores if 1 <= s <= LQB_THRESHOLD)
                    absolute_total_q_sum += sum(q_scores)

        if not records:
             return [os.path.basename(file_path), "NoData", "0", "0", "0", "0", "0", "0", "", "", ""]

        records.sort(key=lambda x: x['len'], reverse=True)
        longest_rec = records[0]
        longest_len = longest_rec['len']
        
        longest_lqb = sum(1 for s in longest_rec['q_scores'] if 1 <= s <= LQB_THRESHOLD)
        
        raw_density = longest_lqb / longest_len if longest_len > 0 else 0
        mask_threshold_used = max(raw_density * DENSITY_MULTIPLIER, MASK_FLOOR)

        data_to_analyze = records if args.all_contigs else [records[0]]

        for rec in data_to_analyze:
            sequence = rec['seq']
            scores = rec['q_scores']
            seq_len = rec['len']

            if seq_len < MIN_CONTIG_LEN:
                continue
            
            scope_total_bases += seq_len
            indices_to_mask = set()
            if not args.no_mask:
                indices_to_mask = identify_lqb_dense_regions(scores, LQB_THRESHOLD, LQB_WINDOW_SIZE, mask_threshold_used)
            
            total_masked_bases += len(indices_to_mask)

            for i, score in enumerate(scores):
                if i in indices_to_mask: continue 
                analyzed_base_count += 1
                
                if 1 <= score <= LQB_THRESHOLD:
                    lqb_filtered_count += 1

            for k in k_list:
                for i in range(seq_len - k + 1):
                    if any((i + offset) in indices_to_mask for offset in range(k)):
                        continue
                    kmer_raw = sequence[i : i+k]
                    if 'N' in kmer_raw: continue
                    
                    kmer = kmer_raw
                    win_scores = scores[i : i+k]
                    kmer_stats[k][kmer][0] += 1 
                    
                    if any(1 <= s <= LQB_THRESHOLD for s in win_scores):
                        kmer_stats[k][kmer][1] += 1 

        avg_q = absolute_total_q_sum / file_total_bases if file_total_bases > 0 else 0
        lqb_raw_str = f"{(file_total_lqb / file_total_bases) * 1_000_000:.2f}" if file_total_bases > 0 else "0.00"
        lqb_filtered_str = f"{(lqb_filtered_count / analyzed_base_count) * 1_000_000:.2f}" if analyzed_base_count > 0 else "Low Data"
        masked_str = f"{total_masked_bases}({(total_masked_bases / scope_total_bases) * 100:.1f}%)" if scope_total_bases > 0 else "0(0.0%)"
        analyzed_pct = (analyzed_base_count / file_total_bases * 100) if file_total_bases > 0 else 0
        contigs_str = f"{len(data_to_analyze)}//{len(records)}"

        if analyzed_base_count > 0:
            k_summaries = {}
            for k in k_list:
                results = []
                for kmer, (tot, broken) in kmer_stats[k].items():
                    if broken == 0: continue
                    p_val = binomtest(broken, tot, 0.0001, alternative='greater').pvalue
                    sig = get_significance_stars(p_val, len(kmer_stats[k]))
                    if sig: 
                        results.append({'kmer': kmer, 'tot': tot, 'broken': broken, 'prop': (broken/tot), 'sig': sig})
                
                results.sort(key=lambda x: x['prop'], reverse=True)
                k_summaries[k] = ",".join([f"{r['kmer']}({r['prop']:.1%})" for r in results[:3]]) if results else "None"

                if args.report:
                    report_name = f"{os.path.splitext(os.path.basename(file_path))[0]}_report.txt"
                    with open(report_name, 'a') as rf:
                        rf.write(f"\nMotif Analysis {k}-mers (Filtered Region)\n" + "-"*40 + "\n")
                        for r in results[:15]: 
                            rf.write(f"{r['kmer']} | Hits:{r['broken']}/{r['tot']} | Rate:{r['prop']:.1%} | {r['sig']}\n")
        else:
            k_summaries = {4: "Low Data", 5: "Low Data", 6: "Low Data"}

        return [
            os.path.basename(file_path), f"{avg_q:.2f}", lqb_raw_str, lqb_filtered_str,
            f"{mask_threshold_used:.4%}", masked_str, contigs_str,      
            f"{analyzed_base_count}({analyzed_pct:.1f}%)", 
            k_summaries[4], k_summaries[5], k_summaries[6]
        ]

    except Exception as e:
        return [os.path.basename(file_path), "ERROR", str(e), "", "", "", "", "", "", "", ""]

def main():
    parser = argparse.ArgumentParser(
        description=f"ALPAQA v{VERSION}", 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("-i", "--input", nargs="+", required=True, help="Input fastq assembly files.")
    parser.add_argument("-o", "--output", default="alpaqa_report.tsv", help="Output TSV filename.")
    parser.add_argument("--report", action="store_true", help="Generate detailed reports on kmer analysis.")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads.")
    parser.add_argument("--all-contigs", action="store_true", help="Analyze ALL contigs > MIN_CONTIG_LEN (Default is longest contig only).")
    parser.add_argument("--no-mask", action="store_true", help="Disable masking of LQB dense regions.")
    parser.add_argument("-v", "--version", action="version", version=f"ALPAQA {VERSION}")

    args = parser.parse_args()
    headers = ["Filename", "AvgQ", "LQB_raw/Mbp", "LQB/Mbp", "MaskThresh", "Bases_Masked", "Contigs_Analyzed", "Bases_Analyzed", "Sig4m", "Sig5m", "Sig6m"]
    all_rows = []
    
    if args.threads > 1:
        pool = multiprocessing.Pool(processes=args.threads)
        tasks = [pool.apply_async(analyze_single_file, (f, args)) for f in args.input]
        for t in tasks:
            all_rows.append(t.get())
        pool.close()
        pool.join()
    else:
        for f in args.input:
            all_rows.append(analyze_single_file(f, args))

    with open(args.output, 'w') as f:
        f.write("\t".join(headers) + "\n")
        for row in all_rows:
            f.write("\t".join(map(str, row)) + "\n")

if __name__ == "__main__":
    main()