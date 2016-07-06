import random

def split_pat_file(pat):
"""Takes pat file location [X].pat and writes two new pat files [X]-A.pat and [X]-B.pat each containing the original header and a randomized half of the original data"""
    original_pat = open(pat, 'r')
    first_half = open(pat[:-4]+ "-A.pat", 'w')
    second_half = open(pat[:-4] + "-B.pat", 'w')
    all_lines = original_pat.readlines()
    numbers = len(all_lines) -2
    first_half.write(all_lines[0])
    first_half.write(all_lines[1])
    second_half.write(all_lines[0])
    second_half.write(all_lines[1])
    
    lines = list(range(numbers))
    random.shuffle(lines)
    for line in lines[:numbers/2]:
        first_half.write(all_lines[line+2])
    for line in lines[numbers/2:]:
        second_half.write(all_lines[line+2])
    return((pat[:-4]+"-A.pat", pat[:-4] + "-B.pat"))
