#!/usr/bin/env python
# coding: utf-8

# In[1]:


# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 22:27:27 2023

@author: MG
"""
import bisect
import numpy as np
import streamlit as st
from Bio import SeqIO
import io

# read fasta files
def read_fasta_file(uploaded_file):
    sequences = []
    content = uploaded_file.read().decode("utf-8")
    for record in SeqIO.parse(io.StringIO(content), "fasta"):
        sequences.append(str(record.seq))
    return sequences


# kmp table builiding
def build_kmp_table(pattern):
    table = [0] * len(pattern)
    j = 0
    for i in range(1, len(pattern)):
        while j > 0 and pattern[i] != pattern[j]:
            j = table[j - 1]
        if pattern[i] == pattern[j]:
            j += 1
        table[i] = j
    return table

#kmp searching 
def kmp_search(text, pattern):
    text = text.upper()
    pattern = pattern.upper()
    m, n = len(pattern), len(text)
    kmp_table = build_kmp_table(pattern)
    i = j = 0
    indices = []

    while i < n:
        if pattern[j] == text[i]:
            i += 1
            j += 1

            if j == m:
                indices.append(i - j)
                j = kmp_table[j - 1]
        else:
            if j != 0:
                j = kmp_table[j - 1]
            else:
                i += 1

    return indices
# approximate matching using levenshtein_distance algs
def levenshtein_distance(s1, s2):
    m, n = len(s1), len(s2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    for i in range(m + 1):
        for j in range(n + 1):
            if i == 0:
                dp[i][j] = j
            elif j == 0:
                dp[i][j] = i
            elif s1[i - 1] == s2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1]
            else:
                dp[i][j] = 1 + min(dp[i - 1][j],      # Deletion
                                   dp[i][j - 1],      # Insertion
                                   dp[i - 1][j - 1])  # Substitution

    return dp[m][n]

def approximate_match(text, pattern, max_distance):
    text = text.upper()
    pattern = pattern.upper()
    indices = []

    for i in range(len(text) - len(pattern) + 1):
        window = text[i:i + len(pattern)]
        distance = levenshtein_distance(window, pattern)
        if distance <= max_distance:
            indices.append(i)

    return indices

# CG count
def CG_count(seq):
    seq = seq.upper()
    length = len(seq)
    c_count = seq.count("C")
    G_count = seq.count("G")
    total = c_count + G_count
    return total / length

# reversed_seq
def reversed_seq(s):
    size = len(s)
    rev = ""
    for i in range(size - 1, -1, -1):
        rev += s[i]
    return rev.upper()

# complement function
def complement(seq):
    seq = seq.upper()
    cm = ""
    for i in range(len(seq)):
        if seq[i] == "A":
            cm += "T"
        elif seq[i] == "T":
            cm += "A"
        elif seq[i] == "C":
            cm += "G"
        elif seq[i] == "G":
            cm += "C"
    return cm

# AMINO ACIDS
def Translation_Table(seq):
    seq = seq.upper()
    dic = {
        "TTT": "F", "CTT": "L", "ATT": "I", "GTT": "V",
        "TTC": "F", "CTC": "L", "ATC": "I", "GTC": "V",
        "TTA": "L", "CTA": "L", "ATA": "I", "GTA": "V",
        "TTG": "L", "CTG": "L", "ATG": "M", "GTG": "V",
        "TCT": "S", "CCT": "P", "ACT": "T", "GCT": "A",
        "TCC": "S", "CCC": "P", "ACC": "T", "GCC": "A",
        "TCA": "S", "CCA": "P", "ACA": "T", "GCA": "A",
        "TCG": "S", "CCG": "P", "ACG": "T", "GCG": "A",
        "TAT": "Y", "CAT": "H", "AAT": "N", "GAT": "D",
        "TAC": "Y", "CAC": "H", "AAC": "N", "GAC": "D",
        "TAA": "*", "CAA": "Q", "AAA": "K", "GAA": "E",
        "TAG": "*", "CAG": "Q", "AAG": "K", "GAG": "E",
        "TGT": "C", "CGT": "R", "AGT": "S", "GGT": "G",
        "TGC": "C", "CGC": "R", "AGC": "S", "GGC": "G",
        "TGA": "*", "CGA": "R", "AGA": "R", "GGA": "G",
        "TGG": "W", "CGG": "R", "AGG": "R", "GGG": "G"
    }
    s = ""
    for i in range(0, len(seq)-2, 3):
        s += dic[seq[i:i+3]]
    return s

# naive match
def naive_match(seq, sub_seq):
    seq = seq.upper()
    sub_seq = sub_seq.upper()
    x = -1
    index = []
    for i in range(len(seq)):
        if sub_seq == seq[i:i+len(sub_seq)]:
            x = i
            index.append(x)
    return index

# Bad character
def Badchars(seq, sub_seq):
    ind = []
    seq = seq.upper()
    sub_seq = sub_seq.upper()
    table = np.zeros([4, len(sub_seq)])
    row = ["A", "C", "G", "T"]
    for i in range(4):
        num = -1
        for j in range(len(sub_seq)):
            if row[i] == sub_seq[j]:
                table[i, j] = -1
                num = -1
            else:
                num += 1
                table[i, j] = num
    x = -1
    i = 0
    while(i < len(seq)-len(sub_seq)+1):
        if sub_seq == seq[i:i+len(sub_seq)]:
            x = i
            ind.append(x)
            # break
        else:
            for j in range(len(sub_seq)-1, -1, -1):
                if seq[i+j] != sub_seq[j]:
                    k = row.index(seq[i+j])
                    i += table[k, j]
                    break
        i = int(i+1)
    return ind

# indexing
def IndexSorted(t, k):
    t = t.upper()
    tl = len(t)
    index = []
    steps = tl - k + 1
    for i in range(0, steps+1):
        pp = t[i: i+k]
        index.append((pp, i))
    index.sort()
    return index

def indexing(t, p, num):
    t = t.upper()
    p = p.upper()
    num = int()
    index = IndexSorted(t, num)
    keys = [r[0] for r in index]
    st = bisect.bisect_left(keys, p[:len(keys[0])])
    en = bisect.bisect(keys, p[:len(keys[0])])
    hits = index[st:en]
    l = [h[1] for h in hits]
    offsets = []
    for i in l:
        if t[i:i+len(p)] == p:
            offsets.append(i)
    return offsets

# suffix array
def suffix(T):
    T = T.upper()
    dec = {
        '$': 0,
        'A': 1,
        'C': 2,
        'G': 3,
        'T': 4
    }
    table = []
    i = 2**0
    n = 0
    while True:
        l = []
        dec2 = {}
        if i > 1:
            for j in range(len(T)):
                if not(table[n-1][j:j+i] in l):
                    l.append(table[n-1][j:j+i])
            l.sort()
            for j in range(len(l)):
                dec2[tuple(l[j])] = j
        row = []
        for j in range(len(T)):
            if i == 1:
                row.append(dec[T[j]])
            else:
                row.append(dec2[tuple(table[n-1][j:j+i])])
        table.append(row)
        flag = 0
        for j in range(len(row)):
            c = 0
            c = row.count(j)
            if c > 1:
                flag = 1
                break
        st.text(row)
        if flag == 0:
            return table
            break
        n += 1
        i = 2**n

def search_pattern(text, pattern):
    table1 = {}
    text.upper()
    pattern = pattern.upper()
    x = suffix(text)
    y = x[-1]
    for i in range(len(text)):
        table1[y[i]] = i
    suffix_array = dict(sorted(table1.items()))

    # Binary search for the pattern
    low, high = 0, len(text)
    while low <= high:
        mid = (low + high) // 2
        sufx = text[suffix_array[mid]:]
        if sufx.startswith(pattern):
            return suffix_array[mid]
        elif pattern < sufx:
            high = mid - 1
        else:
            low = mid + 1

    return -1  # Pattern not found

# assembly
from itertools import permutations

def overlap(a, b, min_length=3):
    start = 0
    while True:
        start = a.find(b[:min_length], start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1

def native_overlap(reads, k):
    olap = {}
    for a, b in permutations(reads, 2):
        olen = overlap(a, b, k)
        if olen > 0:
            olap[(b, a)] = olen
    return olap


st.title("Bio algorithms")
select = st.selectbox("**options**", ["naive_match", "complement", "CG_count", "reversed_seq", "Translation_Table", "Badchars", "indexing", "suffix", "native_overlap", "Read FASTA","kmp_search","approximate_match"])

# Add file uploader for FASTA files
if select == "Read FASTA":
    st.write("Upload a FASTA file:")
    uploaded_file = st.file_uploader("Choose a file", type=["fasta"])

    if uploaded_file is not None:
        # Read and display the sequences from the FASTA file
        sequences = read_fasta_file(uploaded_file)
        st.write("Sequences in the FASTA file:")
        for i, seq in enumerate(sequences):
            st.write(f"Sequence {i + 1}: {seq}")

elif select == "CG_count":
    s = st.text_input("**Enter sequence**")
    if len(s) > 0:
        st.write(f"**CG count in sequence : {CG_count(s)}**")
elif select == "reversed_seq":
    s = st.text_input("**Enter sequence**")
    if len(s) > 0:
        st.write(f"reversed sequence is : {reversed_seq(s)}")

elif select == "complement":
    s = st.text_input("**Enter sequence**")
    if len(s) > 0:
        st.write(f"complement of sequence is : {complement(s)}")
elif select == "Translation_Table":
    s = st.text_input("**Enter sequence**")
    if len(s) > 0:
        st.write(f"{Translation_Table(s)}")

elif select == "naive_match":
    s = st.text_input("**Enter sequence**")
    p = st.text_input("**Enter pattern**")
    x = naive_match(s, p)
    if x == -1:
        st.write(f'**Using naive match, the pattern is not exist in sequence**')
    else:
        st.write(f'**Using naive match, pattern exist in sequence at index  {x}**')

elif select == "Badchars":
    s = st.text_input("**Enter sequence**")
    p = st.text_input("**Enter pattern**")
    x = Badchars(s, p)
    if x == -1:
        st.write(f"**Using Bad character, the pattern is not exist in sequence**")
    else:
        st.write(f'**Using Bad character, pattern exist in sequence at index {x}** ')

elif select == "indexing":
    s = st.text_input("**Enter sequence**")
    p = st.text_input("**Enter pattern**")
    k = int(st.number_input('**length of substrings (k)**'))
    if len(s) > 1 and len(p) > 1:
        x = indexing(s, p, k)
        st.write(f'**Using Indexing pattern exist in sequence at index {x}**')

elif select == "suffix":
    s = st.text_input("**Enter sequence**")
    s1 = s + "$"
    p = st.text_input("***Enter pattern**")
    st.write(f'**suffix array of given sequence is**')
    x = search_pattern(s1, p)
    if x != -1:
        st.write(f'**using suffix array, pattern exist in sequence at index {x}**')
    else:
        st.write(f'**using suffix array, pattern is not exist in sequence**')
        

elif select == "native_overlap":
    s = st.text_input('**Enter reads separated by " , "**')
    k = int(st.number_input('**length of comparisons**'))
    reads = [read.strip() for read in s.split(",")]
    result = native_overlap(reads, k)
    result = {str(key): value for key, value in result.items()}
    st.write("**Overlaps**")
    st.write(result)
    
    
elif select == "kmp_search":
    s = st.text_input("**Enter sequence**")
    p = st.text_input("**Enter pattern**")
    x = kmp_search(s, p)
    if not x:
        st.write(f'**Using KMP algorithm, the pattern is not exist in sequence**')
    else:
        st.write(f'**Using KMP algorithm, pattern exist in sequence at indices  {x}**')
        
        
elif select == "approximate_match":
    s = st.text_input("**Enter sequence**")
    p = st.text_input("**Enter pattern**")
    max_distance = int(st.number_input("**Maximum allowed edit distance**"))
    x = approximate_match(s, p, max_distance)
    if not x:
        st.write(f'**Using Levenshtein distance, no approximate match found in the sequence**')
    else:
        st.write(f'**Using Levenshtein distance, approximate matches found at indices {x}**')


# In[ ]:




