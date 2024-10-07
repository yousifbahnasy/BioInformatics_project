<h1 align="center">Hi 👋, I'm yousif</h1>
<h3 align="center">Future Data Scientist | AI Engineer | ML Engineer
</h3>


<h3 align="left">Connect with me:</h3>
<p align="left">
<a href="https://linkedin.com/in/yousif-bahnasy" target="blank"><img align="center" src="https://raw.githubusercontent.com/rahuldkjain/github-profile-readme-generator/master/src/images/icons/Social/linked-in-alt.svg" alt="yousifbahnasy" height="30" width="40" /></a>
<a href="https://instagram.com/yousif_bahnasy" target="blank"><img align="center" src="https://raw.githubusercontent.com/rahuldkjain/github-profile-readme-generator/master/src/images/icons/Social/instagram.svg" alt="yousif_bahnasy" height="30" width="40" /></a>
</p>

# Bioinformatics Algorithms Tool

This project provides a set of bioinformatics algorithms implemented using Python and Streamlit. These algorithms help analyze DNA/RNA sequences for a variety of tasks such as pattern matching, sequence translation, and computing overlaps between sequences.

## Features
- **FASTA File Reader**: Upload and analyze DNA/RNA sequences from FASTA files.
- **Pattern Matching Algorithms**: Includes Naive Match, KMP, and Approximate Matching using Levenshtein Distance.
- **CG Count Calculation**: Calculate the CG content of a sequence.
- **Sequence Manipulations**: Reverse and complement DNA sequences.
- **Amino Acid Translation**: Translate DNA sequences into their corresponding amino acid chains.
- **Advanced Search Algorithms**: Utilize Suffix Arrays, Bad Character Heuristics, and Indexed Search.
- **Overlap Graph Construction**: Find overlaps between sequences for sequence assembly.

## Usage
To run the project locally, simply run the Streamlit app:

```bash
streamlit run app.py
```
This will open the tool in your browser where you can interact with the different algorithms by selecting from the dropdown menu.

### Supported Operations

- **Read FASTA File**: Upload a `.fasta` file to view sequences.
- **CG Content**: Calculate the CG percentage of a given DNA sequence.
- **Sequence Translation**: Convert DNA sequences to amino acids.
- **Pattern Matching**: Search for a pattern in a sequence using Naive, KMP, or approximate matching.
- **Reversal and Complement**: Get the reverse or complement of a sequence.
- **Search with Suffix Array**: Efficiently search for patterns using suffix arrays.

## Algorithms
This tool includes the following algorithms:

- **Naive Match**: A simple brute-force pattern matching algorithm.
- **KMP Search**: An efficient string-matching algorithm using a partial match table.
- **Levenshtein Distance**: For approximate string matching within a given edit distance.
- **Bad Character Heuristic**: A part of the Boyer-Moore string search algorithm.
- **Indexed Search**: Search using sorted substrings and bisecting for fast lookups.
- **Suffix Array**: Construct a suffix array for efficient pattern searching in large sequences.
- **Overlap Graph**: Identify overlapping regions between sequences.

## Example

Here’s a quick example of how you can use the tool to count the CG content of a DNA sequence:

1. Run the app and select **CG count**.
2. Enter the DNA sequence: `ATCGATCGATCG`.
3. The tool will output the CG percentage.



