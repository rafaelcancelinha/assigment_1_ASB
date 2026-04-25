#!/usr/bin/env python3
"""
SeqConcat: Sequence Concatenation Tool
Concatenates sequence data from multiple fragments in FASTA, Phylip, or Nexus formats.
"""

import argparse
import sys
import os
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple, Set


class SequenceParser:
    """Base class for parsing sequence files"""
    
    @staticmethod
    def parse_fasta(filepath: str) -> Dict[str, str]:
        """Parse a FASTA file and return {sequence_name: sequence}"""
        sequences = {}
        current_name = None
        current_seq = []
        
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.rstrip('\n')
                    if not line:
                        continue
                    if line.startswith('>'):
                        if current_name:
                            sequences[current_name] = ''.join(current_seq)
                        current_name = line[1:].split()[0]
                        current_seq = []
                    else:
                        current_seq.append(line)
                
                if current_name:
                    sequences[current_name] = ''.join(current_seq)
        except IOError as e:
            print(f"Error reading FASTA file {filepath}: {e}", file=sys.stderr)
            sys.exit(1)
        
        return sequences
    
    @staticmethod
    def parse_phylip(filepath: str) -> Dict[str, str]:
        """Parse a Phylip format file and return {sequence_name: sequence}"""
        sequences = {}
        
        try:
            with open(filepath, 'r') as f:
                # Read header line
                header = f.readline().strip().split()
                if len(header) < 2:
                    print(f"Invalid Phylip file format: {filepath}", file=sys.stderr)
                    sys.exit(1)
                
                ntax, nchar = int(header[0]), int(header[1])
                
                # Read sequences
                for line_num, line in enumerate(f, 1):
                    if line.strip():
                        parts = line.split(None, 1)
                        if len(parts) == 2:
                            name, seq = parts
                            sequences[name] = seq.replace(' ', '')
        except (IOError, ValueError) as e:
            print(f"Error reading Phylip file {filepath}: {e}", file=sys.stderr)
            sys.exit(1)
        
        return sequences
    
    @staticmethod
    def parse_nexus(filepath: str) -> Dict[str, str]:
        """Parse a Nexus format file and return {sequence_name: sequence}"""
        sequences = {}
        in_matrix = False
        
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('['):
                        continue
                    
                    if line.lower().startswith('matrix'):
                        in_matrix = True
                        continue
                    
                    if line == ';':
                        if in_matrix:
                            break
                        continue
                    
                    if in_matrix:
                        parts = line.split()
                        if len(parts) >= 2:
                            name = parts[0]
                            seq = ''.join(parts[1:]).replace(' ', '')
                            sequences[name] = seq
        except IOError as e:
            print(f"Error reading Nexus file {filepath}: {e}", file=sys.stderr)
            sys.exit(1)
        
        return sequences


class SequenceWriter:
    """Class for writing sequence files"""
    
    @staticmethod
    def write_fasta(sequences: Dict[str, str], filepath: str, line_width: int = 70) -> None:
        """Write sequences to FASTA file"""
        try:
            with open(filepath, 'w') as f:
                for name in sorted(sequences.keys()):
                    seq = sequences[name]
                    f.write(f">{name}\n")
                    for i in range(0, len(seq), line_width):
                        f.write(seq[i:i+line_width] + '\n')
        except IOError as e:
            print(f"Error writing FASTA file {filepath}: {e}", file=sys.stderr)
            sys.exit(1)
    
    @staticmethod
    def write_phylip(sequences: Dict[str, str], filepath: str) -> None:
        """Write sequences to Phylip format file"""
        try:
            if not sequences:
                print("No sequences to write", file=sys.stderr)
                return
            
            ntax = len(sequences)
            nchar = len(next(iter(sequences.values())))
            
            with open(filepath, 'w') as f:
                f.write(f"{ntax} {nchar}\n")
                for name in sorted(sequences.keys()):
                    seq = sequences[name]
                    f.write(f"{name:<10} {seq}\n")
        except IOError as e:
            print(f"Error writing Phylip file {filepath}: {e}", file=sys.stderr)
            sys.exit(1)
    
    @staticmethod
    def write_nexus(sequences: Dict[str, str], filepath: str, line_width: int = 70) -> None:
        """Write sequences to Nexus format file"""
        try:
            if not sequences:
                print("No sequences to write", file=sys.stderr)
                return
            
            ntax = len(sequences)
            nchar = len(next(iter(sequences.values())))
            
            with open(filepath, 'w') as f:
                f.write("#NEXUS\n")
                f.write(f"BEGIN DATA;\n")
                f.write(f"DIMENSIONS NTAX={ntax} NCHAR={nchar};\n")
                f.write(f"FORMAT DATATYPE=DNA GAP=- MISSING=N;\n")
                f.write(f"MATRIX\n")
                
                for name in sorted(sequences.keys()):
                    seq = sequences[name]
                    f.write(f"{name:<20} ")
                    for i in range(0, len(seq), line_width):
                        if i > 0:
                            f.write(" " * 20)
                        f.write(seq[i:i+line_width] + '\n')
                
                f.write(f";\nEND;\n")
        except IOError as e:
            print(f"Error writing Nexus file {filepath}: {e}", file=sys.stderr)
            sys.exit(1)


class SequenceConcatenator:
    """Main class for concatenating sequences"""
    
    def __init__(self, output_format: str = 'fasta'):
        self.sequences = defaultdict(dict)  # {species: {fragment: sequence}}
        self.fragments = []  # Order of fragments
        self.fragment_lengths = {}  # {fragment: length}
        self.output_format = output_format.lower()
    
    def detect_format(self, filepath: str) -> str:
        """Detect file format based on extension and content"""
        filepath_lower = filepath.lower()
        
        if filepath_lower.endswith('.fasta') or filepath_lower.endswith('.fa'):
            return 'fasta'
        elif filepath_lower.endswith('.phy') or filepath_lower.endswith('.phylip'):
            return 'phylip'
        elif filepath_lower.endswith('.nex') or filepath_lower.endswith('.nexus'):
            return 'nexus'
        
        # Guess by content
        try:
            with open(filepath, 'r') as f:
                first_line = f.readline().strip()
                if first_line.startswith('>'):
                    return 'fasta'
                elif first_line.lower().startswith('#nexus'):
                    return 'nexus'
                else:
                    return 'phylip'
        except:
            return 'fasta'
    
    def load_sequences(self, filepath: str, fragment_name: str = None) -> None:
        """Load sequences from a file"""
        file_format = self.detect_format(filepath)
        
        if file_format == 'fasta':
            sequences = SequenceParser.parse_fasta(filepath)
        elif file_format == 'phylip':
            sequences = SequenceParser.parse_phylip(filepath)
        elif file_format == 'nexus':
            sequences = SequenceParser.parse_nexus(filepath)
        else:
            print(f"Unknown format for file {filepath}", file=sys.stderr)
            sys.exit(1)
        
        # Use filename as fragment name if not provided
        if fragment_name is None:
            fragment_name = Path(filepath).stem
        
        # Store sequences
        if not sequences:
            print(f"Warning: No sequences found in {filepath}", file=sys.stderr)
            return
        
        # Get sequence length (assuming all sequences in file have same length)
        seq_length = len(next(iter(sequences.values())))
        self.fragment_lengths[fragment_name] = seq_length
        
        for species, seq in sequences.items():
            self.sequences[species][fragment_name] = seq
        
        if fragment_name not in self.fragments:
            self.fragments.append(fragment_name)
        
        print(f"✓ Loaded {fragment_name}: {len(sequences)} sequences ({seq_length} bp)")
    
    def concatenate(self, remove_incomplete: bool = False, 
                   add_missing: bool = True) -> Dict[str, str]:
        """Concatenate all sequences"""
        concatenated = {}
        all_species = set(self.sequences.keys())
        
        print(f"\n{'='*60}")
        print(f"CONCATENATION RESULTS")
        print(f"{'='*60}")
        print(f"Total species: {len(all_species)}")
        print(f"Total fragments: {len(self.fragments)}")
        
        # Get complete species
        if remove_incomplete:
            complete_species = {
                sp for sp in all_species 
                if all(frag in self.sequences[sp] for frag in self.fragments)
            }
            incomplete_species = all_species - complete_species
            
            print(f"Complete species: {len(complete_species)}/{len(all_species)}")
            if incomplete_species and len(incomplete_species) <= 10:
                print(f"Removed species: {', '.join(sorted(incomplete_species))}")
            elif incomplete_species:
                print(f"Removed {len(incomplete_species)} incomplete species")
            
            species_to_process = complete_species
        else:
            species_to_process = all_species
        
        # Concatenate sequences
        for species in species_to_process:
            concat_seq = ""
            for fragment in self.fragments:
                if fragment in self.sequences[species]:
                    concat_seq += self.sequences[species][fragment]
                elif add_missing:
                    # Add N for missing sequences
                    gap_length = self.fragment_lengths[fragment]
                    concat_seq += 'N' * gap_length
                else:
                    # This shouldn't happen if remove_incomplete is True
                    gap_length = self.fragment_lengths[fragment]
                    concat_seq += '-' * gap_length
            
            concatenated[species] = concat_seq
        
        # Print fragment info
        print(f"\nFragment information:")
        total_length = 0
        for fragment in self.fragments:
            length = self.fragment_lengths[fragment]
            count = sum(1 for sp in self.sequences.values() if fragment in sp)
            total_length += length
            print(f"  {fragment:<20} {length:6d} bp | {count:3d} species")
        
        print(f"{'─'*60}")
        print(f"  {'TOTAL':<20} {total_length:6d} bp")
        
        return concatenated
    
    def write_output(self, concatenated: Dict[str, str], output_file: str) -> None:
        """Write concatenated sequences to file"""
        if self.output_format == 'fasta':
            SequenceWriter.write_fasta(concatenated, output_file)
        elif self.output_format == 'phylip':
            SequenceWriter.write_phylip(concatenated, output_file)
        elif self.output_format == 'nexus':
            SequenceWriter.write_nexus(concatenated, output_file)
        
        print(f"\n✓ Output written to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='SeqConcat: Concatenate sequence data from multiple fragments',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Simple: just pass the files
  python seqconcat.py file1.fasta file2.fasta file3.fasta
  
  # With output name
  python seqconcat.py *.fasta -o my_output.fasta
  
  # Remove incomplete species
  python seqconcat.py gene1.fasta gene2.fasta gene3.fasta --remove-incomplete
  
  # Output as Nexus format
  python seqconcat.py *.fasta -o output.nexus --output-format nexus
  
  # With custom fragment names
  python seqconcat.py f1.fasta f2.fasta f3.fasta -n ITS LSU RPB1
        """
    )
    
    parser.add_argument('input', nargs='+',
                       help='Input sequence files (FASTA, Phylip, or Nexus format)')
    parser.add_argument('-o', '--output',
                       help='Output file path (default: concatenated.fasta)')
    parser.add_argument('-n', '--names', nargs='+',
                       help='Fragment names (optional, must match number of input files)')
    parser.add_argument('--output-format', choices=['fasta', 'phylip', 'nexus'],
                       default='fasta',
                       help='Output file format (default: fasta)')
    parser.add_argument('--remove-incomplete', action='store_true',
                       help='Remove species not represented in all fragments')
    parser.add_argument('--add-missing', action='store_true', default=True,
                       help='Add 100%% N sequences for missing data (default: True)')
    parser.add_argument('--no-missing', dest='add_missing', action='store_false',
                       help='Use gaps (-) instead of N for missing data')
    
    args = parser.parse_args()
    
    if not args.output:
        args.output = 'concatenated.fasta'

    if args.names and len(args.names) != len(args.input):
        print("Error: Number of fragment names must match number of input files", file=sys.stderr)
        sys.exit(1)
    
    concatenator = SequenceConcatenator(output_format=args.output_format)
    
    print("SeqConcat v1.0 - Sequence Concatenation Tool")
    print("="*60)
    print("\nLoading sequences...")
    
    # Load input files
    for i, input_file in enumerate(args.input):
        if not os.path.exists(input_file):
            print(f"Error: File not found: {input_file}", file=sys.stderr)
            sys.exit(1)
        
        fragment_name = args.names[i] if args.names else None
        concatenator.load_sequences(input_file, fragment_name)
    
    # Concatenate sequences
    concatenated = concatenator.concatenate(
        remove_incomplete=args.remove_incomplete,
        add_missing=args.add_missing
    )
    
    # Write output
    concatenator.write_output(concatenated, args.output)
    print("\n✓ Concatenation completed successfully!")


if __name__ == '__main__':
    main()
