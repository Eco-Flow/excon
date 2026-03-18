#!/usr/bin/env python3

"""
Rescale branch lengths in a phylogenetic tree.
This script reads a Newick format tree and multiplies all branch lengths
by a specified scaling factor to prevent numerical precision issues in 
downstream analyses like CAFE.

This version uses regex parsing instead of ete3 to avoid Python 3.13+ compatibility issues.
"""

import sys
import argparse
import re


def rescale_tree(input_file, output_file, scale_factor=1000):
    """
    Read a tree, rescale branch lengths, and write to output file.
    
    Args:
        input_file: Path to input Newick tree file
        output_file: Path to output Newick tree file
        scale_factor: Factor to multiply branch lengths by (default: 1000)
    """
    try:
        # Read the tree
        with open(input_file, 'r') as f:
            tree_str = f.read().strip()
        
        # Pattern to match branch lengths in Newick format
        # Matches :0.123 or :1.23e-5 etc.
        branch_length_pattern = r':(\d+\.?\d*(?:[eE][+-]?\d+)?)'
        
        rescaled_count = 0
        
        def rescale_match(match):
            """Helper function to rescale a matched branch length"""
            nonlocal rescaled_count
            original_length = float(match.group(1))
            new_length = original_length * scale_factor
            rescaled_count += 1
            return f':{new_length}'
        
        # Replace all branch lengths with rescaled values
        rescaled_tree = re.sub(branch_length_pattern, rescale_match, tree_str)
        
        # Write rescaled tree
        with open(output_file, 'w') as f:
            f.write(rescaled_tree)
            if not rescaled_tree.endswith('\n'):
                f.write('\n')
        
        print(f"✅ Tree rescaling completed successfully!")
        print(f"   Rescaled {rescaled_count} branches by factor of {scale_factor}")
        print(f"   Output written to: {output_file}")
        
        return True
        
    except FileNotFoundError:
        print(f"❌ Error: Input file '{input_file}' not found", file=sys.stderr)
        return False
    except Exception as e:
        print(f"❌ Error rescaling tree: {str(e)}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return False


def main():
    parser = argparse.ArgumentParser(
        description='Rescale branch lengths in a phylogenetic tree',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Rescale by 1000x (default)
  %(prog)s -i input_tree.nwk -o output_tree.nwk
  
  # Rescale by custom factor
  %(prog)s -i input_tree.nwk -o output_tree.nwk -s 500
        """
    )
    
    parser.add_argument('-i', '--input', required=True,
                        help='Input Newick tree file')
    parser.add_argument('-o', '--output', required=True,
                        help='Output Newick tree file')
    parser.add_argument('-s', '--scale', type=float, default=1000,
                        help='Scaling factor for branch lengths (default: 1000)')
    
    args = parser.parse_args()
    
    # Validate scale factor
    if args.scale <= 0:
        print("❌ Error: Scale factor must be positive", file=sys.stderr)
        sys.exit(1)
    
    # Rescale the tree
    success = rescale_tree(args.input, args.output, args.scale)
    
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
