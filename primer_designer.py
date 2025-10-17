"""
DAF-seq Primer Designer - Main Script
Designs PCR primers for DAF-seq treated DNA
"""

import sys
import argparse
from pathlib import Path
from datetime import datetime
import os

from config import PrimerDesignConfig
from sequence_parser import SequenceParser
from primer_generator import PrimerGenerator
from primer_evaluator import PrimerEvaluator
from primer_pairer import PrimerPairer
from primer_specificity import SpecificityChecker


def print_header():
    """Print program header"""
    print("=" * 80)
    print("DAF-seq treated DNA Primer Designer")
    print("=" * 80)
    print()


def print_primers_summary(pairs, num_to_show=10):
    """
    Print a summary of the best primer pairs.
    Positions shown are 1-based for compatibility with standard software.
    
    Args:
        pairs: List of PrimerPair objects
        num_to_show: Number of pairs to display
    """
    print("\n" + "=" * 80)
    print(f"TOP {min(num_to_show, len(pairs))} PRIMER PAIRS")
    print("=" * 80)
    print("(Note: Positions are 1-based)")
    
    for i, pair in enumerate(pairs[:num_to_show], 1):
        print(f"\n--- Primer Set #{i} (Score: {pair.pair_score:.2f}) ---")
        print(f"Amplicon Size: {pair.amplicon_size} bp")
        print(f"Tm Difference: {pair.tm_difference:.2f}°C")
        print(f"Covers Core Region: {'Yes' if pair.covers_core else 'No'}")
        print(f"Primer Dimer Risk: {'Yes' if pair.has_primer_dimer else 'No'}")
        
        print(f"\n  Forward Primer:")
        print(f"    Sequence (5'->3'): {pair.forward.sequence}")
        print(f"    Position: {pair.forward.position + 1}")  # Convert to 1-based
        print(f"    Length: {pair.forward.length} bp")
        print(f"    Tm: {pair.forward.tm:.2f}°C")
        print(f"    GC%: {pair.forward.gc_content:.1f}%")
        print(f"    Number of Gs: {pair.forward.num_guanines}")
        if pair.forward.num_guanines > 0:
            print(f"    G positions from 3' end: {pair.forward.guanine_positions}")
        print(f"    Score: {pair.forward.score:.2f}")
        
        print(f"\n  Reverse Primer:")
        print(f"    Sequence (5'->3'): {pair.reverse.sequence}")
        print(f"    Position: {pair.reverse.position + 1}")  # Convert to 1-based
        print(f"    Length: {pair.reverse.length} bp")
        print(f"    Tm: {pair.reverse.tm:.2f}°C")
        print(f"    GC%: {pair.reverse.gc_content:.1f}%")
        print(f"    Number of Gs: {pair.reverse.num_guanines}")
        if pair.reverse.num_guanines > 0:
            print(f"    G positions from 3' end: {pair.reverse.guanine_positions}")
        print(f"    Score: {pair.reverse.score:.2f}")


def save_results_to_file(pairs, output_file, config, verbose=True):
    """
    Save primer pairs to a tab-delimited file.
    Positions are converted to 1-based for compatibility with standard software.
    
    Args:
        pairs: List of PrimerPair objects
        output_file: Output file path
        config: PrimerDesignConfig object
        verbose: Whether to print save confirmation
    """
    with open(output_file, 'w') as f:
        # Write header with metadata
        f.write(f"# DAF-seq Primer Design Results\n")
        f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# Sequence File: {config.sequence_file}\n")
        f.write(f"# Core Region: {config.target_region.core_start + 1}-{config.target_region.core_end} (1-based)\n")
        f.write(f"# Amplicon Range: {config.target_region.min_amplicon}-{config.target_region.max_amplicon} bp\n")
        f.write(f"# NOTE: All positions are 1-based\n")
        f.write(f"#\n")
        
        # Write column headers
        headers = [
            "Rank",
            "Pair_Score",
            "Amplicon_Size",
            "Tm_Difference",
            "Forward_Sequence",
            "Forward_Position",
            "Forward_Length",
            "Forward_Tm",
            "Forward_GC%",
            "Forward_Num_Gs",
            "Forward_G_Positions",
            "Forward_Score",
            "Reverse_Sequence",
            "Reverse_Position",
            "Reverse_Length",
            "Reverse_Tm",
            "Reverse_GC%",
            "Reverse_Num_Gs",
            "Reverse_G_Positions",
            "Reverse_Score",
            "Primer_Dimer",
            "Covers_Core"
        ]
        f.write("\t".join(headers) + "\n")
        
        # Write data
        for i, pair in enumerate(pairs, 1):
            row = [
                str(i),
                f"{pair.pair_score:.2f}",
                str(pair.amplicon_size),
                f"{pair.tm_difference:.2f}",
                pair.forward.sequence,
                str(pair.forward.position + 1),  # Convert to 1-based
                str(pair.forward.length),
                f"{pair.forward.tm:.2f}",
                f"{pair.forward.gc_content:.1f}",
                str(pair.forward.num_guanines),
                ",".join(map(str, pair.forward.guanine_positions)) if pair.forward.guanine_positions else "None",
                f"{pair.forward.score:.2f}",
                pair.reverse.sequence,
                str(pair.reverse.position + 1),  # Convert to 1-based
                str(pair.reverse.length),
                f"{pair.reverse.tm:.2f}",
                f"{pair.reverse.gc_content:.1f}",
                str(pair.reverse.num_guanines),
                ",".join(map(str, pair.reverse.guanine_positions)) if pair.reverse.guanine_positions else "None",
                f"{pair.reverse.score:.2f}",
                "Yes" if pair.has_primer_dimer else "No",
                "Yes" if pair.covers_core else "No"
            ]
            f.write("\t".join(row) + "\n")
    
    if verbose:
        print(f"\nResults saved to: {output_file}")


def run_primer_design(config_file, output_file=None, verbose_override=None):
    """
    Main function to run the primer design pipeline.
    
    Args:
        config_file: Path to YAML configuration file
        output_file: Optional output file for results
        verbose_override: Override verbose setting from config (True/False/None)
    """
    # Load configuration
    config = PrimerDesignConfig.from_yaml(config_file)
    
    # Override verbose if specified on command line
    if verbose_override is not None:
        config.verbose = verbose_override
    
    if config.verbose:
        print_header()
        print(f"Loading configuration from: {config_file}")
        print(f"✓ Configuration loaded")
        print(f"  Output file: {config.output_file}")
        print(f"  Output dir: {config.output_dir}\n")
    
    # Parse sequence
    if config.verbose:
        print("Step 1: Parsing sequence file")
        print("-" * 80)
    
    parser = SequenceParser(config.sequence_file, verbose=config.verbose)
    
    # Resolve any sequence-based regions to coordinates
    config.target_region = parser.resolve_target_region(config.target_region)
    
    # Get primer search regions
    fwd_region, rev_region = parser.get_primer_search_regions(config.target_region)
    
    if config.verbose:
        print(f"Forward primer search region: {fwd_region[0]}-{fwd_region[1]} ({fwd_region[1]-fwd_region[0]} bp)")
        print(f"Reverse primer search region: {rev_region[0]}-{rev_region[1]} ({rev_region[1]-rev_region[0]} bp)")
        print(f"✓ Sequence parsing complete\n")
    
    # Generate primers
    if config.verbose:
        print("Step 2: Generating candidate primers")
        print("-" * 80)
    
    generator = PrimerGenerator(config.primer_constraints)
    
    fwd_search_seq = parser.get_forward_search_sequence(fwd_region)
    rev_search_seq = parser.get_reverse_search_sequence(rev_region)
    
    fwd_primers = generator.generate_primers(fwd_search_seq, fwd_region[0], is_forward=True)
    rev_primers = generator.generate_primers(rev_search_seq, rev_region[0], is_forward=False, region_end=rev_region[1])
    
    if config.verbose:
        print(f"Generated {len(fwd_primers)} forward primer candidates")
        print(f"Generated {len(rev_primers)} reverse primer candidates")
        
        fwd_stats = generator.get_statistics(fwd_primers)
        rev_stats = generator.get_statistics(rev_primers)
        
        print(f"\nForward primer statistics:")
        print(f"  Average length: {fwd_stats['avg_length']:.1f} bp")
        print(f"  Average GC%: {fwd_stats['avg_gc_content']:.1f}%")
        print(f"  Primers with 0 Gs: {fwd_stats['primers_with_0_g']}")
        print(f"  Primers with 1 G: {fwd_stats['primers_with_1_g']}")
        print(f"  Primers with 2 Gs: {fwd_stats['primers_with_2_g']}")
        
        print(f"\nReverse primer statistics:")
        print(f"  Average length: {rev_stats['avg_length']:.1f} bp")
        print(f"  Average GC%: {rev_stats['avg_gc_content']:.1f}%")
        print(f"  Primers with 0 Gs: {rev_stats['primers_with_0_g']}")
        print(f"  Primers with 1 G: {rev_stats['primers_with_1_g']}")
        print(f"  Primers with 2 Gs: {rev_stats['primers_with_2_g']}")
    
    if len(fwd_primers) == 0 or len(rev_primers) == 0:
        print("\n⚠ WARNING: No primers found with current constraints!")
        print("Consider relaxing constraints (e.g., allow more Gs, wider Tm range)")
        return
    
    if config.verbose:
        print(f"✓ Primer generation complete\n")
    
    # Evaluate primers
    if config.verbose:
        print("Step 3: Evaluating primers (calculating Tm, checking structures, scoring)")
        print("-" * 80)
    
    evaluator = PrimerEvaluator(config.primer_constraints, config.scoring_weights)
    
    fwd_primers = evaluator.evaluate_primers(fwd_primers, verbose=config.verbose)
    rev_primers = evaluator.evaluate_primers(rev_primers, verbose=config.verbose)
    
    # Filter by Tm
    fwd_primers = evaluator.filter_by_tm(fwd_primers)
    rev_primers = evaluator.filter_by_tm(rev_primers)
    
    if config.verbose:
        print(f"\nAfter Tm filtering:")
        print(f"  Forward primers: {len(fwd_primers)}")
        print(f"  Reverse primers: {len(rev_primers)}")
    
    if len(fwd_primers) == 0 or len(rev_primers) == 0:
        print("\n⚠ WARNING: No primers within Tm range!")
        print(f"Tm range: {config.primer_constraints.tm_min}-{config.primer_constraints.tm_max}°C")
        print("Consider widening the Tm range")
        return
    
    if config.verbose:
        print(f"✓ Primer evaluation complete\n")
    
    # Check specificity
    if config.primer_constraints.check_specificity:
        if config.verbose:
            print("Step 3.5: Checking primer specificity")
            print("-" * 80)
        
        specificity_checker = SpecificityChecker(
            parser.sequence,
            min_match_length=config.primer_constraints.specificity_min_match_length,
            allow_mismatches=config.primer_constraints.specificity_allow_mismatches
        )
        
        if config.verbose:
            print("Checking forward primers...")
        fwd_primers_before = len(fwd_primers)
        fwd_primers = specificity_checker.filter_specific_primers(fwd_primers, verbose=config.verbose, debug=False)
        
        if config.verbose:
            print(f"\nChecking reverse primers...")
        rev_primers_before = len(rev_primers)
        rev_primers = specificity_checker.filter_specific_primers(rev_primers, verbose=config.verbose, debug=False)
        
        if len(fwd_primers) == 0 or len(rev_primers) == 0:
            print("\n⚠ WARNING: No specific primers found!")
            print("Consider:")
            print("  - Disabling specificity checking (set check_specificity: false)")
            print("  - Allowing multiple exact matches (set allow_multiple_exact_matches: true)")
            print("  - Increasing specificity_allow_mismatches")
            return
        
        if config.verbose:
            print(f"\n✓ Specificity checking complete")
            print(f"  Forward: {fwd_primers_before} → {len(fwd_primers)} primers")
            print(f"  Reverse: {rev_primers_before} → {len(rev_primers)} primers\n")
    
    # Get best individual primers for pairing
    fwd_primers = evaluator.get_best_primers(fwd_primers, n=min(200, len(fwd_primers)))
    rev_primers = evaluator.get_best_primers(rev_primers, n=min(200, len(rev_primers)))
    
    if config.verbose:
        print(f"Using top {len(fwd_primers)} forward and top {len(rev_primers)} reverse primers for pairing")
    
    # Pair primers
    if config.verbose:
        print("\nStep 4: Pairing primers")
        print("-" * 80)
    
    pairer = PrimerPairer(config.pairing_constraints, config.scoring_weights, config.target_region)
    
    pairs = pairer.create_pairs(fwd_primers, rev_primers, verbose=config.verbose)
    
    if len(pairs) == 0:
        print("\n⚠ WARNING: No valid primer pairs found!")
        print("Possible issues:")
        print("  - Amplicon size constraints too restrictive")
        print("  - Tm difference threshold too strict")
        print("  - Primers don't cover core region")
        return
    
    # Get best pairs
    best_pairs = pairer.get_best_pairs(
        pairs, 
        n=config.num_primer_sets, 
        apply_reuse_filter=True,
        enforce_spacing=True
    )
    
    if len(best_pairs) < config.num_primer_sets and config.verbose:
        print(f"\nNote: Only {len(best_pairs)} diverse primer pairs found with minimum spacing of {config.pairing_constraints.min_primer_spacing} bp")
        print(f"Consider reducing 'min_primer_spacing' if you need more primer sets")
    
    if config.verbose:
        print(f"✓ Primer pairing complete\n")
    
    # Print results
    if config.verbose:
        print_primers_summary(best_pairs, num_to_show=config.num_primer_sets)
    
    # Create output directory if it doesn't exist
    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Determine output file path
    if output_file:
        # Command-line argument takes precedence
        output_path = Path(output_file)
    elif config.output_file:
        # Use config file specification
        output_path = output_dir / config.output_file
    else:
        # Auto-generate filename with timestamp
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        output_path = output_dir / f"primer_results_{timestamp}.txt"
    
    # Save results
    save_results_to_file(best_pairs, str(output_path), config, verbose=config.verbose)
    
    if config.verbose:
        print("\n" + "=" * 80)
        print("Primer design complete!")
        print("=" * 80)
    else:
        # Always print the output file location even in non-verbose mode
        print(f"Primer design complete. Results saved to: {output_path}")


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description='Design PCR primers for DAF-seq treated DNA',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with config file
  python primer_designer.py config.yaml
  
  # Run with config file and specify output
  python primer_designer.py config.yaml -o my_primers.txt
  
  # Generate example config file
  python primer_designer.py --create-example-config
        """
    )
    
    parser.add_argument('config', nargs='?', help='YAML configuration file')
    parser.add_argument('-o', '--output', help='Output file for results')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Print detailed progress information')
    parser.add_argument('-q', '--quiet', action='store_true',
                       help='Run silently (overrides config verbose setting)')
    parser.add_argument('--create-example-config', action='store_true',
                       help='Create an example configuration file and exit')
    
    args = parser.parse_args()
    
    # Handle example config creation
    if args.create_example_config:
        from config import create_example_config
        create_example_config("primer_design_config.yaml")
        return
    
    # Validate arguments
    if not args.config:
        parser.print_help()
        sys.exit(1)
    
    # Check if config file exists
    if not Path(args.config).exists():
        print(f"Error: Configuration file '{args.config}' not found")
        sys.exit(1)
    
    # Determine verbose setting
    verbose_override = None
    if args.verbose:
        verbose_override = True
    elif args.quiet:
        verbose_override = False
    
    # Run primer design
    try:
        run_primer_design(args.config, args.output, verbose_override)
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()