"""
DAF-seq Primer Designer
Configuration and Input Parser Module
"""

import yaml
from dataclasses import dataclass, field
from typing import Optional, Tuple, List
from pathlib import Path


@dataclass
class PrimerConstraints:
    """Constraints for primer design"""
    min_length: int = 18
    max_length: int = 25
    max_guanines: int = 2
    no_gc_3prime_bases: int = 3
    tm_min: float = 58.0
    tm_max: float = 62.0
    min_gc_content: float = 20.0  # percentage
    max_gc_content: float = 60.0  # percentage
    max_poly_base: int = 4  # max consecutive same bases
    
    # Specificity checking
    check_specificity: bool = True  # Check for off-target binding
    allow_multiple_exact_matches: bool = False  # Allow primers with multiple binding sites
    specificity_min_match_length: int = 15  # Minimum length to consider as binding
    specificity_allow_mismatches: int = 2  # Allowed mismatches in off-target sites
    
    # Tm calculation parameters
    Na: float = 0.0  # mM
    K: float = 50.0    # mM
    Mg: float = 1.5   # mM
    dNTPs: float = 0.0  # mM
    primer_conc: float = 300.0  # nM


@dataclass
class PairingConstraints:
    """Constraints for primer pairing"""
    max_tm_difference: float = 5.0
    max_primer_reuse: int = 2
    min_primer_spacing: int = 25  # Minimum distance between primers in final set (for diversity)
    

@dataclass
class ScoringWeights:
    """Weights for scoring primers (higher = more penalty)"""
    guanine_penalty: float = 10.0  # penalty per G
    guanine_3prime_penalty: float = 50.0  # additional penalty for G near 3' end
    tm_deviation_penalty: float = 1.0  # penalty per degree from optimal
    gc_deviation_penalty: float = 0.5  # penalty per % from optimal GC
    length_deviation_penalty: float = 0.5  # penalty per bp from optimal length
    hairpin_penalty: float = 5.0
    self_dimer_penalty: float = 5.0
    poly_base_penalty: float = 3.0  # penalty per base over threshold
    
    # Position-based G penalties (multipliers for bases from 3' end)
    # Position 0 = 3' end, higher positions = closer to 5' end
    def get_g_position_penalty(self, position_from_3prime: int, primer_length: int) -> float:
        """
        Calculate penalty multiplier based on G position from 3' end.
        Exponential decay: positions closer to 3' end get higher penalties.
        """
        if position_from_3prime < 3:
            # Within the no-GC zone at 3' end
            return self.guanine_3prime_penalty
        elif position_from_3prime < 10:
            # Still relatively close to 3' end
            return self.guanine_penalty * (2.0 - position_from_3prime / 10.0)
        else:
            # Closer to 5' end, less problematic
            return self.guanine_penalty * 0.5


@dataclass
class TargetRegion:
    """Define the target amplicon region"""
    # Coordinates (0-based) - will be set either directly or via sequence matching
    core_start: Optional[int] = None  # 0-based
    core_end: Optional[int] = None    # 0-based, exclusive
    min_amplicon: int = 2000
    max_amplicon: int = 7000
    
    # Optional: define regions by sequence instead of coordinates
    core_sequence: Optional[str] = None
    forward_region_sequence: Optional[str] = None
    reverse_region_sequence: Optional[str] = None
    
    # Optional: specific regions where primers can be placed (coordinates)
    forward_region: Optional[Tuple[int, int]] = None  # (start, end)
    reverse_region: Optional[Tuple[int, int]] = None  # (start, end)
    
    def __post_init__(self):
        """Validate the target region"""
        # Check that either coordinates or sequence is provided for core region
        if self.core_start is None and self.core_sequence is None:
            raise ValueError("Must provide either core_start/core_end or core_sequence")
        
        if self.core_start is not None and self.core_end is not None:
            if self.core_start >= self.core_end:
                raise ValueError("core_start must be less than core_end")
            core_size = self.core_end - self.core_start
            if core_size > self.max_amplicon:
                raise ValueError(f"Core region ({core_size} bp) is larger than max_amplicon ({self.max_amplicon} bp)")
        
        if self.min_amplicon > self.max_amplicon:
            raise ValueError("min_amplicon must be less than or equal to max_amplicon")
        
        # Validate that we don't have both coordinates and sequences for the same region
        if self.forward_region is not None and self.forward_region_sequence is not None:
            raise ValueError("Cannot specify both forward_region coordinates and forward_region_sequence")
        if self.reverse_region is not None and self.reverse_region_sequence is not None:
            raise ValueError("Cannot specify both reverse_region coordinates and reverse_region_sequence")


@dataclass
class PrimerDesignConfig:
    """Main configuration for primer design"""
    sequence_file: str
    target_region: TargetRegion
    primer_constraints: PrimerConstraints = field(default_factory=PrimerConstraints)
    pairing_constraints: PairingConstraints = field(default_factory=PairingConstraints)
    scoring_weights: ScoringWeights = field(default_factory=ScoringWeights)
    num_primer_sets: int = 10
    output_file: Optional[str] = None  # Optional output filename
    output_dir: str = "results"  # Directory for output files
    verbose: bool = False  # Print detailed progress information
    
    @classmethod
    def from_yaml(cls, yaml_file: str) -> 'PrimerDesignConfig':
        """Load configuration from YAML file"""
        with open(yaml_file, 'r') as f:
            config_dict = yaml.safe_load(f)
        
        # Parse target region
        target_dict = config_dict['target_region']
        target_region = TargetRegion(
            core_start=target_dict.get('core_start'),
            core_end=target_dict.get('core_end'),
            core_sequence=target_dict.get('core_sequence'),
            min_amplicon=target_dict.get('min_amplicon', 2000),
            max_amplicon=target_dict.get('max_amplicon', 7000),
            forward_region=tuple(target_dict['forward_region']) if 'forward_region' in target_dict else None,
            reverse_region=tuple(target_dict['reverse_region']) if 'reverse_region' in target_dict else None,
            forward_region_sequence=target_dict.get('forward_region_sequence'),
            reverse_region_sequence=target_dict.get('reverse_region_sequence')
        )
        
        # Parse primer constraints
        primer_dict = config_dict.get('primer_constraints', {})
        primer_constraints = PrimerConstraints(**primer_dict)
        
        # Parse pairing constraints
        pairing_dict = config_dict.get('pairing_constraints', {})
        pairing_constraints = PairingConstraints(**pairing_dict)
        
        # Parse scoring weights
        scoring_dict = config_dict.get('scoring_weights', {})
        scoring_weights = ScoringWeights(**scoring_dict)
        
        # Get output settings
        output_file = config_dict.get('output_file')
        output_dir = config_dict.get('output_dir', 'results')
        
        return cls(
            sequence_file=config_dict['sequence_file'],
            target_region=target_region,
            primer_constraints=primer_constraints,
            pairing_constraints=pairing_constraints,
            scoring_weights=scoring_weights,
            num_primer_sets=config_dict.get('num_primer_sets', 10),
            output_file=output_file,
            output_dir=output_dir
        )
    
    def to_yaml(self, yaml_file: str):
        """Save configuration to YAML file (useful for generating templates)"""
        config_dict = {
            'sequence_file': self.sequence_file,
            'target_region': {
                'core_start': self.target_region.core_start,
                'core_end': self.target_region.core_end,
                'min_amplicon': self.target_region.min_amplicon,
                'max_amplicon': self.target_region.max_amplicon,
            },
            'primer_constraints': {
                'min_length': self.primer_constraints.min_length,
                'max_length': self.primer_constraints.max_length,
                'max_guanines': self.primer_constraints.max_guanines,
                'no_gc_3prime_bases': self.primer_constraints.no_gc_3prime_bases,
                'tm_min': self.primer_constraints.tm_min,
                'tm_max': self.primer_constraints.tm_max,
                'min_gc_content': self.primer_constraints.min_gc_content,
                'max_gc_content': self.primer_constraints.max_gc_content,
                'max_poly_base': self.primer_constraints.max_poly_base,
                'check_specificity': self.primer_constraints.check_specificity,
                'allow_multiple_exact_matches': self.primer_constraints.allow_multiple_exact_matches,
                'specificity_min_match_length': self.primer_constraints.specificity_min_match_length,
                'specificity_allow_mismatches': self.primer_constraints.specificity_allow_mismatches,
                'Na': self.primer_constraints.Na,
                'K': self.primer_constraints.K,
                'Mg': self.primer_constraints.Mg,
                'dNTPs': self.primer_constraints.dNTPs,
                'primer_conc': self.primer_constraints.primer_conc,
            },
            'pairing_constraints': {
                'max_tm_difference': self.pairing_constraints.max_tm_difference,
                'max_primer_reuse': self.pairing_constraints.max_primer_reuse,
                'min_primer_spacing': self.pairing_constraints.min_primer_spacing,
            },
            'scoring_weights': {
                'guanine_penalty': self.scoring_weights.guanine_penalty,
                'guanine_3prime_penalty': self.scoring_weights.guanine_3prime_penalty,
                'tm_deviation_penalty': self.scoring_weights.tm_deviation_penalty,
                'gc_deviation_penalty': self.scoring_weights.gc_deviation_penalty,
                'length_deviation_penalty': self.scoring_weights.length_deviation_penalty,
                'hairpin_penalty': self.scoring_weights.hairpin_penalty,
                'self_dimer_penalty': self.scoring_weights.self_dimer_penalty,
                'poly_base_penalty': self.scoring_weights.poly_base_penalty,
            },
            'num_primer_sets': self.num_primer_sets,
            'output_file': self.output_file,
            'output_dir': self.output_dir
        }
        
        if self.target_region.forward_region:
            config_dict['target_region']['forward_region'] = list(self.target_region.forward_region)
        if self.target_region.reverse_region:
            config_dict['target_region']['reverse_region'] = list(self.target_region.reverse_region)
        
        with open(yaml_file, 'w') as f:
            yaml.dump(config_dict, f, default_flow_style=False, sort_keys=False)


def create_example_config(output_file: str = "primer_design_config.yaml"):
    """Create an example configuration file"""
    example_config = PrimerDesignConfig(
        sequence_file="template.fasta",
        target_region=TargetRegion(
            core_start=1000,
            core_end=3000,
            min_amplicon=2000,
            max_amplicon=7000,
            forward_region=(500, 1100),
            reverse_region=(2900, 3500)
        ),
        output_file="my_primers.txt",
        output_dir="results"
    )
    example_config.to_yaml(output_file)
    print(f"Example configuration saved to {output_file}")


if __name__ == "__main__":
    # Generate example config file
    create_example_config()