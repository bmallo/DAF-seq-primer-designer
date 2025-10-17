"""
Primer Generator Module
Generates candidate primers using sliding window approach and filters by basic constraints
"""

from typing import List, Dict, Tuple
from dataclasses import dataclass
from config import PrimerConstraints


@dataclass
class PrimerCandidate:
    """Represents a candidate primer with its properties"""
    sequence: str
    position: int  # Position in the original template (0-based)
    length: int
    is_forward: bool  # True for forward primer, False for reverse
    
    # These will be calculated later
    tm: float = 0.0
    gc_content: float = 0.0
    num_guanines: int = 0
    guanine_positions: List[int] = None  # Positions of Gs from 3' end
    has_gc_at_3prime: bool = False
    num_poly_bases: int = 0  # Max consecutive identical bases
    score: float = float('inf')  # Lower is better
    
    def __post_init__(self):
        if self.guanine_positions is None:
            self.guanine_positions = []


class PrimerGenerator:
    """Generate candidate primers from search regions"""
    
    def __init__(self, constraints: PrimerConstraints):
        """
        Initialize with primer constraints.
        
        Args:
            constraints: PrimerConstraints object
        """
        self.constraints = constraints
    
    def generate_primers(self, sequence: str, region_start: int, 
                        is_forward: bool, region_end: int = None) -> List[PrimerCandidate]:
        """
        Generate all possible primer candidates from a sequence region.
        Uses sliding window approach with pre-filtering.
        
        Args:
            sequence: The DNA sequence to search for primers
            region_start: Start position of this region in the template (for tracking)
            is_forward: True if generating forward primers, False for reverse
            region_end: End position of region (needed for reverse primers to calculate correct positions)
            
        Returns:
            List of PrimerCandidate objects that pass basic filters
        """
        candidates = []
        
        # Generate all possible k-mers
        for length in range(self.constraints.min_length, self.constraints.max_length + 1):
            for i in range(len(sequence) - length + 1):
                kmer = sequence[i:i + length]
                
                # Calculate position in original template
                if is_forward:
                    position = region_start + i
                else:
                    # For reverse primers, the sequence we're searching is the RC of the template region
                    # If region spans [region_start, region_end] on forward strand,
                    # a primer at position i in the RC corresponds to position (region_end - i - length) on forward strand
                    if region_end is None:
                        raise ValueError("region_end must be provided for reverse primers")
                    position = region_end - i - length
                
                # Quick pre-filtering before creating candidate
                if self._passes_basic_filters(kmer):
                    candidate = PrimerCandidate(
                        sequence=kmer,
                        position=position,
                        length=length,
                        is_forward=is_forward
                    )
                    
                    # Calculate basic properties
                    self._calculate_basic_properties(candidate)
                    
                    # Add if it passes all filters
                    if self._passes_all_basic_filters(candidate):
                        candidates.append(candidate)
        
        return candidates
    
    def _passes_basic_filters(self, sequence: str) -> bool:
        """
        Quick pre-filter before creating candidate object.
        Checks the most restrictive constraints first for efficiency.
        
        Args:
            sequence: Primer sequence
            
        Returns:
            True if passes basic filters
        """
        # Check for G/C at 3' end (most restrictive, check first)
        three_prime_end = sequence[-self.constraints.no_gc_3prime_bases:]
        if 'G' in three_prime_end or 'C' in three_prime_end:
            return False
        
        # Count Gs (second most restrictive)
        num_g = sequence.count('G')
        if num_g > self.constraints.max_guanines:
            return False
        
        # Check for excessive poly-bases
        if self._get_max_poly_base(sequence) > self.constraints.max_poly_base:
            return False
        
        return True
    
    def _calculate_basic_properties(self, candidate: PrimerCandidate):
        """
        Calculate basic properties of a primer candidate.
        
        Args:
            candidate: PrimerCandidate object to update
        """
        seq = candidate.sequence
        
        # Count Gs and record positions (from 3' end)
        candidate.num_guanines = seq.count('G')
        candidate.guanine_positions = []
        for i, base in enumerate(reversed(seq)):
            if base == 'G':
                candidate.guanine_positions.append(i)
        
        # Check for G/C at 3' end
        three_prime = seq[-self.constraints.no_gc_3prime_bases:]
        candidate.has_gc_at_3prime = 'G' in three_prime or 'C' in three_prime
        
        # Calculate GC content
        gc_count = seq.count('G') + seq.count('C')
        candidate.gc_content = (gc_count / len(seq)) * 100
        
        # Find longest poly-base run
        candidate.num_poly_bases = self._get_max_poly_base(seq)
    
    def _get_max_poly_base(self, sequence: str) -> int:
        """
        Find the longest run of consecutive identical bases.
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Length of longest poly-base run
        """
        if not sequence:
            return 0
        
        max_run = 1
        current_run = 1
        
        for i in range(1, len(sequence)):
            if sequence[i] == sequence[i-1]:
                current_run += 1
                max_run = max(max_run, current_run)
            else:
                current_run = 1
        
        return max_run
    
    def _passes_all_basic_filters(self, candidate: PrimerCandidate) -> bool:
        """
        Check if candidate passes all basic filters.
        
        Args:
            candidate: PrimerCandidate object
            
        Returns:
            True if passes all filters
        """
        # Check G/C at 3' end
        if candidate.has_gc_at_3prime:
            return False
        
        # Check number of Gs
        if candidate.num_guanines > self.constraints.max_guanines:
            return False
        
        # Check GC content
        if (candidate.gc_content < self.constraints.min_gc_content or 
            candidate.gc_content > self.constraints.max_gc_content):
            return False
        
        # Check poly-bases
        if candidate.num_poly_bases > self.constraints.max_poly_base:
            return False
        
        return True
    
    def filter_by_tm_range(self, candidates: List[PrimerCandidate]) -> List[PrimerCandidate]:
        """
        Filter candidates by Tm range.
        Note: Tm must be calculated first (done in primer_evaluator module).
        
        Args:
            candidates: List of PrimerCandidate objects with Tm calculated
            
        Returns:
            Filtered list of candidates within Tm range
        """
        filtered = []
        for candidate in candidates:
            if self.constraints.tm_min <= candidate.tm <= self.constraints.tm_max:
                filtered.append(candidate)
        return filtered
    
    def get_statistics(self, candidates: List[PrimerCandidate]) -> Dict:
        """
        Get statistics about generated primers.
        
        Args:
            candidates: List of PrimerCandidate objects
            
        Returns:
            Dictionary with statistics
        """
        if not candidates:
            return {
                'total': 0,
                'forward': 0,
                'reverse': 0
            }
        
        forward_primers = [c for c in candidates if c.is_forward]
        reverse_primers = [c for c in candidates if not c.is_forward]
        
        stats = {
            'total': len(candidates),
            'forward': len(forward_primers),
            'reverse': len(reverse_primers),
            'avg_length': sum(c.length for c in candidates) / len(candidates),
            'avg_gc_content': sum(c.gc_content for c in candidates) / len(candidates),
            'avg_num_guanines': sum(c.num_guanines for c in candidates) / len(candidates),
            'primers_with_0_g': len([c for c in candidates if c.num_guanines == 0]),
            'primers_with_1_g': len([c for c in candidates if c.num_guanines == 1]),
            'primers_with_2_g': len([c for c in candidates if c.num_guanines == 2]),
        }
        
        return stats


def test_primer_generator():
    """Test function for PrimerGenerator"""
    from config import PrimerConstraints
    
    # Create test constraints
    constraints = PrimerConstraints(
        min_length=18,
        max_length=22,
        max_guanines=2,
        no_gc_3prime_bases=3,
        min_gc_content=20.0,
        max_gc_content=60.0,
        max_poly_base=4
    )
    
    generator = PrimerGenerator(constraints)
    
    # Test sequence - mix of good and bad primers
    test_seq = "ATATATATATATATATATAT" + "GGGGG" + "ATATATATATAT" + "ATACATATATA"
    #             ^ Good region        ^ Bad (poly-G)  ^ Good       ^ Good
    
    print(f"Test sequence length: {len(test_seq)}")
    print(f"Test sequence: {test_seq}\n")
    
    # Generate primers
    candidates = generator.generate_primers(test_seq, region_start=1000, is_forward=True)
    
    print(f"Generated {len(candidates)} candidate primers")
    
    # Show some examples
    print("\nFirst 5 primers:")
    for i, primer in enumerate(candidates[:5]):
        print(f"{i+1}. Pos {primer.position}: {primer.sequence}")
        print(f"   Length: {primer.length}, Gs: {primer.num_guanines}, GC%: {primer.gc_content:.1f}")
        print(f"   G positions from 3' end: {primer.guanine_positions}")
    
    # Statistics
    stats = generator.get_statistics(candidates)
    print("\nStatistics:")
    for key, value in stats.items():
        if isinstance(value, float):
            print(f"  {key}: {value:.2f}")
        else:
            print(f"  {key}: {value}")


if __name__ == "__main__":
    test_primer_generator()