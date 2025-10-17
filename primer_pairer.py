"""
Primer Pairing Module
Pairs forward and reverse primers and evaluates primer sets
"""

from typing import List, Tuple, Dict
from dataclasses import dataclass
from collections import Counter
from Bio.Seq import Seq
from primer_generator import PrimerCandidate
from config import PairingConstraints, ScoringWeights, TargetRegion


@dataclass
class PrimerPair:
    """Represents a pair of forward and reverse primers"""
    forward: PrimerCandidate
    reverse: PrimerCandidate
    amplicon_size: int
    tm_difference: float
    pair_score: float = 0.0
    has_primer_dimer: bool = False
    primer_dimer_score: int = 0
    covers_core: bool = True
    
    def __str__(self):
        return (f"PrimerPair(amplicon={self.amplicon_size}bp, "
                f"ΔTm={self.tm_difference:.2f}°C, score={self.pair_score:.2f})")


class PrimerPairer:
    """Pair forward and reverse primers into compatible sets"""
    
    def __init__(self, pairing_constraints: PairingConstraints, 
                 scoring_weights: ScoringWeights,
                 target_region: TargetRegion):
        """
        Initialize with pairing constraints.
        
        Args:
            pairing_constraints: PairingConstraints object
            scoring_weights: ScoringWeights object  
            target_region: TargetRegion object for validation
        """
        self.constraints = pairing_constraints
        self.weights = scoring_weights
        self.target_region = target_region
    
    def check_primer_dimer(self, fwd: PrimerCandidate, rev: PrimerCandidate, 
                          min_complementarity: int = 4) -> Tuple[bool, int]:
        """
        Check if forward and reverse primers can form dimers.
        Particularly important to check 3' end complementarity.
        
        Args:
            fwd: Forward primer
            rev: Reverse primer
            min_complementarity: Minimum consecutive complementary bases
            
        Returns:
            Tuple of (has_dimer, max_complementarity)
        """
        fwd_seq = fwd.sequence
        rev_seq_rc = str(Seq(rev.sequence).reverse_complement())
        
        max_comp = 0
        max_3prime_comp = 0
        
        # Check all possible alignments
        for offset in range(-len(fwd_seq) + 1, len(rev_seq_rc)):
            consecutive_matches = 0
            current_run = 0
            three_prime_matches = 0
            
            for i in range(len(fwd_seq)):
                j = i + offset
                if 0 <= j < len(rev_seq_rc):
                    if fwd_seq[i] == rev_seq_rc[j]:
                        current_run += 1
                        consecutive_matches = max(consecutive_matches, current_run)
                        
                        # Check if this involves 3' ends (last 5 bases)
                        if i >= len(fwd_seq) - 5 or j >= len(rev_seq_rc) - 5:
                            three_prime_matches += 1
                    else:
                        current_run = 0
            
            max_comp = max(max_comp, consecutive_matches)
            max_3prime_comp = max(max_3prime_comp, three_prime_matches)
        
        # Primer dimer is problematic if there's significant complementarity,
        # especially at 3' ends
        has_dimer = max_comp >= min_complementarity or max_3prime_comp >= 3
        
        return has_dimer, max(max_comp, max_3prime_comp * 2)  # Weight 3' comp higher
    
    def calculate_amplicon_size(self, fwd: PrimerCandidate, rev: PrimerCandidate) -> int:
        """
        Calculate the amplicon size for a primer pair.
        
        Args:
            fwd: Forward primer
            rev: Reverse primer
            
        Returns:
            Amplicon size in bp
        """
        # Amplicon goes from start of forward primer to end of reverse primer
        amplicon_size = (rev.position + rev.length) - fwd.position
        return amplicon_size
    
    def check_amplicon_covers_core(self, fwd: PrimerCandidate, rev: PrimerCandidate) -> bool:
        """
        Check if the amplicon covers the entire core region.
        
        Args:
            fwd: Forward primer
            rev: Reverse primer
            
        Returns:
            True if amplicon covers core region
        """
        # The amplicon is the region BETWEEN the primers (not including them)
        amplicon_start = fwd.position + fwd.length
        amplicon_end = rev.position
        
        core_start = self.target_region.core_start
        core_end = self.target_region.core_end
        
        return amplicon_start <= core_start and amplicon_end >= core_end
    
    def score_pair(self, pair: PrimerPair) -> float:
        """
        Calculate a score for the primer pair.
        Lower is better.
        
        Args:
            pair: PrimerPair object
            
        Returns:
            Pair score (lower is better)
        """
        # Base score is sum of individual primer scores
        score = pair.forward.score + pair.reverse.score
        
        # Penalty for Tm difference
        score += abs(pair.tm_difference) * self.weights.tm_deviation_penalty * 5  # 5x weight for pair mismatch
        
        # Penalty for primer dimers
        if pair.has_primer_dimer:
            score += self.weights.self_dimer_penalty * pair.primer_dimer_score
        
        # Penalty for amplicon size deviation from optimal
        optimal_amplicon = (self.target_region.min_amplicon + self.target_region.max_amplicon) / 2
        amplicon_deviation = abs(pair.amplicon_size - optimal_amplicon)
        score += amplicon_deviation * 0.01  # Small penalty for size deviation
        
        # Penalty if doesn't cover core
        if not pair.covers_core:
            score += 1000  # Large penalty
        
        return round(score, 2)
    
    def create_pairs(self, forward_primers: List[PrimerCandidate], 
                     reverse_primers: List[PrimerCandidate],
                     verbose: bool = True) -> List[PrimerPair]:
        """
        Create all valid primer pairs from forward and reverse primer lists.
        
        Args:
            forward_primers: List of forward primers
            reverse_primers: List of reverse primers
            verbose: Whether to print progress messages
            
        Returns:
            List of PrimerPair objects
        """
        if verbose:
            print(f"\nCreating primer pairs from {len(forward_primers)} forward and {len(reverse_primers)} reverse primers...")
        
        pairs = []
        total_combinations = len(forward_primers) * len(reverse_primers)
        
        for i, fwd in enumerate(forward_primers):
            for rev in reverse_primers:
                # Calculate amplicon size
                amplicon_size = self.calculate_amplicon_size(fwd, rev)
                
                # Check if amplicon size is within range
                if not (self.target_region.min_amplicon <= amplicon_size <= self.target_region.max_amplicon):
                    continue
                
                # Check if amplicon covers core region
                covers_core = self.check_amplicon_covers_core(fwd, rev)
                if not covers_core:
                    continue
                
                # Calculate Tm difference
                tm_diff = abs(fwd.tm - rev.tm)
                
                # Check if Tm difference is acceptable
                if tm_diff > self.constraints.max_tm_difference:
                    continue
                
                # Check for primer dimers
                has_dimer, dimer_score = self.check_primer_dimer(fwd, rev)
                
                # Create pair
                pair = PrimerPair(
                    forward=fwd,
                    reverse=rev,
                    amplicon_size=amplicon_size,
                    tm_difference=tm_diff,
                    has_primer_dimer=has_dimer,
                    primer_dimer_score=dimer_score,
                    covers_core=covers_core
                )
                
                # Score the pair
                pair.pair_score = self.score_pair(pair)
                
                pairs.append(pair)
            
            # Progress indicator
            if verbose and (i + 1) % 100 == 0:
                print(f"  Processed {i + 1}/{len(forward_primers)} forward primers...")
        
        if verbose:
            print(f"Created {len(pairs)} valid primer pairs from {total_combinations} possible combinations")
        return pairs
    
    def filter_pairs_by_reuse(self, pairs: List[PrimerPair]) -> List[PrimerPair]:
        """
        Filter pairs to ensure no primer is reused too many times.
        Prioritizes best-scoring pairs.
        
        Args:
            pairs: List of PrimerPair objects, should be sorted by score
            
        Returns:
            Filtered list of pairs
        """
        if not pairs:
            return []
        
        # Sort pairs by score (best first)
        sorted_pairs = sorted(pairs, key=lambda p: p.pair_score)
        
        # Track how many times each primer is used
        fwd_usage = Counter()
        rev_usage = Counter()
        
        filtered_pairs = []
        
        for pair in sorted_pairs:
            fwd_id = (pair.forward.position, pair.forward.sequence)
            rev_id = (pair.reverse.position, pair.reverse.sequence)
            
            # Check if we can still use these primers
            if (fwd_usage[fwd_id] < self.constraints.max_primer_reuse and
                rev_usage[rev_id] < self.constraints.max_primer_reuse):
                
                filtered_pairs.append(pair)
                fwd_usage[fwd_id] += 1
                rev_usage[rev_id] += 1
        
        return filtered_pairs
    
    def get_best_pairs(self, pairs: List[PrimerPair], n: int = 10, 
                       apply_reuse_filter: bool = True,
                       enforce_spacing: bool = True) -> List[PrimerPair]:
        """
        Get the best N primer pairs with diversity constraints.
        
        Args:
            pairs: List of PrimerPair objects
            n: Number of pairs to return
            apply_reuse_filter: Whether to limit primer reuse
            enforce_spacing: Whether to enforce minimum spacing between primers
            
        Returns:
            List of best primer pairs, sorted by score
        """
        if not pairs:
            return []
        
        # Sort by score
        sorted_pairs = sorted(pairs, key=lambda p: p.pair_score)
        
        if not enforce_spacing:
            # Simple filtering without spacing constraint
            if apply_reuse_filter:
                sorted_pairs = self.filter_pairs_by_reuse(sorted_pairs)
            return sorted_pairs[:n]
        
        # Apply spacing constraint for diverse primer sets
        selected_pairs = []
        used_forward_positions = []  # Track positions of selected forward primers
        used_reverse_positions = []  # Track positions of selected reverse primers
        
        for pair in sorted_pairs:
            if len(selected_pairs) >= n:
                break
            
            # Check if this primer is too close to any already selected primers
            fwd_too_close = False
            rev_too_close = False
            
            for used_fwd_pos in used_forward_positions:
                if abs(pair.forward.position - used_fwd_pos) < self.constraints.min_primer_spacing:
                    fwd_too_close = True
                    break
            
            for used_rev_pos in used_reverse_positions:
                if abs(pair.reverse.position - used_rev_pos) < self.constraints.min_primer_spacing:
                    rev_too_close = True
                    break
            
            # If primers are sufficiently different from already selected ones, add this pair
            if not fwd_too_close and not rev_too_close:
                selected_pairs.append(pair)
                used_forward_positions.append(pair.forward.position)
                used_reverse_positions.append(pair.reverse.position)
        
        return selected_pairs
    
    def get_pair_statistics(self, pairs: List[PrimerPair]) -> Dict:
        """
        Get statistics about primer pairs.
        
        Args:
            pairs: List of PrimerPair objects
            
        Returns:
            Dictionary with statistics
        """
        if not pairs:
            return {'total_pairs': 0}
        
        return {
            'total_pairs': len(pairs),
            'avg_amplicon_size': sum(p.amplicon_size for p in pairs) / len(pairs),
            'min_amplicon_size': min(p.amplicon_size for p in pairs),
            'max_amplicon_size': max(p.amplicon_size for p in pairs),
            'avg_tm_difference': sum(p.tm_difference for p in pairs) / len(pairs),
            'pairs_with_dimers': sum(1 for p in pairs if p.has_primer_dimer),
            'avg_pair_score': sum(p.pair_score for p in pairs) / len(pairs),
            'pairs_covering_core': sum(1 for p in pairs if p.covers_core),
        }


def test_primer_pairer():
    """Test function for PrimerPairer"""
    from config import PairingConstraints, ScoringWeights, TargetRegion
    from primer_generator import PrimerCandidate
    
    # Create test constraints
    pairing_constraints = PairingConstraints(
        max_tm_difference=2.0,
        max_primer_reuse=3,
        min_primer_spacing=50
    )
    
    scoring_weights = ScoringWeights()
    
    target_region = TargetRegion(
        core_start=2000,
        core_end=3000,
        min_amplicon=2000,
        max_amplicon=5000
    )
    
    pairer = PrimerPairer(pairing_constraints, scoring_weights, target_region)
    
    # Create test primers
    fwd_primers = [
        PrimerCandidate("ATATATATATATATATATAT", 1000, 20, True),
        PrimerCandidate("TATATATATATATATATAAA", 1100, 20, True),
    ]
    
    rev_primers = [
        PrimerCandidate("ATATATATATATATATATAT", 4000, 20, False),
        PrimerCandidate("TATATATATATATATATAAA", 4500, 20, False),
    ]
    
    # Set mock Tm values
    for primer in fwd_primers + rev_primers:
        primer.tm = 60.0
        primer.score = 10.0
    
    # Create pairs
    pairs = pairer.create_pairs(fwd_primers, rev_primers)
    
    print(f"\nCreated {len(pairs)} pairs")
    
    # Show pairs
    for i, pair in enumerate(pairs, 1):
        print(f"\n{i}. {pair}")
        print(f"   Forward: {pair.forward.sequence} at {pair.forward.position}")
        print(f"   Reverse: {pair.reverse.sequence} at {pair.reverse.position}")
        print(f"   Amplicon: {pair.amplicon_size} bp")
        print(f"   Tm diff: {pair.tm_difference:.2f}°C")
        print(f"   Score: {pair.pair_score:.2f}")
    
    # Statistics
    stats = pairer.get_pair_statistics(pairs)
    print("\nStatistics:")
    for key, value in stats.items():
        if isinstance(value, float):
            print(f"  {key}: {value:.2f}")
        else:
            print(f"  {key}: {value}")


if __name__ == "__main__":
    test_primer_pairer()