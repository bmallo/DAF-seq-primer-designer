"""
Primer Evaluator Module
Calculates Tm, checks secondary structures, and scores primers
"""

from typing import List, Tuple
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from primer_generator import PrimerCandidate
from config import PrimerConstraints, ScoringWeights


class PrimerEvaluator:
    """Evaluate primers: calculate Tm, check structures, and score"""
    
    def __init__(self, constraints: PrimerConstraints, weights: ScoringWeights):
        """
        Initialize with constraints and scoring weights.
        
        Args:
            constraints: PrimerConstraints object
            weights: ScoringWeights object
        """
        self.constraints = constraints
        self.weights = weights
    
    def calculate_tm(self, primer: PrimerCandidate) -> float:
        """
        Calculate melting temperature using Primer3 basic method.
        
        Args:
            primer: PrimerCandidate object
            
        Returns:
            Melting temperature in Celsius
        """
        seq = Seq(primer.sequence)
        
        # Use salt-adjusted method (similar to Primer3 default)
        # This accounts for salt concentrations
        tm = mt.Tm_NN(
            seq,
            Na=self.constraints.Na,
            K=self.constraints.K,
            Mg=self.constraints.Mg,
            dNTPs=self.constraints.dNTPs,
            saltcorr=7  # Owczarzy et al. 2008 correction
        )
        
        return round(tm, 2)
    
    def check_hairpin(self, primer: PrimerCandidate, min_stem_length: int = 4,
                     max_loop_size: int = 8) -> Tuple[bool, int]:
        """
        Check for potential hairpin formation.
        A hairpin forms when a primer's 3' end can fold back on itself.
        
        Args:
            primer: PrimerCandidate object
            min_stem_length: Minimum number of complementary bases for hairpin stem
            max_loop_size: Maximum size of hairpin loop
            
        Returns:
            Tuple of (has_hairpin, stem_length)
        """
        seq = primer.sequence
        rev_comp = str(Seq(seq).reverse_complement())
        
        # Check for complementarity that could form hairpins
        # We're looking for regions where the sequence is complementary to itself
        max_stem = 0
        
        # Check different possible loop sizes
        for loop_size in range(3, max_loop_size + 1):
            # For each position, check if it can form a stem
            for i in range(len(seq) - loop_size - min_stem_length):
                # Potential 5' stem region
                stem_5prime = seq[i:i + min_stem_length]
                # Potential 3' stem region (after the loop)
                start_3prime = i + loop_size + min_stem_length
                
                if start_3prime + min_stem_length <= len(seq):
                    stem_3prime = seq[start_3prime:start_3prime + min_stem_length]
                    stem_3prime_rc = str(Seq(stem_3prime).reverse_complement())
                    
                    # Count complementary bases
                    matches = sum(1 for a, b in zip(stem_5prime, stem_3prime_rc) if a == b)
                    
                    if matches >= min_stem_length:
                        max_stem = max(max_stem, matches)
        
        has_hairpin = max_stem >= min_stem_length
        return has_hairpin, max_stem
    
    def check_self_dimer(self, primer: PrimerCandidate, min_complementarity: int = 4) -> Tuple[bool, int]:
        """
        Check for self-dimer formation (primer binding to itself).
        Particularly important to check 3' end complementarity.
        
        Args:
            primer: PrimerCandidate object
            min_complementarity: Minimum consecutive complementary bases to flag
            
        Returns:
            Tuple of (has_self_dimer, max_complementarity)
        """
        seq = primer.sequence
        rev_comp = str(Seq(seq).reverse_complement())
        
        max_comp = 0
        
        # Check for complementarity when primer aligns with its reverse complement
        # This simulates two primers binding to each other
        for offset in range(-len(seq) + 1, len(seq)):
            consecutive_matches = 0
            current_run = 0
            
            for i in range(len(seq)):
                j = i + offset
                if 0 <= j < len(seq):
                    if seq[i] == rev_comp[j]:
                        current_run += 1
                        consecutive_matches = max(consecutive_matches, current_run)
                    else:
                        current_run = 0
            
            max_comp = max(max_comp, consecutive_matches)
        
        has_self_dimer = max_comp >= min_complementarity
        return has_self_dimer, max_comp
    
    def check_three_prime_self_dimer(self, primer: PrimerCandidate, check_length: int = 5) -> Tuple[bool, int]:
        """
        Specifically check if the 3' end can bind to itself or another copy of the primer.
        This is particularly problematic for PCR.
        
        Args:
            primer: PrimerCandidate object
            check_length: Number of bases from 3' end to check
            
        Returns:
            Tuple of (has_3prime_dimer, complementarity_score)
        """
        seq = primer.sequence
        three_prime = seq[-check_length:]
        three_prime_rc = str(Seq(three_prime).reverse_complement())
        
        # Check if 3' end is complementary to any part of the primer
        max_match = 0
        for i in range(len(seq) - check_length + 1):
            region = seq[i:i + check_length]
            matches = sum(1 for a, b in zip(region, three_prime_rc) if a == b)
            max_match = max(max_match, matches)
        
        # High complementarity at 3' end is problematic
        has_3prime_dimer = max_match >= check_length - 1
        return has_3prime_dimer, max_match
    
    def score_primer(self, primer: PrimerCandidate) -> float:
        """
        Calculate a composite score for the primer.
        Lower score is better (0 = perfect primer).
        
        Args:
            primer: PrimerCandidate object with all properties calculated
            
        Returns:
            Score (lower is better)
        """
        score = 0.0
        
        # 1. Guanine penalty (most important for bisulfite primers)
        for g_pos in primer.guanine_positions:
            position_penalty = self.weights.get_g_position_penalty(g_pos, primer.length)
            score += position_penalty
        
        # 2. Tm deviation penalty
        optimal_tm = (self.constraints.tm_min + self.constraints.tm_max) / 2
        tm_deviation = abs(primer.tm - optimal_tm)
        score += tm_deviation * self.weights.tm_deviation_penalty
        
        # 3. GC content deviation penalty
        optimal_gc = (self.constraints.min_gc_content + self.constraints.max_gc_content) / 2
        gc_deviation = abs(primer.gc_content - optimal_gc)
        score += gc_deviation * self.weights.gc_deviation_penalty
        
        # 4. Length deviation penalty
        optimal_length = (self.constraints.min_length + self.constraints.max_length) / 2
        length_deviation = abs(primer.length - optimal_length)
        score += length_deviation * self.weights.length_deviation_penalty
        
        # 5. Poly-base penalty
        if primer.num_poly_bases > self.constraints.max_poly_base:
            excess_poly = primer.num_poly_bases - self.constraints.max_poly_base
            score += excess_poly * self.weights.poly_base_penalty
        
        # 6. Secondary structure penalties (these are stored as attributes)
        # These will be added if structures are found
        
        return round(score, 2)
    
    def evaluate_primers(self, primers: List[PrimerCandidate], verbose: bool = True) -> List[PrimerCandidate]:
        """
        Evaluate all primers: calculate Tm, check structures, and score.
        
        Args:
            primers: List of PrimerCandidate objects
            verbose: Whether to print progress messages
            
        Returns:
            List of evaluated primers (same objects, updated with scores)
        """
        if verbose:
            print(f"Evaluating {len(primers)} primers...")
        
        for i, primer in enumerate(primers):
            # Calculate Tm
            primer.tm = self.calculate_tm(primer)
            
            # Check secondary structures
            has_hairpin, hairpin_stem = self.check_hairpin(primer)
            has_self_dimer, self_dimer_comp = self.check_self_dimer(primer)
            has_3prime_dimer, three_prime_comp = self.check_three_prime_self_dimer(primer)
            
            # Add structure penalties to score
            structure_penalty = 0.0
            if has_hairpin:
                structure_penalty += self.weights.hairpin_penalty * hairpin_stem
            if has_self_dimer:
                structure_penalty += self.weights.self_dimer_penalty * self_dimer_comp
            if has_3prime_dimer:
                structure_penalty += self.weights.self_dimer_penalty * three_prime_comp * 2  # 3' dimers are worse
            
            # Calculate base score
            primer.score = self.score_primer(primer)
            
            # Add structure penalties
            primer.score += structure_penalty
            
            # Progress indicator for large datasets
            if verbose and (i + 1) % 1000 == 0:
                print(f"  Evaluated {i + 1}/{len(primers)} primers...")
        
        if verbose:
            print(f"Evaluation complete!")
        return primers
    
    def filter_by_tm(self, primers: List[PrimerCandidate]) -> List[PrimerCandidate]:
        """
        Filter primers by Tm range.
        
        Args:
            primers: List of PrimerCandidate objects
            
        Returns:
            Filtered list
        """
        return [p for p in primers if self.constraints.tm_min <= p.tm <= self.constraints.tm_max]
    
    def get_best_primers(self, primers: List[PrimerCandidate], n: int = 100) -> List[PrimerCandidate]:
        """
        Get the top N primers by score.
        
        Args:
            primers: List of PrimerCandidate objects
            n: Number of primers to return
            
        Returns:
            List of best primers, sorted by score
        """
        sorted_primers = sorted(primers, key=lambda p: p.score)
        return sorted_primers[:n]


def test_primer_evaluator():
    """Test function for PrimerEvaluator"""
    from config import PrimerConstraints, ScoringWeights
    from primer_generator import PrimerCandidate
    
    # Create test constraints
    constraints = PrimerConstraints(
        min_length=18,
        max_length=25,
        tm_min=58.0,
        tm_max=62.0,
        Na=50.0
    )
    
    weights = ScoringWeights()
    
    evaluator = PrimerEvaluator(constraints, weights)
    
    # Create test primers
    test_primers = [
        PrimerCandidate("ATATATATATATATATATAT", 100, 20, True),  # No Gs, good
        PrimerCandidate("ATATATATGATATATATAT", 200, 19, True),    # 1 G in middle
        PrimerCandidate("ATATATATATATGATATAT", 300, 19, True),    # 1 G closer to 3'
        PrimerCandidate("ATGATGATATATATATATAT", 400, 20, True),   # 2 Gs at 5' end
    ]
    
    # Calculate properties that would normally be done by generator
    from primer_generator import PrimerGenerator
    generator = PrimerGenerator(constraints)
    for primer in test_primers:
        generator._calculate_basic_properties(primer)
    
    # Evaluate primers
    evaluated = evaluator.evaluate_primers(test_primers)
    
    print("\nEvaluated Primers:")
    for i, primer in enumerate(evaluated, 1):
        print(f"\n{i}. {primer.sequence}")
        print(f"   Position: {primer.position}")
        print(f"   Tm: {primer.tm}°C")
        print(f"   GC%: {primer.gc_content:.1f}%")
        print(f"   Gs: {primer.num_guanines} at positions {primer.guanine_positions} (from 3')")
        print(f"   Score: {primer.score}")
    
    # Get best primers
    best = evaluator.get_best_primers(evaluated, n=2)
    print(f"\nBest 2 primers:")
    for i, primer in enumerate(best, 1):
        print(f"{i}. Score {primer.score}: {primer.sequence}")


if __name__ == "__main__":
    test_primer_evaluator()