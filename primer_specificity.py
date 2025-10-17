"""
Primer Specificity Checker Module
Checks for off-target binding sites within the template sequence
"""

from typing import List, Tuple, Dict
from Bio.Seq import Seq
from primer_generator import PrimerCandidate


class SpecificityChecker:
    """Check primer specificity against template sequence"""
    
    def __init__(self, template_sequence: str, min_match_length: int = 15,
                 allow_mismatches: int = 2):
        """
        Initialize specificity checker.
        
        Args:
            template_sequence: The full template DNA sequence
            min_match_length: Minimum length of match to consider as potential binding
            allow_mismatches: Number of mismatches allowed in binding site
        """
        self.template = template_sequence.upper()
        self.template_rc = str(Seq(template_sequence).reverse_complement())
        self.min_match_length = min_match_length
        self.allow_mismatches = allow_mismatches
    
    def find_exact_matches(self, primer_seq: str) -> List[int]:
        """
        Find all exact matches of primer sequence in template.
        
        Args:
            primer_seq: Primer sequence to search for
            
        Returns:
            List of positions where exact matches are found
        """
        positions = []
        primer_seq = primer_seq.upper()
        
        # Search forward strand
        start = 0
        while True:
            pos = self.template.find(primer_seq, start)
            if pos == -1:
                break
            positions.append(pos)
            start = pos + 1
        
        # Search reverse complement strand
        start = 0
        while True:
            pos = self.template_rc.find(primer_seq, start)
            if pos == -1:
                break
            # Convert RC position back to forward strand position
            rc_pos = len(self.template) - pos - len(primer_seq)
            positions.append(rc_pos)
            start = pos + 1
        
        return sorted(set(positions))  # Remove duplicates and sort
    
    def find_partial_matches(self, primer_seq: str, expected_position: int) -> List[Tuple[int, int, int]]:
        """
        Find partial matches allowing mismatches.
        Only searches for matches that are NOT at the expected position.
        
        Args:
            primer_seq: Primer sequence to search for
            expected_position: The intended binding position (to exclude from off-targets)
            
        Returns:
            List of tuples: (position, num_mismatches, match_length)
        """
        partial_matches = []
        primer_seq = primer_seq.upper()
        primer_len = len(primer_seq)
        
        # Check both strands
        for strand_seq, is_rc in [(self.template, False), (self.template_rc, True)]:
            # Slide primer across template
            for pos in range(len(strand_seq) - self.min_match_length + 1):
                # Skip if this is the expected position on forward strand
                if not is_rc and pos == expected_position:
                    continue
                
                # Check for match with mismatches
                match_len = min(primer_len, len(strand_seq) - pos)
                if match_len < self.min_match_length:
                    continue
                
                template_segment = strand_seq[pos:pos + match_len]
                primer_segment = primer_seq[:match_len]
                
                # Count mismatches
                mismatches = sum(1 for a, b in zip(primer_segment, template_segment) if a != b)
                
                # If within allowed mismatches, record as potential binding site
                if mismatches <= self.allow_mismatches:
                    # Convert RC position back to forward strand if necessary
                    actual_pos = len(self.template) - pos - match_len if is_rc else pos
                    partial_matches.append((actual_pos, mismatches, match_len))
        
        return partial_matches
    
    def check_3prime_specificity(self, primer_seq: str, expected_position: int,
                                 check_length: int = 10) -> List[Tuple[int, int]]:
        """
        Check if the 3' end of the primer matches elsewhere in the template.
        This is critical because PCR extension starts from the 3' end.
        
        Args:
            primer_seq: Primer sequence
            expected_position: The intended binding position
            check_length: Length of 3' end to check
            
        Returns:
            List of tuples: (position, num_mismatches) for off-target 3' matches
        """
        three_prime_end = primer_seq[-check_length:].upper()
        off_target_matches = []
        
        # Check both strands
        for strand_seq, is_rc in [(self.template, False), (self.template_rc, True)]:
            start = 0
            while True:
                pos = strand_seq.find(three_prime_end, start)
                if pos == -1:
                    break
                
                # Convert position if reverse complement
                actual_pos = len(self.template) - pos - check_length if is_rc else pos
                
                # Skip if this is the expected position
                if not is_rc and actual_pos == expected_position + len(primer_seq) - check_length:
                    start = pos + 1
                    continue
                
                # Check how well the full primer matches here
                full_primer_start = max(0, actual_pos - (len(primer_seq) - check_length))
                template_segment = strand_seq[pos - (len(primer_seq) - check_length):pos + check_length]
                
                if len(template_segment) >= len(primer_seq):
                    mismatches = sum(1 for a, b in zip(primer_seq, template_segment) if a != b)
                    off_target_matches.append((actual_pos, mismatches))
                
                start = pos + 1
        
        return off_target_matches
    
    def is_specific(self, primer: PrimerCandidate, allow_multiple_exact: bool = False) -> Tuple[bool, Dict]:
        """
        Determine if a primer is specific to its intended binding site.
        
        Args:
            primer: PrimerCandidate object
            allow_multiple_exact: If False, reject primers with multiple exact matches
            
        Returns:
            Tuple of (is_specific, details_dict)
        """
        details = {
            'exact_matches': [],
            'partial_matches': [],
            'three_prime_matches': [],
            'is_specific': True,
            'reason': ''
        }
        
        # For reverse primers, we need to search for the reverse complement
        # because the primer sequence is 5'->3' but binds to the minus strand
        search_seq = primer.sequence
        if not primer.is_forward:
            search_seq = str(Seq(primer.sequence).reverse_complement())
        
        # Find exact matches
        exact_matches = self.find_exact_matches(search_seq)
        details['exact_matches'] = exact_matches
        
        # For reverse primers, the expected position should match where the RC of the primer binds
        # The primer.position is where the primer binds on the forward strand
        expected_position = primer.position
        
        # Check if expected position is in exact matches
        expected_in_matches = expected_position in exact_matches
        
        # For debugging - if expected position not found, try to be more lenient
        # Check if any match is within primer length of expected position
        close_matches = [pos for pos in exact_matches if abs(pos - expected_position) < len(primer.sequence)]
        
        if not allow_multiple_exact and len(exact_matches) > 1:
            # Check if the "extra" matches are just the same site detected differently
            if len(close_matches) <= 1:
                # All other matches are far away, this is truly non-specific
                details['is_specific'] = False
                details['reason'] = f"Multiple exact matches found at positions: {exact_matches}"
                return False, details
        
        # If we have exact matches but none near the expected position, that's a problem
        if len(exact_matches) > 0 and len(close_matches) == 0:
            details['is_specific'] = False
            details['reason'] = f"Exact matches found but not at expected position {expected_position}: {exact_matches}"
            return False, details
        
        # Find partial matches (off-targets)
        partial_matches = self.find_partial_matches(search_seq, expected_position)
        details['partial_matches'] = partial_matches
        
        if len(partial_matches) > 0:
            details['is_specific'] = False
            details['reason'] = f"Found {len(partial_matches)} potential off-target binding sites"
            return False, details
        
        # Check 3' end specificity (most critical)
        three_prime_matches = self.check_3prime_specificity(search_seq, expected_position)
        details['three_prime_matches'] = three_prime_matches
        
        if len(three_prime_matches) > 0:
            details['is_specific'] = False
            details['reason'] = f"3' end matches found at {len(three_prime_matches)} off-target sites"
            return False, details
        
        # Primer is specific
        details['reason'] = "Specific to intended binding site"
        return True, details
    
    def filter_specific_primers(self, primers: List[PrimerCandidate],
                               verbose: bool = True, debug: bool = False) -> List[PrimerCandidate]:
        """
        Filter primers to keep only those that are specific.
        
        Args:
            primers: List of PrimerCandidate objects
            verbose: Print statistics
            debug: Print detailed debug information for first few primers
            
        Returns:
            Filtered list of specific primers
        """
        if verbose:
            print(f"Checking specificity for {len(primers)} primers...")
        
        specific_primers = []
        rejection_reasons = {
            'multiple_exact': 0,
            'partial_matches': 0,
            'three_prime_matches': 0,
            'unexpected_position': 0
        }
        
        for i, primer in enumerate(primers):
            is_specific, details = self.is_specific(primer)
            
            # Debug first few primers
            if debug and i < 3:
                print(f"\n  Debug primer {i+1}:")
                print(f"    Sequence: {primer.sequence}")
                print(f"    Position: {primer.position}")
                print(f"    Is forward: {primer.is_forward}")
                print(f"    Exact matches: {details['exact_matches']}")
                print(f"    Is specific: {is_specific}")
                print(f"    Reason: {details['reason']}")
            
            if is_specific:
                specific_primers.append(primer)
            else:
                # Track rejection reasons
                if 'Multiple exact matches' in details['reason']:
                    rejection_reasons['multiple_exact'] += 1
                elif 'unexpected position' in details['reason']:
                    rejection_reasons['unexpected_position'] += 1
                elif '3\' end matches' in details['reason']:
                    rejection_reasons['three_prime_matches'] += 1
                elif 'off-target' in details['reason']:
                    rejection_reasons['partial_matches'] += 1
            
            # Progress indicator
            if verbose and (i + 1) % 500 == 0:
                print(f"  Checked {i + 1}/{len(primers)} primers...")
        
        if verbose:
            print(f"\nSpecificity filtering results:")
            print(f"  Specific primers: {len(specific_primers)}/{len(primers)}")
            print(f"  Rejected due to multiple exact matches: {rejection_reasons['multiple_exact']}")
            print(f"  Rejected due to unexpected position: {rejection_reasons['unexpected_position']}")
            print(f"  Rejected due to partial matches: {rejection_reasons['partial_matches']}")
            print(f"  Rejected due to 3' end off-targets: {rejection_reasons['three_prime_matches']}")
        
        return specific_primers


def test_specificity_checker():
    """Test function for SpecificityChecker"""
    from primer_generator import PrimerCandidate
    
    # Create test sequence with a repeated region
    unique_region = "ATCGATCGATCGATCG"
    repeated_region = "GATTACAGATTACA"
    
    template = unique_region * 10 + repeated_region + "NNNN" * 50 + repeated_region + unique_region * 10
    
    print(f"Template length: {len(template)} bp")
    print(f"Repeated region: {repeated_region}")
    
    checker = SpecificityChecker(template)
    
    # Test primer with unique sequence
    unique_primer = PrimerCandidate(
        sequence=unique_region + "ATCG",
        position=0,
        length=20,
        is_forward=True
    )
    
    # Test primer with repeated sequence  
    repeated_primer = PrimerCandidate(
        sequence=repeated_region,
        position=160,
        length=len(repeated_region),
        is_forward=True
    )
    
    print("\n--- Testing unique primer ---")
    is_specific, details = checker.is_specific(unique_primer)
    print(f"Is specific: {is_specific}")
    print(f"Reason: {details['reason']}")
    print(f"Exact matches: {len(details['exact_matches'])}")
    
    print("\n--- Testing repeated primer ---")
    is_specific, details = checker.is_specific(repeated_primer)
    print(f"Is specific: {is_specific}")
    print(f"Reason: {details['reason']}")
    print(f"Exact matches at positions: {details['exact_matches']}")
    
    # Test filtering
    primers = [unique_primer, repeated_primer]
    specific = checker.filter_specific_primers(primers)
    print(f"\n{len(specific)}/{len(primers)} primers passed specificity filter")


if __name__ == "__main__":
    test_specificity_checker()