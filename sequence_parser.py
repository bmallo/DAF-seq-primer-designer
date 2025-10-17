"""
Sequence Parser Module
Handles reading FASTA files and extracting primer search regions
"""

from typing import Tuple, Optional, List
from Bio import SeqIO
from Bio.Seq import Seq
from config import TargetRegion


class SequenceParser:
    """Parse and prepare sequences for primer design"""
    
    def __init__(self, fasta_file: str, verbose: bool = True):
        """
        Initialize with a FASTA file.
        
        Args:
            fasta_file: Path to FASTA file containing the template sequence
            verbose: Whether to print progress messages
        """
        self.fasta_file = fasta_file
        self.sequence = None
        self.sequence_id = None
        self.sequence_length = None
        self.verbose = verbose
        self._load_sequence()
    
    def _load_sequence(self):
        """Load the sequence from FASTA file"""
        try:
            # Read the first sequence from the FASTA file
            record = next(SeqIO.parse(self.fasta_file, "fasta"))
            self.sequence = str(record.seq).upper()
            self.sequence_id = record.id
            self.sequence_length = len(self.sequence)
            if self.verbose:
                print(f"Loaded sequence '{self.sequence_id}' ({self.sequence_length} bp)")
        except StopIteration:
            raise ValueError(f"No sequences found in {self.fasta_file}")
        except Exception as e:
            raise ValueError(f"Error reading FASTA file: {e}")
    
    def find_sequence_positions(self, query_seq: str, allow_multiple: bool = False) -> List[int]:
        """
        Find all positions where a query sequence appears in the template.
        Searches both forward strand and reverse complement.
        
        Args:
            query_seq: The sequence to search for
            allow_multiple: If False, raises error if multiple matches found
            
        Returns:
            List of start positions (0-based) where sequence is found
            
        Raises:
            ValueError: If sequence not found or multiple found when not allowed
        """
        query_seq = query_seq.upper().strip()
        positions = []
        
        # Search forward strand
        start = 0
        while True:
            pos = self.sequence.find(query_seq, start)
            if pos == -1:
                break
            positions.append(pos)
            start = pos + 1
        
        # Search reverse complement
        query_rc = str(Seq(query_seq).reverse_complement())
        start = 0
        while True:
            pos = self.sequence.find(query_rc, start)
            if pos == -1:
                break
            # Only add if not already found on forward strand
            if pos not in positions:
                positions.append(pos)
            start = pos + 1
        
        # Validate results
        if len(positions) == 0:
            raise ValueError(f"Sequence not found in template: {query_seq[:50]}...")
        
        if len(positions) > 1 and not allow_multiple:
            raise ValueError(f"Sequence found at multiple positions ({len(positions)} matches): {query_seq[:50]}...\n"
                           f"Positions: {positions}\n"
                           f"Please provide coordinates instead or use a more specific sequence.")
        
        return sorted(positions)
    
    def resolve_target_region(self, target_region: TargetRegion) -> TargetRegion:
        """
        Resolve any sequence-based regions to coordinates.
        Updates the TargetRegion object with actual coordinates.
        
        Args:
            target_region: TargetRegion that may contain sequences instead of coordinates
            
        Returns:
            Updated TargetRegion with all coordinates resolved
        """
        # Resolve core region
        if target_region.core_sequence is not None:
            positions = self.find_sequence_positions(target_region.core_sequence)
            target_region.core_start = positions[0]
            target_region.core_end = positions[0] + len(target_region.core_sequence)
            print(f"Core region '{target_region.core_sequence[:30]}...' found at position {target_region.core_start}-{target_region.core_end}")
        
        # Resolve forward primer region
        if target_region.forward_region_sequence is not None:
            positions = self.find_sequence_positions(target_region.forward_region_sequence, allow_multiple=True)
            # Find the position that makes most sense (closest to but before core region)
            best_pos = None
            for pos in positions:
                seq_end = pos + len(target_region.forward_region_sequence)
                if seq_end <= target_region.core_start:
                    if best_pos is None or pos > best_pos:
                        best_pos = pos
            
            if best_pos is None:
                raise ValueError(f"Forward region sequence found, but not in valid position (must be upstream of core region)")
            
            target_region.forward_region = (best_pos, best_pos + len(target_region.forward_region_sequence))
            print(f"Forward primer region found at position {target_region.forward_region[0]}-{target_region.forward_region[1]}")
        
        # Resolve reverse primer region
        if target_region.reverse_region_sequence is not None:
            positions = self.find_sequence_positions(target_region.reverse_region_sequence, allow_multiple=True)
            # Find the position that makes most sense (closest to but after core region)
            best_pos = None
            for pos in positions:
                if pos >= target_region.core_end:
                    if best_pos is None or pos < best_pos:
                        best_pos = pos
            
            if best_pos is None:
                raise ValueError(f"Reverse region sequence found, but not in valid position (must be downstream of core region)")
            
            target_region.reverse_region = (best_pos, best_pos + len(target_region.reverse_region_sequence))
            print(f"Reverse primer region found at position {target_region.reverse_region[0]}-{target_region.reverse_region[1]}")
        
        return target_region
    
    def get_primer_search_regions(self, target_region: TargetRegion) -> Tuple[Tuple[int, int], Tuple[int, int]]:
        """
        Determine the regions where forward and reverse primers should be searched.
        
        Args:
            target_region: TargetRegion object defining the amplicon
            
        Returns:
            Tuple of (forward_region, reverse_region) where each is (start, end) in 0-based coords
        """
        core_start = target_region.core_start
        core_end = target_region.core_end
        min_amplicon = target_region.min_amplicon
        max_amplicon = target_region.max_amplicon
        
        # If specific regions are provided, use those
        if target_region.forward_region and target_region.reverse_region:
            fwd_region = target_region.forward_region
            rev_region = target_region.reverse_region
            
            # Validate that these regions make sense
            if fwd_region[1] > core_start:
                raise ValueError("Forward primer region overlaps with core region")
            if rev_region[0] < core_end:
                raise ValueError("Reverse primer region overlaps with core region")
                
            return fwd_region, rev_region
        
        # Otherwise, calculate based on amplicon size constraints
        # Forward primer region: from (core_start - max_amplicon) to (core_start - min_amplicon + core_size)
        # This ensures the amplicon will be between min and max size
        
        core_size = core_end - core_start
        
        # Forward primer search region
        fwd_start = max(0, core_start - max_amplicon + core_size)
        fwd_end = max(0, core_start - min_amplicon + core_size)
        
        # Make sure we have some space for forward primers
        if fwd_end <= fwd_start:
            fwd_end = core_start
            fwd_start = max(0, fwd_end - 500)  # Default to 500bp region
        
        # Reverse primer search region  
        rev_start = min(self.sequence_length, core_end + min_amplicon - core_size)
        rev_end = min(self.sequence_length, core_end + max_amplicon - core_size)
        
        # Make sure we have some space for reverse primers
        if rev_end <= rev_start:
            rev_start = core_end
            rev_end = min(self.sequence_length, rev_start + 500)  # Default to 500bp region
        
        return (fwd_start, fwd_end), (rev_start, rev_end)
    
    def get_forward_search_sequence(self, region: Tuple[int, int]) -> str:
        """
        Extract the sequence for forward primer search.
        
        Args:
            region: (start, end) in 0-based coordinates
            
        Returns:
            Sequence string for forward primer search
        """
        start, end = region
        if start < 0 or end > self.sequence_length:
            raise ValueError(f"Region ({start}, {end}) is out of bounds for sequence length {self.sequence_length}")
        return self.sequence[start:end]
    
    def get_reverse_search_sequence(self, region: Tuple[int, int]) -> str:
        """
        Extract the sequence for reverse primer search.
        Reverse primers bind to the reverse complement, so we need to return
        the reverse complement of the region.
        
        Args:
            region: (start, end) in 0-based coordinates on the forward strand
            
        Returns:
            Reverse complement sequence for reverse primer search
        """
        start, end = region
        if start < 0 or end > self.sequence_length:
            raise ValueError(f"Region ({start}, {end}) is out of bounds for sequence length {self.sequence_length}")
        
        forward_seq = self.sequence[start:end]
        # Get reverse complement
        seq_obj = Seq(forward_seq)
        reverse_comp = str(seq_obj.reverse_complement())
        return reverse_comp
    
    def get_amplicon_sequence(self, fwd_primer_pos: int, rev_primer_pos: int, 
                             fwd_primer_length: int, rev_primer_length: int) -> str:
        """
        Get the amplicon sequence given primer positions.
        
        Args:
            fwd_primer_pos: Start position of forward primer (0-based)
            rev_primer_pos: Start position of reverse primer (0-based, on forward strand)
            fwd_primer_length: Length of forward primer
            rev_primer_length: Length of reverse primer
            
        Returns:
            Amplicon sequence from start of forward primer to end of reverse primer
        """
        amplicon_start = fwd_primer_pos
        amplicon_end = rev_primer_pos + rev_primer_length
        return self.sequence[amplicon_start:amplicon_end]
    
    def get_amplicon_size(self, fwd_primer_pos: int, rev_primer_pos: int,
                         fwd_primer_length: int, rev_primer_length: int) -> int:
        """
        Calculate amplicon size.
        
        Args:
            fwd_primer_pos: Start position of forward primer (0-based)
            rev_primer_pos: Start position of reverse primer (0-based, on forward strand)  
            fwd_primer_length: Length of forward primer
            rev_primer_length: Length of reverse primer
            
        Returns:
            Amplicon size in bp
        """
        return (rev_primer_pos + rev_primer_length) - fwd_primer_pos
    
    def validate_amplicon_covers_core(self, fwd_primer_pos: int, rev_primer_pos: int,
                                     fwd_primer_length: int, rev_primer_length: int,
                                     target_region: TargetRegion) -> bool:
        """
        Check if the amplicon covers the entire core region.
        
        Args:
            fwd_primer_pos: Start position of forward primer (0-based)
            rev_primer_pos: Start position of reverse primer (0-based, on forward strand)
            fwd_primer_length: Length of forward primer
            rev_primer_length: Length of reverse primer
            target_region: TargetRegion object
            
        Returns:
            True if amplicon covers the core region, False otherwise
        """
        amplicon_start = fwd_primer_pos + fwd_primer_length
        amplicon_end = rev_primer_pos
        
        core_start = target_region.core_start
        core_end = target_region.core_end
        
        return amplicon_start <= core_start and amplicon_end >= core_end


def test_sequence_parser():
    """Test function for SequenceParser"""
    # This is just for testing - you'd provide your actual FASTA file
    from config import TargetRegion
    
    # Create a test FASTA file with a recognizable sequence
    test_fasta = "test_sequence.fasta"
    # Create a sequence with a specific core region we can find
    core_seq = "GATTACAGATTACA"
    before_core = "ATCGATCG" * 250  # 2000 bp
    after_core = "GCTAGCTA" * 250   # 2000 bp
    test_seq = before_core + core_seq + after_core
    
    with open(test_fasta, 'w') as f:
        f.write(">test_sequence\n")
        f.write(test_seq)
    
    # Test the parser
    parser = SequenceParser(test_fasta)
    print(f"Sequence length: {parser.sequence_length}")
    
    # Test with coordinates
    print("\n--- Testing with coordinates ---")
    target1 = TargetRegion(
        core_start=2000,
        core_end=2014,
        min_amplicon=2000,
        max_amplicon=5000
    )
    
    fwd_region, rev_region = parser.get_primer_search_regions(target1)
    print(f"Forward primer region: {fwd_region}")
    print(f"Reverse primer region: {rev_region}")
    
    # Test with sequence search
    print("\n--- Testing with sequence search ---")
    target2 = TargetRegion(
        core_sequence=core_seq,
        min_amplicon=2000,
        max_amplicon=5000
    )
    
    target2 = parser.resolve_target_region(target2)
    fwd_region2, rev_region2 = parser.get_primer_search_regions(target2)
    print(f"Forward primer region: {fwd_region2}")
    print(f"Reverse primer region: {rev_region2}")
    
    # Test sequence extraction
    fwd_seq = parser.get_forward_search_sequence(fwd_region2)
    rev_seq = parser.get_reverse_search_sequence(rev_region2)
    print(f"Forward search sequence length: {len(fwd_seq)}")
    print(f"Reverse search sequence length: {len(rev_seq)}")


if __name__ == "__main__":
    test_sequence_parser()