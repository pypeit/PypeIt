"""Tests for fiber block-slit extraction logic."""
import numpy as np
import pytest
from pypeit.spectrographs.mmt_binospec import MMTBINOSPECIFUSpectrograph


class TestFiberBlockConfig:
    """Tests for fiber block configuration methods."""

    @classmethod
    def setup_class(cls):
        cls.spec = MMTBINOSPECIFUSpectrograph()

    def test_get_fiber_blocks_det01(self):
        """DET01 should have 21 blocks: 5 sky (8 fibers) + 16 science (20)."""
        blocks = self.spec.get_fiber_blocks(1)
        assert len(blocks) == 21
        # Sky blocks: 1, 6, 11, 16, 21
        for b in [0, 5, 10, 15, 20]:  # 0-indexed
            assert blocks[b]['type'] == 'sky'
            assert blocks[b]['nfibers'] == 8
        # Science blocks
        for b in [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14, 16, 17, 18, 19]:
            assert blocks[b]['type'] == 'science'
            assert blocks[b]['nfibers'] == 20

    def test_get_fiber_blocks_det02(self):
        """DET02 should have 21 live blocks (dead fibers in block -1 are excluded).

        DET02 has 4 dead fibers assigned to block_id -1 in the reference profile.
        These are excluded from the returned block list, leaving 21 standard blocks.
        Some science blocks have 19 fibers instead of 20 because of those dead fibers.
        """
        blocks = self.spec.get_fiber_blocks(2)
        assert len(blocks) == 21
        # Sky blocks are at the same positions as DET01
        for b in [0, 5, 10, 15, 20]:  # 0-indexed
            assert blocks[b]['type'] == 'sky'
            assert blocks[b]['nfibers'] == 8
        # Science blocks have 19 or 20 fibers (some have dead fibers replaced)
        for b in [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14, 16, 17, 18, 19]:
            assert blocks[b]['type'] == 'science'
            assert blocks[b]['nfibers'] in (19, 20)

    def test_block_fiber_positions(self):
        """Fiber positions within blocks should be monotonically ordered."""
        blocks = self.spec.get_fiber_blocks(1)
        for block in blocks:
            positions = block['fiber_positions']
            assert np.all(np.diff(positions) > 0), \
                f"Block {block['block_id']} positions not monotonic"

    def test_block_gaps(self):
        """Inter-block gaps should be > 50 pixels."""
        blocks = self.spec.get_fiber_blocks(1)
        for i in range(len(blocks) - 1):
            gap = blocks[i+1]['min_pix'] - blocks[i]['max_pix']
            assert gap > 50, f"Gap between blocks {i} and {i+1} is only {gap:.1f} px"


class TestFiberMatchingInBlocks:
    """Tests for identifying fibers within block-slits."""

    @classmethod
    def setup_class(cls):
        cls.spec = MMTBINOSPECIFUSpectrograph()

    def test_identify_fibers_in_block(self):
        """Given fiber peak positions within a block, identify each fiber."""
        blocks = self.spec.get_fiber_blocks(1)
        # Use the reference positions as "detected" positions for block 2 (science, 20 fibers)
        block = blocks[1]  # block 2 (0-indexed = 1)
        detected_positions = block['fiber_positions']
        result = self.spec.identify_fibers_in_block(
            det=1, block_idx=1, detected_positions=detected_positions)
        assert len(result['fiber_id']) == 20
        assert len(result['fiber_name']) == 20
        # All should match since we used reference positions
        assert np.all(result['fiber_id'] > 0)

    def test_identify_fibers_with_offset(self):
        """Fibers should still match with a small position offset."""
        blocks = self.spec.get_fiber_blocks(1)
        block = blocks[1]
        # Add 2-pixel offset (simulating flexure)
        detected_positions = block['fiber_positions'] + 2.0
        result = self.spec.identify_fibers_in_block(
            det=1, block_idx=1, detected_positions=detected_positions)
        assert np.all(result['fiber_id'] > 0)
