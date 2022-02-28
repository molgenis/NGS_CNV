#!/usr/bin/env python
class FpOverlap:
    def __init__(self, fpregion):
        self.region = fpregion
        self.overlaps = {}

    def add_overlap(self, overlapregion, overlappercentage):
        self.overlaps[overlapregion] = overlappercentage

    def get_overlap(self, overlapregion):
        if overlapregion in self.overlaps:
            return [overlapregion, self.overlaps[overlapregion]]
        return []

    def has_region(self, overlapregion):
        return overlapregion in self.overlaps
