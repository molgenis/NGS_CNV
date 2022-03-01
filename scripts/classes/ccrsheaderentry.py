#!/usr/bin/env python
class CcrsHeaderEntry:
    def __init__(self):
        self.entry_data = {}

    def get_entry(self, entrykey):
        if entrykey in self.entry_data:
            return self.entry_data[entrykey]
        return None
