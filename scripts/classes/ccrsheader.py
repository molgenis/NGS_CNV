#!/usr/bin/env python
from collections import OrderedDict

class CcrsHeader:
    def __init__(self, ccrs_header):
        self.header_data = ccrs_header

    def add_header_entry(self):
        

    def field_in_header(self, headerfield):
        return headerfield in self.header_data

    def get_header_data(self, headerfield):
        if headerfield in self.header_data:
            return self.header_data[headerfield]
        return {}
