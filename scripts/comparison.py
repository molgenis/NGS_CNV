import os
import argparse

# Import required classes
from classes.arraycnv import ArrayCnv
from classes.conifercall import ConiferCall
from classes.gatkcall import GatkCall

# Import comparison scripts
import comparison.comparison as comcom

# Import parameter scripts
import parameters.parameters as parpar

# Import util scripts
import utils.filereader as ufr
import utils.filewriter as ufw

#Make some parameter defining variables
TOOL_CHOICES = ["gatk4_conifer", "gatk4_exomedepth", "conifer_exomedepth"]
REQUIRED_PARAMS = {"": ["", ""]}
OPTIONAL_PARAMS = {}
PARAM_TYPES = {}
TOOL_USAGE = {}