#!/usr/bin/env python

import sys
import os
from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup

def configuration(parent_package='',top_path=None):
    config = Configuration('admom',parent_package,top_path)
    #config.add_extension('_admom', ['admom/ad_momi.f'])
    config.add_extension('_admomf', ['admom.f'])
    return config

setup(configuration=configuration)


