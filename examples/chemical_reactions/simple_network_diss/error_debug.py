# -*- coding: utf-8 -*-
#
# Author: Ian Eure <http://github.com/ieure>, <http://atomized.org>
#

"""Enter the debugger on exceptions.

example:

from __future__ import with_statement
from error_debug import debug

with debug():
    raise Exception("Just testing")
"""

from contextlib import contextmanager

@contextmanager
def debug(use_pdb=True):
    """When use_pdb is True, enter the debugger if an exception is raised."""
    try:
        yield
    except Exception, e:
        if not use_pdb:
            raise
        import sys
        import traceback
        import pdb
        info = sys.exc_info()
        traceback.print_exception(*info)
        pdb.post_mortem(info[2])
