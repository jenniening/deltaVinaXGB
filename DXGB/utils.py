import os
import sys
import errno
import distutils.spawn
import DXGB.path_const as pc


def get_tool(name):
    executable = distutils.spawn.find_executable(name)
    if executable:
        return executable
    else:
        try:
            tool = eval('pc.' + name.split('.')[0].upper())
        except:
            tool = None
            raise ValueError('Cannot find ' + name)
    return tool 