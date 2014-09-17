#  Copyright (C) 2014
#      Jakub Krajniak (jkrajniak at gmail.com)
#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
#
#  This file is part of ESPResSo++.
#
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


import datetime
import os
import sys



class TerminalFileLogger(object):
  """Enables the pipeling messages both to the file and to the terminal.

  Args:
    root_dir: The optional root dir where the log file will be stored.
    log_name: The name of the log file.

  Usage:
    terminal_logger = TerminalFileLogger()
    terminal_logger.enable()

    print 'Test'  # This will be printed out on the console and saved in the run_log.log file

    terminal_logger.disable()
    terminal_logger.close()
  """

  def __init__(self, root_dir, log_name='run_log.log'):
    if log_name is None:
      log_name = '%s.log' % sys.argv[0]

    self.terminal = sys.stdout
    self.log_path = os.path.join(root_dir, log_name)
    self.log = open(self.log_path, 'w')
    self.log.write('LOG FILE date=%s\n' % datetime.datetime.now())
    self.log.write('PWD: %s\n' % os.environ.get('PWD'))
    self.log.write('USER: %s\n====================\n' % os.environ.get('USER'))

  def close(self):
    """Closes the logging file."""
    self.log.write('\n====================\nLOG file closed date=%s' % datetime.datetime.now())
    self.log.close()
    self.disable()

  def write(self, message):
    """Write message both to the logging file and send out to the console.

    Args:
      message: The message to print out.
    """
    self.terminal.write(message)
    self.log.write(message)

  def writelines(self, lines):
    """Write multiple messages.

    Args:
      lines: The array with messages.
    """
    for l in lines:
      self.terminal.write(l)
      self.log.write(l)

  def enable(self):
    """Enable logger."""
    sys.stdout = self

  def disable(self):
    """Disable logger."""
    sys.stdout = self.terminal

