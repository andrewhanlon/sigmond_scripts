#!/usr/bin/env python

import sys
import getpass
import argparse
import logging

import launcher

def main():

  # Get command line arguments
  parser = argparse.ArgumentParser(description="Run SigMonD Analysis")
  parser.add_argument(
      "-c", "--configs", type=str, required=False, metavar='CONFIG-FILE(S)', nargs='+',
      help="Run sigmond by using a configuration file specified by CONFIG-FILE")
  parser.add_argument(
      "-v", "--verbose", action="store_const", dest="loglevel", const=logging.INFO,
      default=logging.WARNING, help="Sets logging level to INFO")
  parser.add_argument(
      "-d", "--debug", action="store_const", dest="loglevel", const=logging.DEBUG,
      default=logging.WARNING, help="Sets logging level to DEBUG")

  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s: %(message)s', handlers=[ExitOnExceptionHandler()], level=args.loglevel)

  try:
    launcher.launch(args.configs)
  except KeyboardInterrupt:
    print("\n\nGoodbye, {}".format(getpass.getuser()))

  return


# Thanks be to https://stackoverflow.com/a/48201163/191474
class ExitOnExceptionHandler(logging.StreamHandler):

  def emit(self, record):
    super().emit(record)
    if record.levelno in (logging.ERROR, logging.CRITICAL):
      raise SystemExit(-1)


if __name__ == "__main__":
  main()
