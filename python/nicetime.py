#!/usr/bin/env python3

import sys

#############################################################################
def NiceTime(s):
  """ Express time in human readable format. """
  if s<60: return f"{s}s"
  m,s = divmod(s, 60)
  if m<60: return f"{m}m:{s:02d}s"
  h,m = divmod(m, 60)
  if h<24: return f"{h}h:{m:02d}m:{s:02d}s"
  d,h = divmod(h, 24)
  return f"{d}d:{h:02d}h:{m:02d}m:{s:02d}s"

#############################################################################
if __name__=="__main__":
  if len(sys.argv)!=2:
    print(f"ERROR: syntax {sys.argv[0]} <SECONDS>")
    sys.exit(1)

  s = int(sys.argv[1])
  print(NiceTime(s))

