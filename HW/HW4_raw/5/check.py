#encoding=utf-8

import os
import sys

# if len(sys.argv) != 3:
# 	sys.stdout.write("Usage: %s <program1> <program2> <datamaker>\n" % sys.argv[0])
# 	sys.stdout.write("    Example on Linux: python %s ./a \"python b.py\" \"python3 datamaker.py\"\n" % sys.argv[0])
# 	exit(1)

cnt = 0
while True:
		cnt += 1
		sys.stdout.write("Running Case #%s ... " % cnt)
		if 0 != os.system("%s > input.txt" % sys.argv[2]):
			sys.stdout.write("Interrupted by keyboard or Run datamaker failed\n")
			break
		if 0 != os.system("%s < input.txt" % sys.argv[1]):
			sys.stdout.write("Interrupted by keyboard or Wrong Result!\n")
			break
		sys.stdout.write("OK\n")
