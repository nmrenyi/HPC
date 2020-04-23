import os
import sys
cnt  =0
while True:
	cnt += 1
	sys.stdout.write("Running Case #%s ... " % cnt)
	# if 0 != os.system("%s > input.txt" % sys.argv[2]):
	# 	sys.stdout.write("Run datamaker failed\n")
	# 	break
	# print sys.argv[1]
	# print sys.argv[2]

	print "%s  < input.txt > output2.txt" % sys.argv[1]
	if 0 != os.system("%s  < input.txt " % sys.argv[1]):
		sys.stdout.write("Run program2 failed\n")
		break
