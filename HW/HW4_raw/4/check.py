#encoding=utf-8

import os
import sys

def checkAns(outputFile1, outputFile2):
	with open(outputFile1, "r") as f:
		output1 = f.readlines()
	with open(outputFile2, "r") as f:
		output2 = f.readlines()
	return "".join(output1) == "".join(output2)


cnt = 0
while True:
	cnt += 1
	sys.stdout.write("Running Case #%s ... " % cnt)
	if 0 != os.system("%s > input.txt" % sys.argv[3]):
		sys.stdout.write("Run datamaker failed\n")
		break
	
	if 0 != os.system("%s < input.txt > output1.txt" % sys.argv[1]):
		sys.stdout.write("Run program1 failed\n")
		break
	# print 'serial complete'
	if 0 != os.system("%s  < input.txt > output2.txt" % sys.argv[2]):
		sys.stdout.write("Run program2 failed\n")
		break
	# print 'parallel complete'
	# if 0 != os.system("%s  < input.txt" % sys.argv[2]):
	# 	sys.stdout.write("Run program2 failed\n")
	# 	break

	if not checkAns("output1.txt", "output2.txt"):
		sys.stdout.write("Wrong\n")
		break
	sys.stdout.write("OK\n")

	# if cnt == 100:
	# 	break