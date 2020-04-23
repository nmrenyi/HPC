#encoding=utf-8

import os
import sys


execute = 'srun'
cnt = 0
p_n = [1, 4, 9, 16, 25]
m_n = [1080, 2160, 4320, 7200, 14400]
for p in p_n:
	for m in m_n:
		print "Running %s_%s" % (str(p), str(m))
		with open('tmp.txt', 'w') as f:
			f.write(str(m) + '\n')
		
		if 0 != os.system("%s -n %d parallel < tmp.txt > ./result/%s_%s.txt" % (execute, p, str(p), str(m))):
			sys.stdout.write("Run program1 failed\n")
		