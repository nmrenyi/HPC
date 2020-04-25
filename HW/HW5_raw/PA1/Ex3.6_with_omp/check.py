#encoding=utf-8

import os
import sys


execute = 'srun'
cnt = 0
p_n = [4, 9, 16]
m_n = [2160, 4320, 7200]
thread_count = 10
for p in p_n:
	for m in m_n:
		print "Running %s_%s" % (str(p), str(m))
		with open('tmp.txt', 'w') as f:
			f.write(str(m) + '\n')
		
		if 0 != os.system("%s -n %d parallel %d < tmp.txt > ./result/%s_%s.txt" % (execute, p, thread_count, str(p), str(m))):
			sys.stdout.write("Run program1 failed\n")
		