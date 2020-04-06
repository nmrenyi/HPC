#encoding=utf-8

import random
import sys


f = sys.stdout
data = ''
for i in range(40):
	data += str(random.randint(0, 10)) + ' '
f.write(data + '\n')
