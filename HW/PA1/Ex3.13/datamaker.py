#encoding=utf-8

import random
import sys

def dataMaker(length_range, data_range):
	f = sys.stdout
	order = random.randint(1, length_range)
	f.write(str(order) + "\n")
	v1 = ''
	v2 = ''
	for i in range(order):
		v1 += str(random.random() * data_range) + ' '
		v2 += str(random.random() * data_range) + ' '
	f.write(v1 + '\n')
	f.write(v2 + '\n')



dataMaker(200000, 100)
