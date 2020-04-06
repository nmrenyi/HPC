#encoding=utf-8

import random
import sys

f = sys.stdout
data = ''
threads = random.randint(1, 20)
n = random.randint(0, 90)
data += str(threads) + '\n' + str(n) + '\n'
f.write(data)
