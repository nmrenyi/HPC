#encoding=utf-8

import random
import sys

f = sys.stdout
data = ''

TASKS_COUNT = 15

tasks_list = []
# append different tasks into the list according to the ratio
for i in range(int(TASKS_COUNT * 4 / 15)):
	tasks_list.append('0')
for i in range(int(TASKS_COUNT * 2 / 15)):
	tasks_list.append('1')
for i in range(int(TASKS_COUNT * 8 / 15)):
	tasks_list.append('2')
for i in range(int(TASKS_COUNT * 1 / 15)):
	tasks_list.append('3')

random.shuffle(tasks_list) # shuffule the task list
data += ' '.join(tasks_list) + '\n'

one_number = tasks_list.count('0')
one_cnt = 0
insert_list = random.sample(range(0, 2 * one_number), one_number)
insert_list = [str(x) for x in insert_list]

input_seq = ''
# give different number to the different task
for x in tasks_list:
	if x == '0':
		input_seq += insert_list[one_cnt] + '\n'
		one_cnt += 1
	elif x == '1':
		if random.random() > 0.3:
			input_seq += random.choice(insert_list) + '\n'
		else:
			input_seq += str(random.randint(2 * one_number, 3 * one_number)) + '\n'
	elif x == '2':
		if random.random() > 0.3:
			input_seq += random.choice(insert_list) + '\n'
		else:
			input_seq += str(random.randint(2 * one_number, 3 * one_number)) + '\n'

data += input_seq

f.write(data)
