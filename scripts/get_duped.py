import sys

with open(sys.argv[1], 'r') as fi:

	lines = fi.readlines()

	curr_header = ""

	for line in lines:

		if line.startswith('>'):
			curr_header = line
			header_prefix = curr_header.split(':')[0]
			header_suffix = curr_header.split(':')[-1]
			quant = int(header_prefix.split('.')[-1])

		else:

			for i in range(0, quant):
				print(header_prefix + "." + str(i) + ":" + header_suffix, end='')
				print(line, end='')


