import sys

output_handler = open(sys.argv[2], 'w')
size_map = dict()
for line in open(sys.argv[1]):
    line = line.strip().split("\t")
    if line[0] not in size_map:
        size_map[line[0]] = 0
    size_map[line[0]] += int(line[2])

for key in size_map:
    output_handler.write(key + "\t" + str(size_map[key]) + "\n")
