import sys

seq_header = ""
N_count = 0

output = open(sys.argv[2], "w")
for line in open(sys.argv[1]):
    if line.startswith(">"):
        if seq_header != "":
            output.write(seq_header + "\t" + str(N_count) + "\n")
        N_count = 0
        seq_header = line[1:].strip()
    else:
        N_count += line.upper().count("N")

output.write(seq_header + "\t" + str(N_count) + "\n")