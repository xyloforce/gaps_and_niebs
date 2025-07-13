import sys
import statistics

microsat_file = open(sys.argv[2], "w")

def get_gc(microsat):
  microsat = microsat.upper()
  c = microsat.count("C")
  g = microsat.count("G")
  len_microsat = len(microsat)
  return (c + g) / len_microsat
    
print("getting gc values")
gc_values = [get_gc(line.split("\t")[3])
             for line in open(sys.argv[1])
             if not line.startswith("#")]
quantiles = statistics.quantiles(gc_values)

print("annotating lines")
for line in open(sys.argv[1]):
  if not line.startswith("#"):
    line = line.strip().split("\t")
    gc_level = ""
    gc_value = get_gc(line[3])
    if gc_value <= quantiles[0]:
      gc_level = "Q1"
    elif gc_value <= quantiles[1]:
      gc_level = "Q2"
    elif gc_value <= quantiles[2]:
      gc_level = "Q3"
    else:
      gc_level = "Q4"
    
    line[3] = gc_level
    microsat_file.write("\t".join(line) + "\n")