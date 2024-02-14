for line in open("sample_lists"):
  if "include" not in line: continue
  line = line.strip().replace('#include "','').replace('"','') 
  
  print(line)

