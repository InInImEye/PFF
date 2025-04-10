import numpy as np

target_string = "(1,20) 22\n"
print(target_string)

results = target_string.replace("("," ").replace(")"," ").replace(","," ")

print(results.split())
