import numpy as np
import matplotlib.pyplot as plt

filename = "problem2.txt"
with open(filename, 'r') as infile:
    x = []
    u = []
    for line in infile:
        try:
            xi = float(line.split()[0])
            x.append(xi)
            ui = float(line.split()[1])
            u.append(ui)
        except:
            continue

fig = plt.figure(figsize=(5,4))
plt.plot(x,u)
plt.title("Analytical solution $u(x)$")
plt.xlabel("$x$"); plt.ylabel("$u$")
plt.savefig("problem2.pdf")
plt.show()