"""
FYS4150 project 1
problem 7, 8
plots approximations v for different number of steps
plots also relative and absolute errors compared to exact solution
"""

import numpy as np
import matplotlib.pyplot as plt


n_list = ["10","100","1000"]#,"10000"]
problem = input("problem 7 or 8?\n")

if problem == "7":
    fig1, ax1 = plt.subplots(figsize=(5,4))

if problem == "8":
    fig1, ax1 = plt.subplots(figsize=(5,4))
    fig2, ax2 = plt.subplots(figsize=(5,4))

for n in n_list:
    filename = "problem" + problem + "_" + n + ".txt"
    with open(filename, 'r') as infile:
        x = []
        u = []
        delta = []
        eps = []
        for line in infile:
            try:
                if problem=="7":
                    xi = float(line.split()[0])
                    x.append(xi)
                    ui = float(line.split()[1])
                    u.append(ui)
                
                if problem=="8" and len(line.split())==4:
                    xi = float(line.split()[0])
                    x.append(xi)
                    ui = float(line.split()[1])
                    u.append(ui)
                    deltai = float(line.split()[2])
                    delta.append(deltai)
                    epsi = float(line.split()[3])
                    eps.append(epsi)

            except:
                continue
    if problem=="7":
        ax1.plot(x,u, label=n)
    if problem=="8":
        ax1.plot(x, delta, label=n)
        ax2.plot(x, eps, label=n)

if problem=="7":
    x = np.linspace(x[0], x[-1], int(n))
    u = 1 - (1-np.exp(-10))*x - np.exp(-10*x)
    ax1.plot(x, u, label="exact", linestyle="--")
    fig1.suptitle("Approximations for different $n_{steps}$")
    ax1.set_xlabel("$x$"); ax1.set_ylabel("$u$")
    ax1.legend()
    fig1.savefig("problem7.pdf")

if problem=="8":
    fig1.suptitle("Absolute error")
    ax1.set_xlabel("$x$"); ax1.set_ylabel("|u - v|")
    ax1.set_yscale("log")
    ax1.legend()
    fig1.savefig("problem8a.pdf")
    
    fig2.suptitle("Relative error")
    ax2.set_xlabel("$x$"); ax2.set_ylabel("|1 - v/u|")
    ax2.set_yscale("log")
    ax2.legend()
    fig2.savefig("problem8b.pdf")

plt.show()