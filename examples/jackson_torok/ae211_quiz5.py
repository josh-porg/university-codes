"""Jackson Torok, AE 211 MATLAB quiz 5 (``Matlab_Quiz5_Jackson_Torok``) and a scratch plot (``Untitled4``).

The quiz is a conversation driven by ``input``: name, favourite colour (a
colour starting with "g" is the program's favourite too) and sibling count.
The answers are options here. ``Untitled4`` plotted ``X.^coth(Y)`` on a
grid, which is complex for negative X; its real part is plotted.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

p = argparse.ArgumentParser()
p.add_argument("--name", default="Jackson")
p.add_argument("--color", default="green")
p.add_argument("--siblings", type=int, default=2)
a = p.parse_args()

print(f"Hi, {a.name}, it is nice to meet you! What is your favorite color?")
print(f"I really like {a.color}" + (", and it is actually my favorite!" if a.color[:1] == "g" else ", but my favorite is green!"))
print(f"Wow! {a.siblings} is {'a lot!' if a.siblings >= 2 else 'not that many.'} I have 6 siblings!")

X, Y = np.meshgrid(np.arange(-1, 1.0001, 0.025), np.arange(-1, 1.0001, 0.025))
with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
    z = np.power(X.astype(complex), 1 / np.tanh(Y))
plt.plot(z.real)
plt.axis([-10, 10, -10, 10])
plt.show()
