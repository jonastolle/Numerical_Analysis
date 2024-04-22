# Numerical Analysis

**Jonas Tölle**

*Lecture notes for MS-C1650 Numerical Analysis, Aalto University*

*Last updated: 22.4.2024*

## Floating-point numbers

We refer to
-  [\[Fundamentals of Numerical Computation, Floating-point numbers\]](https://fncbook.github.io/fnc/intro/floating-point.html)
- ...

> The set of real numbers $\mathbb{R}$ is infinite in two ways: it is unbounded and continuous. In most practical computing, the second kind of infiniteness is more consequential than the first kind, so we turn our attention there first.

Instead of $\mathbb{R}$, we shall introduce the set of (binary) *floating-point numbers* $\mathbb{F}$ whose members are zero and all numbers of the form

$$\pm(1+f)\cdot 2^e$$,

where $e$ is an integer called the *exponent*, and $1+f$ is called the *mantissa* or *significand*, in which,

$$f=\sum_{i=1}^d b_i 2^{-i},\quad b_i\in \{0,1\},$$

for a fixed integer $d$, which we call *precision*.

A general floating point number $x$ has the representation
$$x=\pm (d_0.d_1 d_2 \ldots d_9)_k \cdot k^e$$

<!--stackedit_data:
eyJoaXN0b3J5IjpbLTE1MDExNzkyNzUsLTE4NjcxNzYxNzVdfQ
==
-->