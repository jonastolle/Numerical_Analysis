# Numerical Analysis

**Jonas Tölle**

*Lecture notes for MS-C1650 Numerical Analysis, Aalto University*

*Last updated: 28.4.2024*

## Floating-point numbers

We refer to
-  [\[Tobin A. Driscoll and Richard J. Braun, Fundamentals of Numerical Computation, Floating-point numbers\]](https://fncbook.github.io/fnc/intro/floating-point.html)
- ...

> The set of real numbers $\mathbb{R}$ is infinite in two ways: it is unbounded and continuous. In most practical computing, the second kind of infiniteness is more consequential than the first kind, so we turn our attention there first.

Instead of $\mathbb{R}$, we shall introduce the set of (binary) *floating-point numbers* $\mathbb{F}$.

A general floating point number $x$ has the representation
$$x=\pm (d_0.d_1 d_2 \ldots d_p)_k \cdot k^e$$

Binary floating point numbers are zero and all numbers of the form
$$\pm(1+f)\cdot 2^e$$,
where $e$ is an integer called the *exponent*, and $1+f$ is called the *mantissa* or *significand*, in which,
$$f=\sum_{i=1}^d b_i 2^{-i},\quad b_i\in \{0,1\},$$
for a fixed integer $d$, which we call *precision*. Here, $f\in [1,2)$.

## Condition number and stability

### Numerical differentiation

Recall Taylor's theorem, for a twice differentiable function $f:\mathbb{R}\to\mathbb{R}$,
$$f(x)=f(x_0)+f'(x_0)(x-x_0)+\frac{1}{2}(x-x_0)^2 f''(\xi),$$
for any $x_0,x\in\mathbb{R}$, where $\xi\in [x_0,x]$.

>What is that ξ (= "xi")? Under certain assumptions elementary functions have their series expansions. If the series is truncated, we have the Taylor polynomial. However, the residual has an explicit expression but due to application of an intermediate value theorem, the exact location of the point ξ is not known, i.e. "generic".

By setting $x:=z+h$, $x_0:=z$, we obtain the useful equivalent formula
$$f(z+h)=f(z)+f'(z)h+\frac{1}{2} f''(\xi)h^2,$$
for every $z,h\in\mathbb{R}$, $\xi\in [z,z+h]$.

## Convergence rate

Let $(x_k)_{n\in\mathbb{N}}$ be an infinite sequence of real numbers.  Let $s_k:=\sup_{l\ge k}x_l$, $k\in\mathbb{N}$, be the *supremum* (i.e., the lowest upper bound) of the *tail* of $(x_k)$. Define the $\limsup$ (*limes superior*) as
$$\limsup_{k\to\infty}x_k:=\lim_{k\to\infty}s_k\in[-\infty,+\infty].$$
Other than a limit, it always exists, but can be $\pm\infty$, if $(x_k)$ is bounded, the $\limsup$ is the largest limit of a converging subsequence.

>The *rate of convergence* can be used interchangeably with the *order of convergence*. However, there is some caution necessary, as different authors use different terminology here. Usually, the order of convergence always refers to the same thing, namely, the α-exponent in the denominator of the limit defining the order of convergence. Most confusingly, some authors call the order of convergence "rate of convergence", as e.g. [here](https://www.math-cs.gordon.edu/courses/ma342/handouts/rate.pdf). The English [Wikipedia article](https://en.wikipedia.org/wiki/Rate_of_convergence) calls it the order of convergence, whereas here the rate of convergence is the constant in the definition, which also determines the speed of convergence, together with the order of convergence. So, please always check the context, as the use of the terminology should be clear from it. If there is no definition, try to figure out what is meant in each text. As a rule of thumb: The "order of convergence" is a unique terminology in numerical analysis. The "rate of convergence" can mean at least two different things. I will use both words for the same thing, but will try to make clear what I mean from case to case. In any case, to be sure, use "order of convergence". My PhD advisor usually said that in mathematics "it's all hollow words" (meaning that one should check the definition).

## Finding roots of functions

Let $f:\mathbb{R}\to\mathbb{R}$ be continuous. We are interested in methods for finding *zeros*, that is, *roots* of $f$, in other words, $x\in\mathbb{R}$, such that $f(x)=0$.

### Bisection method

### Secant method

### Newton's method
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTI5NzQwMjY0NiwtNDE0NjE2MDIwLDEwOD
E1NjAyNjUsLTE1NTkyMzgzODIsLTEwMzAyOTkzNTksLTE4NzAx
OTU2MTMsLTE1MDExNzkyNzUsLTE4NjcxNzYxNzVdfQ==
-->