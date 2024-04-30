# Numerical Analysis

**Jonas Tölle**[^1]

[^1]: [jonas.tolle@aalto.fi](mailto:jonas.tolle@aalto.fi)

*Lecture notes for MS-C1650 Numerical Analysis, Aalto University*

*Last updated: 30.4.2024*

> **Intended learning outcomes.** After the course, the student will be able to...
> - explain the fundamental concepts of numerical analysis, like condition number, stablilty, and convergence rate;
> - construct the floating point numbers;
> - discuss and employ basic numerical algorithms like Newton's method;
> - use the Monte-Carlo method in basic problems in analysis and geometry;
> - apply different methods of interpolation polynomials and numerical quadrature rules;
> - understand the Euler scheme and linear multi-step methods for solving ordinary differential equations.

## Floating-point numbers
> The set of real numbers $\mathbb{R}$ is infinite in two ways: it is unbounded and continuous. In most practical computing, the second kind of infiniteness is more consequential than the first kind, so we turn our attention there first.

Instead of $\mathbb{R}$, we shall introduce the set of *floating-point numbers* (*floats*) $\mathbb{F}$. They come with different bases, precisions and exponent ranges, and other features. The basic representation is
$$x=\pm (d_0. d_1 d_2 \ldots d_p)_k\cdot k^e.$$
$k\in\mathbb{N}\setminus\{1\}$ is called *base* or *radix*, $p\in\mathbb{N}_0:=\mathbb{N}\cup\{0\}$ is called the precision, $d_i$, $i\in\{0,\ldots,p\}$, and the sequence of numbers
$$(d_0. d_1 d_2 \ldots d_p)_k:=\sum_{i=0}^p d_i k^{-i}$$
is called *mantissa* or *significand*. The exponent $e\in\mathbb{Z}$ is bounded $m\le e\le M$, where $m,M\in\mathbb{Z}$.
If $k=10$, we can read the usual decimal commas from the mantissa:
$$(1.01)_{10}=1\cdot 10^0+0\cdot 10^{-1}+1\cdot 10^{-2}=1.01.$$
If $k=2$, we have binary floats. In the binary case, we observe that we can always choose $d_0=1$, hence saving one bit, which can be expressed by $e$. We refer to this as *normalization*. The mantissa is always contained in the interval $[1,2)$.

>**Example.** (Toy floating point system). Binary floats of the type
>$$(1.b_1b_2)_2$$ with exponents $e=-1,0,1$.
> Hence $(1.00)_2=1$, $(1.01)_2=\frac{5}{4}$, $(1.10)_2=\frac{3}{2}$, and $(1.11)_2=\frac{7}{4}$. By multiplying with the exponents $2^{-1}=\frac{1}{2}$, $2^0=1$, $2^1=2$, we get the whole set:
> $e$ |    |      |    |
> -|-|-|-|-|-
> $1$ | $\frac{5}{4}$ | $\frac{3}{2}$  | $\frac{7}{4}$
> $2$ | $\frac{5}{2}$ | $3$ | $\frac{7}{2}$
> $\frac{1}{2}$ | $\frac{5}{8}$ | $\frac{3}{4}$ | $\frac{7}{8}$
> Important quantity, $(1.01)_2-1=\frac{1}{4}$, the so-called *machine epsilon*.

Define the machine epsilon by $\varepsilon:=2^{-p}=(1.00\ldots 01)_2-1$.

### Rounding

For the rounding function $\text{round}:\mathbb{R}\to\mathbb{F}$, we have 5 alternative definitions:
- Rounding to nearest (*default*) 
- Rounding to $+\infty$
- Rounding to $-\infty$
- Rounding to $0$
- Rounding away from $0$

Here's an easy-to-follow [video](https://youtu.be/p8u_k2LIZyo?si=Gfi9TId6x6BAcmpo) explaining floating point numbers (and a specific version of Newton's algorithm).

## Condition number and stability

### Stability of an algorithm

**Definition.** An algorithm or numerical process is called *stable* if small changes in the input produce small changes in the output. It is called *unstable* if large changes in the output are produced by small changes in the input.

An algorithm is stable, if every step is well-conditioned (i.e. has a uniformly bounded condition number). It is unstable if any step is ill-conditioned (i.e. the condition number may become arbitrarily large).

>**Example.** Consider evaluating $f(x)=\sqrt{1+x}-1$ for $x$ close to zero. The relative condition number 

### Numerical differentiation

Recall Taylor's theorem, for a twice differentiable function $f:\mathbb{R}\to\mathbb{R}$,
$$f(x)=f(x_0)+f'(x_0)(x-x_0)+\frac{1}{2}(x-x_0)^2 f''(\xi),$$
for any $x_0,x\in\mathbb{R}$, where $\xi\in [x_0,x]$.

>What is that $\xi$ (= "xi")? Under certain assumptions elementary functions have their series expansions. If the series is truncated, we have the Taylor polynomial. However, the residual has an explicit expression but due to application of an intermediate value theorem, the exact location of the point $\xi$ is not known, i.e. "generic".

By setting $x:=z+h$, $x_0:=z$, we obtain the useful equivalent formula
$$f(z+h)=f(z)+f'(z)h+\frac{1}{2} f''(\xi)h^2,$$
for every $z,h\in\mathbb{R}$, $\xi\in [z,z+h]$.



## Rate of convergence ($Q$-convergence)

Let $(x_k)$ be an infinite sequence of real numbers.  Let $s_k:=\sup_{l\ge k}x_l$, $k\in\mathbb{N}$, be the *supremum* (i.e., the lowest upper bound) of the *tail* (that is, large indicies $l\ge k$) of $(x_k)$. Define the $\limsup$ (*limes superior*) as
$$\limsup_{k\to\infty}x_k:=\lim_{k\to\infty}s_k\in[-\infty,+\infty].$$
Other than a limit, it always exists, but can be $\pm\infty$. If $(x_k)$ is bounded, the $\limsup$ is the largest limit of a converging subsequence. If $\lim_{k\to\infty} x_k\in (-\infty,\infty)$ exists, then $\lim_{k\to\infty}x_k=\limsup_{k\to\infty}x_k$. The opposite is not true.

>**Examples**. (1) $x_k:=(-1)^k$, then $s_k=1$, and $\limsup_{k\to\infty}x_k=1$.
> (2) $x_k=\sin(k)$, then $s_k=1$, and $\limsup_{k\to\infty}x_k=1$.

Assume that $\lim_{k\to\infty}x_k=x$ and that there exists some large index $M\in\mathbb{N}$ such that $x_k\not=x$ for all $k\ge M$. Then we define the following quantity for $p\ge 0$
$$C(p):=\limsup_{k\to\infty}\frac{|x_{k+1}-x|}{|x_k-x|^p}.$$
We observe that $C(p^*)<\infty$ for some $p^*> 0$ implies $C(p)=0$ for every $0\le p<p^*$. If $C(p^*)>0$ for some $p^*> 0$ then $C(p)=\infty$ for any $p>p^*$.
>**Proof**. By the submultiplicative property of $\limsup$,
>$$\begin{split}C(p)=&\limsup_{k\to\infty}\frac{|x_{k+1}-x|}{|x_k-x|^p}=\limsup_{k\to\infty}\left[\frac{|x_{k+1}-x|}{|x_k-x|^{p^*}}|x_k-x|^{p^*-p}\right]\\\le&\limsup_{k\to\infty}\frac{|x_{k+1}-x|}{|x_k-x|^{p^*}}\limsup_{k\to\infty}|x_k-x|^{p^*-p}=C(p^*)\cdot \begin{cases}0&\text{if}\;\;p<p^*,\\\infty&\text{if}\;\;p>p^*.\end{cases}\end{split}$$

Thus, there exists a (possibly infinite) $p^*$ such that
$$C(p)=\begin{cases}0&\text{if}\;\;0\le p<p^*,\\C(p^*)&\text{if}\;\;p=p^*,\\\infty &\text{if}\;\;p>p^*.\end{cases}$$
The number $p^*$ is called *order of convergence* for the sequence $(x_k)$ and determines the *rate of convergence* as follows:

- If $p^*=1$ and $C(1)=1$ then we say the convergence is *sublinear*. 
- If $p^*=1$ and $1>C(1)>1$ then we say the convergence is *linear*.
- If $p^*>1$ or $C(1)=0$ then we say the convergence is *superlinear*.
- If $p^*=2$ then we say the convergence is *quadratic*.
- If $p^*=3$ then we say the convergence is *cubic*, etc.

When working with convergence estimates it is often useful to use the following approximation:
$$|x_{k+1}-x|\approx C|x_k-x|^{p^*}$$
for some constant $C>0$, not necessarily $C(p^*)$.

>The *rate of convergence* can be used interchangeably with the *order of convergence*. However, there is some caution necessary, as different authors use different terminology here. Usually, the order of convergence always refers to the same thing, namely, the $p^*$-exponent in the denominator of the limit defining the order of convergence. Most confusingly, some authors call the order of convergence "rate of convergence", as e.g. [here](https://www.math-cs.gordon.edu/courses/ma342/handouts/rate.pdf). The English [Wikipedia article](https://en.wikipedia.org/wiki/Rate_of_convergence) calls it the order of convergence, whereas here the rate of convergence is the constant in the definition, which also determines the speed of convergence, together with the order of convergence. So, please always check the context, as the use of the terminology should be clear from it. If there is no definition, try to figure out what is meant in each text. As a rule of thumb: The "order of convergence" is a unique terminology in numerical analysis. The "rate of convergence" can mean at least two different things. I will use both words for the same thing, but will try to make clear what I mean from case to case. In any case, to be sure, use "order of convergence". My PhD advisor usually said that in mathematics "it's all hollow words" (meaning that one should check the definition).

## Little $o$-notation



## Finding roots of functions

Let $f:\mathbb{R}\to\mathbb{R}$ be continuous. We are interested in methods for finding *zeros*, that is, *roots* of $f$, in other words, $x\in\mathbb{R}$, such that $f(x)=0$.

### Bisection method

### Secant method

### Newton's method

## Literature
1. Anne Greenbaum and Tim P. Chartier.  [Numerical Methods: Design, Analysis, and Computer Implementation of Algorithms](https://press.princeton.edu/books/hardcover/9780691151229/numerical-methods), Princeton University Press, 2012.
2. L. Ridgway Scott.  [Numerical Analysis](https://people.cs.uchicago.edu/~ridg/newna/natwo.pdf), Princeton University Press, 2011.
3.  Qingkai Kong, Timmy Siauw, and Alexandre Bayen. [Python Programming and Numerical Methods. A Guide for Engineers and Scientists](https://pythonnumericalmethods.studentorg.berkeley.edu/notebooks/Index.html).  Academic Press, 2020.
4. Tobin A. Driscoll and Richard J. Braun, [Fundamentals of Numerical Computation](https://fncbook.github.io/fnc/intro/floating-point.html), SIAM, 2017.
5. Ernst Hairer, Gerhard Wanner, Syvert P. Nørsett.  Solving Ordinary Differential Equations I: Nonstiff Problems. Springer, 2nd ed., 1993.
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTIzMDUwMDY2OCwxMTI5MzY4MzIwLC0xOT
M0MjY5ODI3LC0xNDAwMDAwNzY0LDE0MDA3MjY2NjQsNDIxMjgy
MzY3LDc2NDY0MDE5MiwxNjc2MDIyMzA5LDkxMzg3MTAxNywxNj
UyNzYxMzA3LDE2NzkxMTIyNjAsMTY2MzM2NTcwNywxMjQ3OTcy
NjEyLC0xOTg3ODU3MDI2LC0xMTE4MDAwNjY4LC00MTQ2MTYwMj
AsMTA4MTU2MDI2NSwtMTU1OTIzODM4MiwtMTAzMDI5OTM1OSwt
MTg3MDE5NTYxM119
-->