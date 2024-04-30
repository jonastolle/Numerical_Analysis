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
>Hence $(1.00)_2=1$, $(1.01)_2=\frac{5}{4}$, $(1.10)_2=\frac{3}{2}$, and >$(1.11)_2=\frac{7}{4}$. By multiplying with the exponents $2^{-1}=\frac{1}{2}$, $2^0=1$, $2^1=2$, we get the whole set:
> $e$ |    |   |    
>---|--|--|--
>$1$ | $\frac{5}{4}$ | $\frac{3}{2}$  | $\frac{7}{4}$
>$2$ | $\frac{5}{2}$ | $3$ | $\frac{7}{2}$
>$\frac{1}{2}$ | $\frac{5}{8}$ | $\frac{3}{4}$ | $\frac{7}{8}$
> Important quantity: $(1.01)_2-1=\frac{1}{4}$, the so-called *machine epsilon*.

Define the machine epsilon by $\varepsilon:=2^{-p}=(1.00\ldots 01)_2-1$.

### Rounding

For the rounding function $\text{round}:\mathbb{R}\to\mathbb{F}$, we have 5 alternative definitions:
- Rounding to nearest (*default*) 
- Rounding to $+\infty$
- Rounding to $-\infty$
- Rounding to $0$
- Rounding away from $0$

It holds that $\text{round}(x)=x(1+\delta)$, where $|\delta|<\frac{\varepsilon}{2}$, where $\varepsilon$ denotes the machine epsilon. Note that usually $\delta$ depends on $x$.
There is a way to define the standard arithmetic operations on $\mathbb{F}$ such that
$$a\oplus b=\text{round}(a+b)=(a+b)(1+\delta_1),$$
$$a\ominus b=\text{round}(a-b)=(a-b)(1+\delta_2),$$
$$a\otimes b=\text{round}(ab)=ab(1+\delta_3),$$
$$a\oslash b=\text{round}\left(\frac{a}{b}\right)=\frac{a}{b}(1+\delta_4),\quad b\not= 0.$$
Here, generally $\delta_i\not=\delta_j$, $i\not=j$.

### IEEE 754 "Double precision"

$k=2$, 64 bits, where:
- The sign: 1 bit;
- The exponent field $0\le\text{ex}\le 2047$: 11 bits, where $e=\text{ex}-1023$, and $-1022\le e\le  1023$, where $\text{ex}=0$ and $\text{ex}=2047$ are special cases.
- The mantissa 52 bits, precision $p=52$.

Exponent field | Number | Type of number
--|---|------
$00\ldots 00=0$ | $\pm (0.b_1 b_2\ldots b_{52})_2 \cdot 2^{-1022}$ | $0$ or *subnormal*
$00\ldots 01=1$ | $\pm (1.b_1 b_2\ldots b_{52})_2 \cdot 2^{-1022}$ | 
$00\ldots 10=2$ | $\pm (1.b_1 b_2\ldots b_{52})_2 \cdot 2^{-1021}$ | 
$\ldots$ | $\ldots$ | 
$01\ldots 11=1023$ | $\pm (1.b_1 b_2\ldots b_{52})_2 \cdot 2^{0}$ |
$\ldots$ | $\ldots$ | 
$11\ldots 10=2046$ | $\pm (1.b_1 b_2\ldots b_{52})_2 \cdot 2^{1023}$ | 
$11\ldots 11=2047$ | $\pm \infty$ if $b_1=b_2=\ldots=b_{52}=0$, otherwise *NaN* |*exception*

Thus, there are two zeros, two infinities and NaN which denotes "*not a number*". The smallest positive normalized number is:
$$(1.0)_2\cdot 2^{-1022}\approx 2.2\cdot 10^{-308}.$$
The largest positive number is:
$$(1.1\ldots 1)_2\cdot 2^{1023}\approx 1.8\cdot 10^{308}$$
The machine epsilon is:
$$2^{-52}\approx 2.22\cdot 10^{-16}.$$

Here's an easy-to-follow [video](https://youtu.be/p8u_k2LIZyo?si=Gfi9TId6x6BAcmpo) explaining floating point numbers (and a specific version of Newton's algorithm).

## Condition number and stability

### Conditioning of problems

Assume that $f:\mathbb{R}\to\mathbb{R}$ "solution map" of the problem, input numbers $x$, $\hat{x}$, close in value, e.g. $\hat{x}=\text{round}(x)$. Set $y:=f(x)$, $\hat{y}:=f(\hat{x})$.

**Definition.** The *absolute condition number* $C(x)$ is defined by the relation
$$|y-\hat{y}|\approx C(x)|x-\hat{x}|.$$
The *relative condition number* $K(x)$ is defined by the relation
$$\left|\frac{y-\hat{y}}{y}\right|\approx K(x)\left|\frac{x-\hat{x}}{x}\right|$$

By the normalization, we guarantee that
$$\text{(relative error in the output)}\approx K(x)\times \text{(relative error in the input)}.$$
Now,
$$y-\hat{y}=f(x)-f(\hat{x})=\underbrace{\frac{f(x)-f(\hat{x})}{x-\hat{x}}}_{\approx f'(x)\;\text{as}\;\hat{x}\to x}(x-\hat{x})$$
Thus, $C(x)=|f'(x)|$.
Furthermore,
$$\frac{y-\hat{y}}{y}=\frac{f(x)-f(\hat{x})}{f(x)}=\underbrace{\frac{f(x)-f(\hat{x})}{x-\hat{x}}}_{\approx f'(x)\;\text{as}\;\hat{x}\to x}\frac{x-\hat{x}}{x}\frac{x}{f(x)}$$
Thus, $K(x)=\left|\frac{x f'(x)}{f(x)}\right|$.

> **Example.** $f(x)=2x$, $f'(x)=2$. Thus, $C(x)=2$, $K(x)=\left|\frac{2x}{2x}\right|=1$. This is a well-conditioned problem.

> **Example.** $g(x)=\sqrt{x}$, $g'(x)=\frac{1}{2\sqrt{x}}$. Thus, $C(x)$ is becomes unbounded for $x>0$ close to zero, e.g. $x\approx 10^{-8}$ yields $C(x)\approx 10^4$. On the other hand, $K(x)=\left|\frac{x}{2\sqrt{x}\sqrt{x}}\right|=\frac{1}{2}$.

### Stability of algorithms

**Definition.** An algorithm or numerical process is called *stable* if small changes in the input produce small changes in the output. It is called *unstable* if large changes in the output are produced by small changes in the input.

An algorithm is stable, if every step is well-conditioned (i.e. has a uniformly bounded condition number). It is unstable if any step is ill-conditioned (i.e. the condition number may become arbitrarily large).

*Forward error analysis (FEA)* is asking:
"How far are we from the true solution?" 

*Backward error analysis (BEA)* is asking:
"Given the answer, what was the problem?"
> **Example.**
> Set:
> $$\text{fl}(x+y):=\text{round}(x)\oplus\text{round}(y)=((x(1+\delta_1)+y(1+\delta_2))(1+\delta_3),$$
> where $|\delta_i|<\frac{\varepsilon}{2}$, $i=1,2,3$.
> FEA:
> $$\text{fl}(x+y)=x+y+x(\delta_1+\delta_3+\delta_1\delta_3)+y(\delta_2+\delta_3+\delta_2\delta_3).$$
> The absolute error is
> $$|\text{fl}(x+y)-(x+y)|\le(|x|+|y|)\left(\varepsilon+\frac{\varepsilon^2}{4}\right).$$
> The relative error is:
>  $$\left|\frac{\text{fl}(x+y)-(x+y)}{x+y}\right|\le\frac{(|x|+|y|)\left(\varepsilon+\frac{\varepsilon^2}{4}\right)}{|x+y|}.$$
> BEA:
> $$\text{fl}(x+y)=x(1+\delta_1)(1+\delta_3)+y(1+\delta_2)(1+\delta_3).$$
> Thus the relative error for each term is less or equal to $\varepsilon+\frac{\varepsilon^2}{4}$.
> Hence the sum of two floating point numbers is backwards stable.

Well-conditioned problems may have unstable algorithms. For stability, each step has to be well conditioned. Some ill-conditioned steps produce an unstable algorithm. Ill-conditioned problems cannot be reliably solved with a stable algorithm.

>**Example.** Consider evaluating $f(x)=\sqrt{1+x}-1$ for $x$ close to zero. The relative condition number is:
>$$K(x)=\left|\frac{xf'(x)}{f(x)}\right|=\frac{x}{2\sqrt{1+x}(\sqrt{1+x}-1)}$$
>$$=\frac{x(\sqrt{1+x}+1)}{2\sqrt{1+x}(\sqrt{1+x}-1)(\sqrt{1+x}+1)}=\frac{\sqrt{1+x}+1}{2\sqrt{1+x}},$$
>and $K(0)=1$.
>Consider the following 3 steps;
>1. $t_1:=1+x$, well-conditioned, $x$ close to $0$.
>2. $t_2:=\sqrt{t_1}$, relatively well-conditioned, also absolutely well conditioned, because $t_1$ is close to $1$.
>3. $t_3:=t_2-1$, ill-conditioned, relative condition number of this step: $K_3(t_2)=\left|\frac{t_2}{t_2-1}\right|$, which becomes unbounded for $t_2$ close to $1$!
>On the other hand, the problem is well-conditioned. Solve it by writing:
>$$f(x)=\sqrt{1+x}-1=\frac{(\sqrt{1+x}+1)(\sqrt{1+x}-1)}{\sqrt{1+x}+1}=\frac{x}{\sqrt{1+x}+1},$$
>which can be evaluated directly close to zero.

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
We observe that $C(p^{*})<\infty$ for some $p^{*}> 0$ implies $C(p)=0$ for every $0\le p<p^{*}$. If $C(p^{*})>0$ for some $p^{*}> 0$ then $C(p)=\infty$ for any $p>p^{*}$.
>**Proof**. By the submultiplicative property of $\limsup$,
>$$C(p)=\limsup_{k\to\infty}\frac{|x_{k+1}-x|}{|x_k-x|^p}=\limsup_{k\to\infty}\left[\frac{|x_{k+1}-x|}{|x_k-x|^{p^{*}}}|x_k-x|^{p^{*}-p}\right]$$
>$$\le\limsup_{k\to\infty}\frac{|x_{k+1}-x|}{|x_k-x|^{p^{*}}}\limsup_{k\to\infty}|x_k-x|^{p^{*}-p}=C(p^{*})\cdot \begin{cases}0&\text{if}\;\;p<p^{*},\\\infty&\text{if}\;\;p>p^{*}.\end{cases}$$

Thus, there exists a (possibly infinite) $p^{*}$ such that
$$C(p)=\begin{cases}0&\text{if}\;\;0\le p<p^{*},\\C(p^{*})&\text{if}\;\;p=p^{*},\\\infty &\text{if}\;\;p>p^{*}.\end{cases}$$
The number $p^{*}$ is called *order of convergence* for the sequence $(x_k)$ and determines the *rate of convergence* as follows:

- If $p^{*}=1$ and $C(1)=1$ then we say the convergence is *sublinear*. 
- If $p^{*}=1$ and $1>C(1)>1$ then we say the convergence is *linear*.
- If $p^{*}>1$ or $C(1)=0$ then we say the convergence is *superlinear*.
- If $p^{*}=2$ then we say the convergence is *quadratic*.
- If $p^{*}=3$ then we say the convergence is *cubic*, etc.

When working with convergence estimates it is often useful to use the following approximation:
$$|x_{k+1}-x|\approx C|x_k-x|^{p^{*}}$$
for some constant $C>0$, not necessarily $C(p^{*})$.
Here, it is useful to look at the logarithmic behavior:
$$\log(|x_{k+1}-x|)\approx\log\left(C|x_k-x|^{p^{*}}\right)=\log(C)+\log\left(|x_k-x|^{p^{*}}\right)=\log(C)+p^{*}\log(|x_k-x|).$$

>The *rate of convergence* can be used interchangeably with the *order of convergence*. However, there is some caution necessary, as different authors use different terminology here. Usually, the order of convergence always refers to the same thing, namely, the $p^{*}$-exponent in the denominator of the limit defining the order of convergence. Most confusingly, some authors call the order of convergence "rate of convergence", as e.g. [here](https://www.math-cs.gordon.edu/courses/ma342/handouts/rate.pdf). The English [Wikipedia article](https://en.wikipedia.org/wiki/Rate_of_convergence) calls it the order of convergence, whereas here the rate of convergence is the constant in the definition, which also determines the speed of convergence, together with the order of convergence. So, please always check the context, as the use of the terminology should be clear from it. If there is no definition, try to figure out what is meant in each text. As a rule of thumb: The "order of convergence" is a unique terminology in numerical analysis. The "rate of convergence" can mean at least two different things. I will use both words for the same thing, but will try to make clear what I mean from case to case. In any case, to be sure, use "order of convergence". My PhD advisor usually said that in mathematics "it's all hollow words" (meaning that one should check the definition).

## Landau's little $o$- and big $O$-notation

Copied from [Wikibooks](https://en.wikibooks.org/wiki/Real_Analysis/Landau_notation) under a Creative Commons BY-SA 4.0 license.

The Landau notation is an amazing tool applicable in all of real analysis. The reason it is so convenient and widely used is because it underlines a key principle of real analysis, namely ''estimation''. Loosely speaking, the Landau notation introduces two operators which can be called the "order of magnitude" operators, which essentially compare the magnitude of two given functions.

### The "little-$o$"

The "little-$o$" provides a function that is of lower order of magnitude than a given function, that is the function $o(g(x))$ is of a lower order than the function $g(x)$. Formally,

**Definition.**
Let $A\subseteq\mathbb{R}$ and let $c\in\mathbb{R}$.
Let $f,g:A\to\mathbb{R}$.

If $\lim_{x\to c}\frac{f(x)}{g(x)}=0$ then we say that

"As $x\to c$, $f(x)=o(g(x))$"

**Examples.**
- As $x\to\infty$, (and $m<n$), $x^m=o(x^n)$;
- As $x\to\infty$, (and $n\in\mathbb{N}$), $\log x=o(x^n)$;
- As $x\to 0$, $\sin x=o(1)$.

### The "Big-$O$"
The "Big-$O$" provides a function that is at most the same order as that of a given function, that is the function $O(g(x))$ is at most the same order as the function $g(x)$. Formally,

**Definition.**
Let $A\subseteq\mathbb{R}$ and let $c\in\mathbb{R}$

Let $f,g:A\to\mathbb{R}$

If there exists $M>0$ such that $\lim_{x\to c}\left| \frac{f(x)}{g(x)}\right| <M$ then we say that

"As $x\to c$, $f(x)=O(g(x))$"

**Examples ===
* As <math>x\to 0</math>, <math>\sin x=O(x)</math>
* As <math>x\to \tfrac{\pi}{2}</math>, <math>\sin x=O(1)</math>

**Applications.**
We will now consider few examples which demonstrate the power of this notation.
Differentiability ===

Let <math>f: \mathcal{U} \subseteq \mathbb{R} \to\mathbb{R}</math> and <math> x_0 \in \mathcal{U}</math>.

Then <math>f</math> is differentiable at <math>x_0</math> if and only if

There exists a <math>\lambda \in\mathbb{R}</math> such that as <math>x\to x_0</math>, <math>f(x) = f(x_0) + \lambda(x-x_0)+o\left( |x-x_0|\right)</math>.

=== Mean Value Theorem ===
Let <math>f:[a,x]\to\mathbb{R}</math> be differentiable on <math>[a,b]</math>. Then,

As <math>x\to a</math>, <math>f(x)=f(a)+O(x-a)</math>

=== Taylor's Theorem ===
Let <math>f:[a,x]\to\mathbb{R}</math> be ''n''-times differentiable on <math>[a,b]</math>. Then,

As <math>x\to a</math>, <math>f(x)=f(a)+\tfrac{(x-a)f'(a)}{1!}+\tfrac{(x-a)^2f''(a)}{2!}+\ldots+\tfrac{(x-a)^{n-1}f^{(n-1)}(a)}{(n-1)!}+O\left( (x-a)^n\right)</math>


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
eyJoaXN0b3J5IjpbLTUwNTk0MjMyLC0xMjE5NDQ5NzQsLTEyMT
k0NDk3NCwtMjk5ODgxODA2LC0xMjIxNzE2NzY5LDYxMDY2MzM2
NiwtMTIyMTcxNjc2OSw4ODI5NzA5NzksMTYzOTUxMTAwOSwtMj
EyNzI1MDk2LDExMjkzNjgzMjAsLTE5MzQyNjk4MjcsLTE0MDAw
MDA3NjQsMTQwMDcyNjY2NCw0MjEyODIzNjcsNzY0NjQwMTkyLD
E2NzYwMjIzMDksOTEzODcxMDE3LDE2NTI3NjEzMDcsMTY3OTEx
MjI2MF19
-->