# Numerical Analysis

**Jonas Tölle**[^1]

[^1]: [jonas.tolle@aalto.fi](mailto:jonas.tolle@aalto.fi)

*Lecture notes for MS-C1650 Numerical Analysis, Aalto University*

*Last updated: 27.5.2024*

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
$$y-\hat{y}=f(x)-f(\hat{x})=\underbrace{\frac{f(x)-f(\hat{x})}{x-\hat{x}}}_{\to f'(x)\;\text{as}\;\hat{x}\to x}(x-\hat{x})$$
Thus, $C(x)=|f'(x)|$.
Furthermore,
$$\frac{y-\hat{y}}{y}=\frac{f(x)-f(\hat{x})}{f(x)}=\underbrace{\frac{f(x)-f(\hat{x})}{x-\hat{x}}}_{\to f'(x)\;\text{as}\;\hat{x}\to x}\frac{x-\hat{x}}{x}\frac{x}{f(x)}$$
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
>$$\le\left(\limsup_{k\to\infty}\frac{|x_{k+1}-x|}{|x_k-x|^{p^{*}}}\right)\cdot\left(\limsup_{k\to\infty}|x_k-x|^{p^{*}-p}\right)=C(p^{*})\cdot \begin{cases}0&\text{if}\;\;p<p^{*},\\\infty&\text{if}\;\;p>p^{*}.\end{cases}$$
>$\Box$

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

**Examples.**
- As $x\to 0$, $\sin x=O(x)$;
- As $x\to \tfrac{\pi}{2}$, $\sin x=O(1)$.

### Applications

We will now consider few examples which demonstrate the power of this notation.

#### Differentiability

Let $f: \mathcal{U} \subseteq \mathbb{R} \to\mathbb{R}$ and $x_0 \in \mathcal{U}$.

Then $f$ is differentiable at $x_0$ if and only if

There exists a $\lambda \in\mathbb{R}$ such that as $x\to x_0$, $f(x) = f(x_0) + \lambda(x-x_0)+o\left( |x-x_0|\right)$.

#### Mean Value Theorem
Let $f:[a,x]\to\mathbb{R}$ be differentiable on $[a,b]$. Then,

As $x\to a$, $f(x)=f(a)+O(x-a)$.

#### Taylor's Theorem
Let $f:[a,x]\to\mathbb{R}$ be $n$-times differentiable on $[a,b]$. Then,

As $x\to a$, $f(x)=f(a)+\tfrac{(x-a)f'(a)}{1!}+\tfrac{(x-a)^2f''(a)}{2!}+\ldots+\tfrac{(x-a)^{n-1}f^{(n-1)}(a)}{(n-1)!}+O\left( (x-a)^n\right)$.


## Finding roots of functions and fixed points

Let $f:\mathbb{R}\to\mathbb{R}$ be continuous. We are interested in methods for finding *zeros*, that is, *roots* of $f$, in other words, $x\in\mathbb{R}$, such that $f(x)=0$.

**Definition.** If $f:\mathbb{R}^n\to\mathbb{R}$, $n\in\mathbb{N}$ is $k$-times continuously differentiable, then we write $f\in C^k(\R^n)$ or just $f\in C^k$. For $k=0$, $f\in C^0$ or $f\in C(\R^n)$ or $f\in C([a,b])$ just means that $f$ is assumed to be continuous.

### Bisection method

The intermediate value theorem for continuous functions implies that $x_1<x<x_2$ with $f(x)=0$ exists if $f(x_1) f(x_2)<0$, i.e., there is a sign change. The *bisection method* is based on halving the interval such that the sign condition is preserved. Note that, in principle, we have to look for intervals $[x_1,x_2]$.

Let us analyze the convergence rate. Let $[a,b]$ be an interval. After $k$ steps the interval of analysis has length $\frac{b-a}{2^k}$ which converges to zero for $k\to\infty$. Let us look in a neighborhood of radius $\delta>0$, so that
$$\frac{b-a}{2^k}\le 2\delta\quad\Leftrightarrow\quad 2^{k+1}\ge \frac{b-a}{\delta}\quad\Leftrightarrow\quad k\ge \log_2\left(\frac{b-a}{\delta}\right)-1.$$
Every step reduces the error by factor $\frac{1}{2}$. The convergence rate is thus linear.

### Newton's method

Assume that $f\in C^1$. For an intial value $x_0$, consider the iteration
$$x_{k+1}=x_k-\frac{f(x_k)}{f'(x_k)}\quad k=0,1,\ldots.$$
> **Heuristics.** If $f(x_\ast)=0$ and $f\in C^2$, by the Taylor's expansion,
> $$0=f(x_\ast)=f(x_0)+(x_\ast-x_0)f'(x_0)+\frac{(x_\ast-x_0)^2}{2}f''(\xi)$$
> for $\xi\in [x_0,x_\ast]$, and upon neglecting the 2nd order term,
> $$x_\ast \approx x_0-\frac{f(x_0)}{f'(x_0)}.$$

**Theorem.** If $f\in C^2$ and $x_0$ is sufficiently good (i.e. close to the root $x_\ast$) and if $f'(x_\ast)\not=0$, then Newton's method converges quadratically.

**Proof.** By Taylor's expansion, it follows that
$$x_\ast=x_k-\frac{f(x_k)}{f'(x_k)}-\frac{(x_\ast-x_k)^2}{2}\frac{f''(\xi_k)}{f'(x_k)}$$
for some $\xi_k\in [x_\ast,x_k]$. Take $x_{k+1}$ from the method and substract,
$$x_{k+1}-x_\ast=(x_\ast-x_k)^2\underbrace{\frac{f''(\xi_k)}{2f'(x_k)}}_{\le D}.$$
In other words,
$$|x_{k+1}-x_\ast|\le D |x_k-x^\ast|^2,$$
as $k\to\infty$ and thus $x_k\to x_\ast$.
Hence the method is quadratic. Note that $f'(x_k)$ does not vanish by continuity if $x_k$ is close to $x_\ast$. $\Box$

*What happens if $f'(x_\ast)=0$?*
$$x_{k+1}-x_\ast=(x_\ast-x_k)^2\frac{f''(\xi_k)}{2\underbrace{f'(x_k)}_{\to0}}$$
as $k\to\infty$.
By Taylor's expansion,
$$f'(x_k)=\underbrace{f'(x_\ast)}_{=0}+(x_k-x_\ast)f''(\eta_k)=(x_k-x_\ast)f''(\eta_k)$$
for some $\eta_k\in [x_\ast,x_k]$, and hence
$$x_{k+1}-x_\ast=f''(\eta_k)=(x_k-x_\ast)f''(\eta_k)(x_k-x_\ast).$$
**The method has degenerated to a linear method!**

> **Example.** $f(x)=x^2$, $f'(x)=2x$. Newton:
> $$x_{k+1}=x_k-\frac{x^2_k}{2x_k}=\frac{1}{2}x_k.$$

### Secant method

Sometimes it can be difficult or computationally expensive to compute the derivative $f'(x_k)$. Newton 's method can be adapted by approximating the derviative by the differential quotient. The secant method is the following two-step recursive algorithm.
$$x_{k+1}=x_k-\frac{f(x_k)(x_k-x_{k-1})}{f(x_k)-f(x_{k-1})},\quad k=1,2,\ldots$$
with two distinct starting points $x_0\not=x_1$. The convergence rate is $\frac{1+\sqrt{5}}{2}\approx 1.62$, the *golden ratio*. 

### Fixed point iteration

**Definition.** A point $x\in\mathbb{R}$ is called a fixed point of $\varphi:\mathbb{R}\to\mathbb{R}$ if $\varphi(x)=x$.

We could for instance use Newton's method to find fixed points by setting $f(x):=\varphi(x)-x$.

**Banach's Fixed Point Theorem.** Suppose that $\varphi$ is a *contraction*, that is, there exists a constant $L<1$ such that
$$|\varphi(x)-\varphi(y)|\le L|x-y|$$
for all $x,y\in\mathbb{R}$. Then there exists a unique fixed point $x_\ast\in\mathbb{R}$ of $\varphi$, i.e., $\varphi(x_\ast)=x_\ast$, and the fixed point iteration $\varphi_n:=\underbrace{\varphi\circ\ldots\circ\varphi}_{n\text{-times}}$ satisfies $\lim_{n\to\infty}\varphi_n(x_0)=x_{\ast}$ for any starting point $x_0\in\mathbb{R}$. The convergence rate is at least linear.

**Proof.**  We prove that the  sequence $\{\varphi_k(x_0)\}_{k=0}^{\infty}=\{x_k\}_{k=0}^{\infty}$ is a Cauchy sequence. Let $k>j$. Then, by the triangle inequality,
$$|x_k-x_j|\le\underbrace{|x_k-x_{k-1}|+|x_{k-1}-x_{k-2}|+\ldots+|x_{j+1}-x_j|}_{(k-j)\text{-summands}}.$$
Furthermore,
$$|x_m-x_{m-1}|=|\varphi(x_{m-1})-\varphi(x_{m-2})|\le L|x_{m-1}-x_{m-2}|\le L^{m-1}|x_1-x_0|.$$
Hence, by the geometric series,
$$|x_k-x_j|\le L^j\frac{1-L^{k-j}}{1-L}|x_1-x_0|.$$
If $k>N$, $j>N$, then
$$|x_k-x_j|\le L^N\frac{1}{1-L}|x_1-x_0|\to 0$$
as $N\to\infty$, which proves that $\{x_k\}$ is a Cauchy sequence. The linear convergence rate follows also from this estimate.
The existence of a fixed point follows from the continuity of $\varphi$ (as contractions are uniformly continuous, in fact, even Lipschitz continuous) as follows.
$$x_\ast=\lim_{k\to\infty}x_k=\lim_{k\to\infty}x_{k+1}=\lim_{k\to\infty}\varphi(x_k)=\varphi(\lim_{k\to\infty}x_k)=\varphi(x_\ast).$$
$\Box$

**Theorem**. Assume that $\varphi\in C^p$ for $p\ge 1$. Furthermore, assume that has a fixed point $x_\ast$ and assume that
$$\varphi'(x_\ast)=\varphi''(x_\ast)=\ldots=\varphi^{(p-1)}(x_\ast)=0$$ for $p\ge 2$ and
$$G'(x_\ast)<1$$ if $p=1$. Then  the fixed point sequence $\{\varphi_k(x_0)\}$ converges to $x_\ast$ at least with rate $p$, provided that the starting point $x_0$ is sufficiently close to $x_\ast$. If, in addition, $\varphi^{(p)}(x_\ast)\not=0$, then the rate of convergence is precisely $p$.

**Proof.** First note that by Banach's fixed point theorem the limit indeed converges to $x_\ast$ for suitable starting points $x_0$. By the Taylor expansion,
$$x_{k+1}-x_\ast=\varphi(x_k)-\varphi(x_\ast)=\sum_{l=1}^{p-1}\frac{\varphi^{(l)}(x_\ast)}{l!}(x_k-x_\ast)^l+\frac{G^{(p)}(\xi_k)}{p!}(x_k-x_\ast)^p$$
for some $\xi_k$ between $x_\ast$ and $x_k$. The sum will be left empty for the case $p=1$. Since $\varphi^{(l)}(x_\ast)=0$ for $1\le l\le p-1$, we get that
$$|x_{k+1}-x_\ast|=\frac{|\varphi^{(p)}(\xi_k)|}{p!} |x_k-x_\ast|^p$$.
By continuity, there exists $C>0$ (with $C<1$ for $p=1$) such that
$$\frac{|\varphi^{(p)}(\xi_k)|}{p!}\le C$$ for $\xi_k$ sufficiently close to $x_\ast$ (that is, for sufficiently large $k$). Thus,
$$|x_{k+1}-x_\ast|\le C|x_k-x_\ast|^p$$
for large $k$, and thus the rate of convergence is at least $p$. Note that  for $p=1$, this also proves convergence by
$$|x_{k+1}-x_\ast|<\underbrace{|\varphi(\xi_k)|}_{<1}|x_k-x_\ast|.$$
If $\varphi^{(p)}\not=0$, then by continuity, there exists $K>0$ such that
$$\frac{|\varphi^{(p)}(\xi_k)|}{p!}\ge K$$
for $\xi_K$ sufficiently close to $x_\ast$. Thus
$$|x_{k+1}-x_\ast|\ge K|x_k-x_\ast|^p$$
which implies that the rate of convergence cannot be higher than $p$. Thus the rate of convergence is percisely $p$. $\Box$

> **Note.** From the above proof, we expect that close to the fixed point $x_\ast$
> $$|x_{k+1}-x_\ast|\approx \frac{|\varphi^{(p)}(x_\ast)|}{p!} |x_k-x_\ast|^p,$$
> when
> $$\varphi'(x_\ast)=\varphi''(x_\ast)=\ldots=\varphi^{(p-1)}(x_\ast)=0,$$
> but $\varphi^{(p)}(x_\ast)\not=0$.

## Polynomial interpolation

**Idea.** Approximate a function $f:\mathbb{R}\to\mathbb{R}$ over $[a,b]$ by a polynomial $p$ such that in **distinct** data points $(x_i,y_i)$, $i=0,1,\ldots, n$, the approximation is *exact*, that is,
$$f(x_i)=y_i=p(x_i),\quad\text{for all}\quad  i=0,1,\ldots, n.$$
We may call $x_i$ *node* and $y_i$ *value*.
We need at least $2$ data points. We usually just assume that $x_i\not=x_j$ for $i\not=j$.

> **Note.** Interpolation polynomials are not per se unique, for instance the data $\{(-1,1),(1,1)\}$ can be interpolated by
> $p(x)=1$, $q(x)=x^2$, or $r(x)=x^4-x^2+1$. However, we will see later that $p$ is the unique interpolation polynomial with $\deg p\le 1= n$.

**Example.** $(1,2)$, $(2,3)$, $(3,6)$, as data set $\{(x_i,y_i)\;\colon\;i=0,1,2\}$ on the interval $[1,3]$.
We are looking for a polynomial $p_2(x)=\sum_{j=0}^2 c_j x^j$, which is chosen to be 2nd order, because we have $3$ data points and $3$ unknown coeffiecients.
We can formulate the problem in matrix form:
$$\begin{pmatrix}1 & x_0 & x_0^2\\ 1 & x_1 & x_1^2\\ 1 &x_2 & x_2^2\end{pmatrix}\cdot \begin{pmatrix}c_0\\c_1\\c_2\end{pmatrix}=\begin{pmatrix}y_0\\y_1\\y_2\end{pmatrix},$$
which is a Vandermonde matrix which has determinant
$$\det\begin{pmatrix}1 & x_0 & x_0^2\\ 1 & x_1 & x_1^2\\ 1 &x_2 & x_2^2\end{pmatrix}=\prod_{i<j}(x_j-x_i)\not=0,$$ and is thus invertible. Here,
$$\begin{pmatrix}1 & 1 & 1\\ 1 & 2 & 4\\ 1 &3 & 9\end{pmatrix}\cdot \begin{pmatrix}c_0\\c_1\\c_2\end{pmatrix}=\begin{pmatrix}2\\3\\6\end{pmatrix}.$$
As a result, $c_0=3$, $c_1=-2$, and $c_2=1$, and thus,
$$p_2(x)=x^2-2x+3.$$

The computational complexity of solving the linear system is $O(n^3)$. We used the natural basis for the polynomials.
*What would be the ideal basis?*
**Definition.** (Lagrange basis polynomials) Suppose that $x_i\not=x_j$ if $i\not=j$. We call
$$\phi_i(x):=\prod_{\substack{j=0\\ i\not=j}}^n\frac{x-x_j}{x_i-x_j}$$
the $i$th *Lagrange basis polynomial*.
The *Lagrange interpolation polynomial* is given by
$$p(x):=\sum_{i=0}^n y_i \varphi_i(x).$$

Clearly,
$$\varphi_i(x_j)=\delta_{i,j}:=\begin{cases}1\text{\;\;if\;\;}i=j,\\0\text{\;\;if\;\;}i\not=j.\end{cases}$$

> **Example.** $(1,2)$, $(2,3)$, $(3,6)$:
> $$\varphi_0(x)=\frac{(x-2)(x-3)}{(1-2)(1-3)}$$
> $$\varphi_1(x)=\frac{(x-1)(x-3)}{(2-1)(2-3)}$$
> $$\varphi_2(x)=\frac{(x-1)(x-2)}{(3-1)(3-2)}$$
> $$p_2(x)=2\varphi_0(x)+3\varphi_1(x)+6\varphi_2(x)=x^2-2x+3.$$

Evaluating the Lagrange polynomials has the computational complextity $O(n^2)$.

### Newton's interpolation

**Idea.** Extend the natural basis:
$$1,\quad x-x_0,\quad (x-x_0)(x-x_1),\quad\ldots,\quad \prod_{j=0}^{n-1}(x-x_j).$$

**Definition.** Define Newton's interpolation polynomials by
$$p_n(x)=a_0+a_1(x-x_0)+\ldots+a_n\prod_{j=0}^{n-1} (x-x_j)$$
in such a way that $p_n(x_i)=y_i$.
Clearly,
$$p(x_0)=y_0\quad\Rightarrow\quad a_0=y_0,$$
and
$$p(x_1)=a_0+a_1(x_1-x_0)=y_1\quad \Rightarrow \quad a_1=\frac{y_1-a_0}{x_1-x_0}.$$
More generally, we have the lower triangular linear system
$$\begin{pmatrix}
1& 0 &\cdots&&&&\\
1& x_1-x_0 & 0&\cdots&&&\\
1& x_1-x_0 & (x_2-x_0)(x_2-x_1)&0&\cdots&&\\
\vdots &\vdots & \vdots&\vdots&&&\\
1& x_1-x_0 & (x_2-x_0)(x_2-x_1)&&\cdots&&\prod_{j=0}^{n-1}(x_n-x_j)
\end{pmatrix}\begin{pmatrix}a_0\\a_1\\a_2\\\vdots\\a_n\end{pmatrix}=\begin{pmatrix}y_0\\y_1\\y_2\\\vdots\\y_n\end{pmatrix}.$$

> **Example.** $$p_2(x)=a_0+a_1(x-1)+a_2(x-1)(x-2)$$ with the system
> $$\begin{pmatrix}1&0&0\\ 1& 1&0\\ 1&2&2\end{pmatrix}\cdot\begin{pmatrix}a_0\\a_1\\a_2\end{pmatrix}=\begin{pmatrix}2\\3\\6\end{pmatrix}$$
> and hence $a_0=2$, $a_1=1$, and $a_2=1$ which yields $p_2(x)=x^2-2x+3$.

### Uniqueness

**Theorem.** Interpolation polynomials with $n+1$ nodes $x_i$, $i=0,1,\ldots, n$ are unique in the class of polynomials $q$ with $\deg q\le n$.

**Proof.** (Idea). $p_n$ has at most $n$ roots. Let $p_n$ and $q_n$ be two in

## Literature
1. Anne Greenbaum and Tim P. Chartier.  [Numerical Methods: Design, Analysis, and Computer Implementation of Algorithms](https://press.princeton.edu/books/hardcover/9780691151229/numerical-methods), Princeton University Press, 2012.
2. L. Ridgway Scott.  [Numerical Analysis](https://people.cs.uchicago.edu/~ridg/newna/natwo.pdf), Princeton University Press, 2011.
3.  Qingkai Kong, Timmy Siauw, and Alexandre Bayen. [Python Programming and Numerical Methods. A Guide for Engineers and Scientists](https://pythonnumericalmethods.studentorg.berkeley.edu/notebooks/Index.html).  Academic Press, 2020.
4. Tobin A. Driscoll and Richard J. Braun, [Fundamentals of Numerical Computation](https://fncbook.github.io/fnc/frontmatter.html), SIAM, 2017.
5. Ernst Hairer, Gerhard Wanner, Syvert P. Nørsett.  Solving Ordinary Differential Equations I: Nonstiff Problems. Springer, 2nd ed., 1993.
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTExODM4MDgzNTUsLTIxMDYxOTk1ODMsMT
Y0NTIwMzc4NSwtMzM2NTc1MjYxLC0xOTkyNDMzNzQ5LDEwMTc4
Mjk2NzQsLTEyODYzNjk3OCwzOTY1NDYwMzAsLTE2NDc3NzM0MT
UsLTE4MTAxMzIxOTMsLTg5MjU4ODg2MSwtMTM5MzAxMDI2OSwt
MTY0MzE5MTc2NywtMTIxOTQ0OTc0LC0xMjE5NDQ5NzQsLTI5OT
g4MTgwNiwtMTIyMTcxNjc2OSw2MTA2NjMzNjYsLTEyMjE3MTY3
NjksODgyOTcwOTc5XX0=
-->