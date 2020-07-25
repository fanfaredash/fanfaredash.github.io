---
layout: post
title: AGC034F RNG and XOR
subtitle: 
date: 2020-03-14 00:00:00
author: fpdqwq
header-img: /img/home.png
mathjax: true
catalog: true
tags: 算法 题解
public: true
---

## 前言

### 前置知识

- 期望方程的设立
- 线性代数（矩阵乘法、矩阵的逆）
- 快速沃尔什变换（FWT）
- 一些其它反演

以上为阅读本文章的前置知识，而并非解决本题必备知识。

### 相关阅读

- [yyb前辈的FWT学习笔记](https://www.cnblogs.com/cjyyb/p/9065615.html)
- VFK的《炫酷反演魔术》

本文取材自上述两篇博客/讲义/文章。

## 题解

### 记号与概念

- 反演：已知 $f(n)=\sum_{k}{a_{n,k}g(k)}$, 由 $f$ 逆向推演出 $g$ 的过程。
- 记列向量 $\alpha=(f(0),f(1),…f(n-1))$, $\beta=(g(0),g(1),…g(n-1))$, 矩阵 $A=(a_{i,j})_{n,n}$
- 显然 $\alpha = A\beta$。我们知道 $A, \alpha$ 求 $\beta$ 可以使用高斯消元求逆矩阵，但是复杂度较高，难以承受。
- 一般我们考察的反演都极具规律性，也就是 $A, A^{-1}$ 不仅已知，而且易于计算（与向量的）乘法。例如二项式反演、莫比乌斯反演、快速傅里叶变换等。
- 记号：在本文中，$ \alpha \times \beta $ 不再表示叉乘，而表示向量每个维度分别相乘得到的向量。
- 记号：在本文中，记 $A$ 为 xor_FWT 的反演矩阵，即用 $A\alpha$ 来代指 $FWT(\alpha)$。

### 题解

- 列期望方程如下：
$$
\begin{cases}
E_0 = 0 \\
E_i = 1 + \sum_j{E_j \ p_{i\ xor\ j}}
\end{cases}
$$
其中 $p_i = A_i/S$
- 像我们预期的那样，高斯消元效率过低。
- 希望能够利用 FWT 对运算进行加速，将方程配出 FWT_xor 的形式：
$$
\begin{cases}
E_0 + \sum_j{E_j \ p_j} &= \sum_j{E_j \ p_j} \\
E_i - 1 &= \sum_j{E_j \ p_{i\ xor\ j}}
\end{cases}
$$
- 得到如下方程：
$$
\epsilon - \beta = A^{-1}((A\epsilon)\times(A\alpha))
$$
其中 $\epsilon=(E_0, E_1, \dots E_{2^N-1})$, $\alpha = (p_0, p_1, \dots p_{2^N-1})$,$\beta=(-\sum_j{E_j \ p_j}, 1,\dots1)$
- 移项，相乘可得：
$$
A\beta = (A\epsilon)\times(A\alpha - (1, 1,\dots 1))
$$
- 注意到非常关键的一点（作者在这里 WA 了很久）：$A\alpha-1$ 的第一个维度为 0。这由反演矩阵 $A$ 本身和 $\sum p_i=1$ 决定。也就是说，$A\beta$ 的第一个维度亦为 0。这意味着 $A\beta$ 变为已知，求解方法如下——
- 由于反演是线性变换，设 $\beta$ 的第一个维度为 $x$, 那么经过 FWT 变换之后，向量的每个维度应该是关于 $x$ 的线性函数；我们只要把取这个线性函数在某点，使得 $(A\beta)_0=0$ 即可。
- 到这里我们几乎要完成这道题了，因为 $A\beta, A\alpha$ 除了第一个维度以外都已知，由此我们可以计算出 $A\epsilon$ 除了第一个维度以外的所有值。由 $E_0=0$，用和求解 $A\beta$ 一样的方法可以求出 $\epsilon$。

有趣的是，$0x=0$ 在其他题目中是烦人的 edge case，在本题中反而成为了突破口。

### 代码实现

#### 代码

> [AGC034F: RNG and XOR](https://atcoder.jp/contests/agc034/tasks/agc034_f)

```cpp
#include <bits/stdc++.h>

using namespace std;

typedef pair <int, int> pi;
typedef long long LL;   

const int P = 998244353;
const int INF = 0x3f3f3f3f;
const LL INFLL = 0x3f3f3f3f3f3f3f3fll;
const int N = 1e6 + 1e5;

int n, m, inv_n;
LL a[N], b[N], c[N], S, inv_S;

inline LL ksm(LL base, LL x) {
    LL res = 1;
    while(x) {
        if(x & 1) res *= base, res %= P;
        x >>= 1, base *= base, base %= P;
    }
    return res;
}

inline void fwt_xor(LL a[], int type) {
    for(int i = 1; i < n; i <<= 1)
        for(int j = 0; j < n; j += i << 1)
            for(int k = j; k < j + i; k++)
                a[k] = (a[k] + a[k + i]) % P,
                a[k + i] = (a[k] - a[k + i] * 2 + P * 2) % P;
    if(type == -1)
        for(int i = 0; i < n; i++) a[i] = (a[i] * inv_n) % P;
}

void init() {
    int len;
    cin >> len;
    n = 1 << len;
    inv_n = ksm(n, P - 2);
    for(int i = 0; i < n; i++) scanf("%lld", &a[i]), S = (S + a[i]) % P;
    inv_S = ksm(S, P - 2);
    for(int i = 0; i < n; i++) a[i] = (a[i] * inv_S) % P;
    fwt_xor(a, 1);
    for(int i = 0; i < n; i++) b[i] = c[i] = 1;
    b[0] = 0;
    fwt_xor(b, 1), fwt_xor(c, 1);

    for(int i = 0; i < n; i++) c[i] = (c[i] - b[i] + P) % P;
    LL x = (P - b[0]) * ksm(c[0], P - 2) % P;
    for(int i = 0; i < n; i++) c[i] = (b[i] + c[i] * x) % P;
    for(int i = 0; i < n; i++) {
        LL tmp = ksm((1 - a[i] + P) % P, P - 2);
        c[i] = (c[i] * tmp) % P;
    }
    for(int i = 0; i < n; i++) b[i] = c[i];
    c[0] = 1;
    fwt_xor(b, -1), fwt_xor(c, -1);

    
    for(int i = 0; i < n; i++) c[i] = (c[i] - b[i] + P) % P;
    x = (P - b[0]) * ksm(c[0], P - 2) % P;
    for(int i = 0; i < n; i++) c[i] = (b[i] + c[i] * x) % P;
    for(int i = 0; i < n; i++) printf("%lld\n", c[i]);
}

int main() {
    init();
    return 0;
}

```

Author : Fanfaredash

![71257925_p0.png](https://i.loli.net/2020/03/15/IHgm7JyA5xCjpsO.png)
