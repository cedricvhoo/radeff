# Theory

## Radiation efficiency for rectangular geometries

The geometrical radiation efficiency can be calculated from:
```{math}
\sigma(k_\mathrm{f},\phi) = \mathrm{Re}\left\{\dfrac{\mathrm{i}k_0}{S}\int_S\int_SG(\mathbf{x},\mathbf{x}')\mathrm{e}^{\mathrm{i}k_\mathrm{f}\left[(x-x')\cos\phi+(y-y')\sin\phi\right]}\mathrm{d}\mathbf{x}\mathrm{d}\mathbf{x}'\right\},
```
with {math}`\mathbf{x} = [x,y]` and {math}`\mathbf{x}'= [x',y']` and the Green's function equals:
```{math}
G(\mathbf{x},\mathbf{x}') = \dfrac{\mathrm{e}^{-\mathrm{i}k_0||\mathbf{x}-\mathbf{x}'||}}{2\pi||\mathbf{x}-\mathbf{x}'||}.
```
The heading averaged geometrical radiation efficiency then equals:
```{math}
\bar{\sigma}(k_\mathrm{f}) = \dfrac{1}{2\pi}\int_0^{2\pi}\sigma(k_\mathrm{f},\phi)\mathrm{d}\phi = \mathrm{Re}\left\{\dfrac{\mathrm{i}k_0}{S}\int_S\int_SG(||\mathbf{x}-\mathbf{x}'||)\,\mathrm{J}_0(k_\mathrm{f}||\mathbf{x}-\mathbf{x}'||)\mathrm{d}\mathbf{x}\mathrm{d}\mathbf{x}'\right\}.
```
An order reduction of this radiation efficiency is achieved in [1] for rectangular geometries:
```{math}
\bar{\sigma}(k_\mathrm{f}) = \,&\mathrm{Re}\left\{\dfrac{2\mathrm{i}k_0}{\pi S}\int_0^{L_y}\left(\dfrac{L_xL_y\pi}{2}-L_yR-L_xR+\dfrac{R^2}{2}\right)\mathrm{J}_0(k_\mathrm{f}R)\mathrm{e}^{-\mathrm{i}k_0R}\mathrm{d}R\right\} \nonumber \\
&+\mathrm{Re}\left\{\dfrac{2\mathrm{i}k_0}{\pi S}\int_{L_y}^{L_x}\left(L_xL_y\arcsin\dfrac{L_y}{R}-\dfrac{L_y^2}{2}+L_x\sqrt{R^2-L_y^2}-L_xR\right)\mathrm{J}_0(k_\mathrm{f}R)\mathrm{e}^{-\mathrm{i}k_0R}\mathrm{d}R\right\} \nonumber \\
&+\mathrm{Re}\left\{\dfrac{2\mathrm{i}k_0}{\pi S}\int_{L_x}^{\sqrt{L_x^2+L_y^2}}\left[L_xL_y\left(\arcsin\dfrac{L_y}{R}-\arccos\dfrac{L_x}{R}\right)+L_y\sqrt{R^2-L_x^2}+L_x\sqrt{R^2-L_y^2}\right.\right. \nonumber \\
&\qquad\left.\left.-\dfrac{L_x^2+L_y^2+R^2}{2}\right]\mathrm{J}_0(k_\mathrm{f}R)\mathrm{e}^{-\mathrm{i}k_0R}\mathrm{d}R\right\}.
```
The aim of this text is to use this order reduction expression for polygonal geometries with right angles (e.g., an L shape), in order to obtain an efficient calculation of the radiation efficiency for these kind of structures.


## Radiation efficiency for non-rectangular geometries
The idea here is to split the non-rectangular geometry into rectangular subdomains. By dividing the geometry into rectangular subdomains with area {math}`S_i`, the radiation efficiency becomes:
```{math}
\bar{\sigma}(k_\mathrm{f}) &= \sum\limits_i\sum\limits_j \mathrm{Re}\left\{\dfrac{\mathrm{i}k_0}{S}\int_{S_i}\int_{S_j}M(k_\mathrm{f},||\mathbf{x}-\mathbf{x}'||)\mathrm{d}\mathbf{x}\mathrm{d}\mathbf{x}'\right\} \\
&= \sum\limits_i \mathrm{Re}\left\{\dfrac{\mathrm{i}k_0}{S}\int_{S_i}\int_{S_i}M(k_\mathrm{f},||\mathbf{x}-\mathbf{x}'||)\mathrm{d}\mathbf{x}\mathrm{d}\mathbf{x}'\right\}+\sum\limits_i\sum\limits_{j\neq i} \mathrm{Re}\left\{\dfrac{\mathrm{i}k_0}{S}\int_{S_i}\int_{S_j}M(k_\mathrm{f},||\mathbf{x}-\mathbf{x}'||)\mathrm{d}\mathbf{x}\mathrm{d}\mathbf{x}'\right\} \\
&= \sum\limits_i \mathrm{Re}\left\{\dfrac{\mathrm{i}k_0}{S}\int_{S_i}\int_{S_i}M(k_\mathrm{f},||\mathbf{x}-\mathbf{x}'||)\mathrm{d}\mathbf{x}\mathrm{d}\mathbf{x}'\right\}+2\sum\limits_i\sum\limits_{j>i} \mathrm{Re}\left\{\dfrac{\mathrm{i}k_0}{S}\int_{S_i}\int_{S_j}M(k_\mathrm{f},||\mathbf{x}-\mathbf{x}'||)\mathrm{d}\mathbf{x}\mathrm{d}\mathbf{x}'\right\},
```
where {math}`M(k_\mathrm{f},||\mathbf{x}-\mathbf{x}'||) = G(||\mathbf{x}-\mathbf{x}'||)\,\mathrm{J}_0(k_\mathrm{f}||\mathbf{x}-\mathbf{x}'||)`.



### Zero order approximation
As the Green's function and the zero order Bessel function in these equations decrease with distance, the correlation for distant points is very small and can be neglected. Therefore the second term could be neglected, as this concerns the correlation between points on a different subdomain:
```{math}
\bar{\sigma}^{(0)}(k_\mathrm{f}) = \sum\limits_i \mathrm{Re}\left\{\dfrac{\mathrm{i}k_0}{S}\int_{S_i}\int_{S_i}M(k_\mathrm{f},||\mathbf{x}-\mathbf{x}'||)\mathrm{d}\mathbf{x}\mathrm{d}\mathbf{x}'\right\} = \sum\limits_i \dfrac{S_i}{S}\bar{\sigma}_{i}(k_\mathrm{f}),
```
where
```{math}
\label{eq:radeffisol}
\bar{\sigma}_{i}(k_\mathrm{f}) = \mathrm{Re}\left\{\dfrac{\mathrm{i}k_0}{S_i}\int_{S_i}\int_{S_i}M(k_\mathrm{f},||\mathbf{x}-\mathbf{x}'||)\mathrm{d}\mathbf{x}\mathrm{d}\mathbf{x}'\right\}
```
is the radiation efficiency for subdomain {math}`i` in isolation.



### First order approximation
Alternatively, only adjacent subdomains {math}`j\in \epsilon_i` can be considered:
```{math}
\bar{\sigma}^{(1)}(k_\mathrm{f}) &= \sum\limits_i \mathrm{Re}\left\{\dfrac{\mathrm{i}k_0}{S}\int_{S_i}\int_{S_i}M(k_\mathrm{f},||\mathbf{x}-\mathbf{x}'||)\mathrm{d}\mathbf{x}\mathrm{d}\mathbf{x}'\right\}+2\sum\limits_i\sum\limits_{\substack{j\in \epsilon_i \\ j>i}} \mathrm{Re}\left\{\dfrac{\mathrm{i}k_0}{S}\int_{S_i}\int_{S_j}M(k_\mathrm{f},||\mathbf{x}-\mathbf{x}'||)\mathrm{d}\mathbf{x}\mathrm{d}\mathbf{x}'\right\} \\
\label{eq:radeffinterm}
&= \sum\limits_i \dfrac{S_i}{S}\bar{\sigma}_{i}(k_\mathrm{f})+2\sum\limits_i\sum\limits_{\substack{j\in \epsilon_i \\ j>i}} \mathrm{Re}\left\{\dfrac{\mathrm{i}k_0}{S}\int_{S_i}\int_{S_j}M(k_\mathrm{f},||\mathbf{x}-\mathbf{x}'||)\mathrm{d}\mathbf{x}\mathrm{d}\mathbf{x}'\right\},
```
Consider now two adjacent subdomains {math}`i` and {math}`j`. The total radiation efficiency equals:
```{math}
\bar{\sigma}_{ij}(k_\mathrm{f}) &= \mathrm{Re}\left\{\dfrac{\mathrm{i}k_0}{S_i+S_j}\int_{S_i+S_j}\int_{S_i+S_j}M(k_\mathrm{f},||\mathbf{x}-\mathbf{x}'||)\mathrm{d}\mathbf{x}\mathrm{d}\mathbf{x}'\right\} \\
&= \mathrm{Re}\left\{\dfrac{\mathrm{i}k_0}{S_i+S_j}\int_{S_i}\int_{S_i}M(k_\mathrm{f},||\mathbf{x}-\mathbf{x}'||)\mathrm{d}\mathbf{x}\mathrm{d}\mathbf{x}'\right\}+\mathrm{Re}\left\{\dfrac{\mathrm{i}k_0}{S_i+S_j}\int_{S_j}\int_{S_j}M(k_\mathrm{f},||\mathbf{x}-\mathbf{x}'||)\mathrm{d}\mathbf{x}\mathrm{d}\mathbf{x}'\right\} \nonumber \\
&\qquad+2\mathrm{Re}\left\{\dfrac{\mathrm{i}k_0}{S_i+S_j}\int_{S_i}\int_{S_j}M(k_\mathrm{f},||\mathbf{x}-\mathbf{x}'||)\mathrm{d}\mathbf{x}\mathrm{d}\mathbf{x}'\right\}
```
The first two integrals can be recognized:
```{math}
\bar{\sigma}_{ij}(k_\mathrm{f}) = \dfrac{S_i}{S_i+S_j}\bar{\sigma}_{i}(k_\mathrm{f})+\dfrac{S_j}{S_i+S_j}\bar{\sigma}_{j}(k_\mathrm{f})+2\mathrm{Re}\left\{\dfrac{\mathrm{i}k_0}{S_i+S_j}\int_{S_i}\int_{S_j}M(k_\mathrm{f},||\mathbf{x}-\mathbf{x}'||)\mathrm{d}\mathbf{x}\mathrm{d}\mathbf{x}'\right\},
```
and therefore:
```{math}
2\mathrm{Re}\left\{\dfrac{\mathrm{i}k_0}{S}\int_{S_i}\int_{S_j}M(k_\mathrm{f},||\mathbf{x}-\mathbf{x}'||)\mathrm{d}\mathbf{x}\mathrm{d}\mathbf{x}'\right\} &= 2\dfrac{S_i+S_j}{S}\mathrm{Re}\left\{\dfrac{\mathrm{i}k_0}{S_i+S_j}\int_{S_i}\int_{S_j}M(k_\mathrm{f},||\mathbf{x}-\mathbf{x}'||)\mathrm{d}\mathbf{x}\mathrm{d}\mathbf{x}'\right\} \nonumber \\
&= \dfrac{S_i+S_j}{S}\bar{\sigma}_{ij}(k_\mathrm{f})-\dfrac{S_i}{S}\bar{\sigma}_{i}(k_\mathrm{f})-\dfrac{S_j}{S}\bar{\sigma}_{j}(k_\mathrm{f})
```
The first order approximation of the radiation efficiency therefore becomes:
```{math}
\bar{\sigma}^{(1)}(k_\mathrm{f}) &= \sum\limits_i \dfrac{S_i}{S}\bar{\sigma}_{i}(k_\mathrm{f})+\sum\limits_i\sum\limits_{\substack{j\in \epsilon_i \\ j>i}} \left[\dfrac{S_i+S_j}{S}\bar{\sigma}_{ij}(k_\mathrm{f})-\dfrac{S_i}{S}\bar{\sigma}_{i}(k_\mathrm{f})-\dfrac{S_j}{S}\bar{\sigma}_{j}(k_\mathrm{f})\right].
```
This expression can be further simplified by noting that in the latter summation, the term {math}`S_i/S\bar{\sigma}_{i}(k_\mathrm{f})` appears {math}`N_\epsilon^i` times, with {math}`N_\epsilon^i` the number of adjacent elements {math}`\epsilon_i` for subdomain {math}`i`. Therefore, this expression becomes:
```{math}
\bar{\sigma}^{(1)}(k_\mathrm{f}) &= \sum\limits_i \dfrac{S_i}{S}\bar{\sigma}_{i}(k_\mathrm{f})+\sum\limits_i\sum\limits_{\substack{j\in \epsilon_i \\ j>i}} \dfrac{S_i+S_j}{S}\bar{\sigma}_{ij}(k_\mathrm{f})-\sum_iN_\epsilon^i\dfrac{S_i}{S}\bar{\sigma}_{i}(k_\mathrm{f}) \\
&= \sum\limits_i\sum\limits_{\substack{j\in \epsilon_i \\ j>i}} \dfrac{S_i+S_j}{S}\bar{\sigma}_{ij}(k_\mathrm{f})-\sum_i(N_\epsilon^i-1)\dfrac{S_i}{S}\bar{\sigma}_{i}(k_\mathrm{f})
```
It therefore follows that the radiation efficiency can be obtained from the radiation efficiencies {math}`\bar{\sigma}_{i}(k_\mathrm{f})` of the rectangular subdomains and the radiation efficiencies {math}`\bar{\sigma}_{ij}(k_\mathrm{f})` of two adjacent subdomains. These can be obtained efficiently using the efficient integral of Yu and Hopkins.



## References
[1] Y. Yu and C. Hopkins. Reduced order integration for the radiation eciency of a rectangular plate. JASA Express Letters, 1(6):062801, 2021.

