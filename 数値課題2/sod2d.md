# Sod問題ソルバを2次元に拡張するための解説（Roe-FDS + TVD RK2, 構造格子）

> 目的：いまの **1D（$x$方向のみ）** の衝撃波管ソルバを、**2D（$x,y$）** の Euler 方程式へ拡張する。  
> 手法はそのまま：**有限体積法（FVM） + Roe-FDS（数値流束） + TVD RK2（時間積分）**。

---

## 1. 2次元 Euler 方程式（保存形）

2D では保存量が4成分になる：

$$
U=
\begin{pmatrix}
\rho\\
\rho u\\
\rho v\\
E
\end{pmatrix}
$$

- $\rho$：密度  
- $u$：$x$方向速度  
- $v$：$y$方向速度  
- $E$：全エネルギー密度

圧力は（理想気体）：

$$
p=(\gamma-1)\left(E-\frac{1}{2}\rho(u^2+v^2)\right)
$$

支配方程式は「$x$方向の流束 $F$」と「$y$方向の流束 $G$」で

$$
\frac{\partial U}{\partial t}
+\frac{\partial F(U)}{\partial x}
+\frac{\partial G(U)}{\partial y}
=0
$$

---

## 2. 2Dの物理流束 $F(U),G(U)$

### 2.1 $x$方向流束（縦面を通る）
$$
F(U)=
\begin{pmatrix}
\rho u\\
\rho u^2 + p\\
\rho u v\\
u(E+p)
\end{pmatrix}
$$

### 2.2 $y$方向流束（横面を通る）
$$
G(U)=
\begin{pmatrix}
\rho v\\
\rho u v\\
\rho v^2 + p\\
v(E+p)
\end{pmatrix}
$$

---

## 3. 有限体積法（2D版の更新式）

セルを $(i,j)$ とする（$i$が$x$、$j$が$y$）。

セルの更新は「4つの面からの出入り」で決まる：

$$
\frac{dU_{i,j}}{dt}
=
-\frac{\hat F_{i+1/2,j}-\hat F_{i-1/2,j}}{\Delta x}
-\frac{\hat G_{i,j+1/2}-\hat G_{i,j-1/2}}{\Delta y}
$$

ここで必要なのはセル中心の $F(U_{i,j})$ ではなく、
**界面を通る数値流束** $\hat F,\hat G$。

---

## 4. 「1DのRoe-FDS」を2Dにどう使うか（超重要）

2D でも、界面は「法線方向」を持つので、
**界面に垂直な方向の1D問題**として Roe を使える。

構造格子（直交格子）の場合：

- $x$方向界面（縦面）：法線は $(1,0)$ → **法線速度は $u$**
- $y$方向界面（横面）：法線は $(0,1)$ → **法線速度は $v$**

つまり

- $\hat F_{i+1/2,j}$ は **$x$方向Roe流束**
- $\hat G_{i,j+1/2}$ は **$y$方向Roe流束**

をそれぞれ計算すればよい。

---

## 5. 2D版 Roe-FDS（軸方向界面の式だけ覚えればOK）

2D Euler の $x$方向（縦面）では、波は4本：

- 音波（左向き） $u-a$
- 接触波（エントロピー波） $u$
- せん断波（横速度の段差が運ばれる） $u$
- 音波（右向き） $u+a$

つまり固有値は

$$
\lambda_1=\tilde u-\tilde a,\quad
\lambda_2=\tilde u,\quad
\lambda_3=\tilde u,\quad
\lambda_4=\tilde u+\tilde a
$$

ここで $\tilde{\cdot}$ は Roe 平均。

---

### 5.1 Roe平均（2D版）

左右状態（primitive）を $W_L=(\rho_L,u_L,v_L,p_L)$、$W_R=(\rho_R,u_R,v_R,p_R)$。

全エンタルピー：

$$
H=\frac{E+p}{\rho}
$$

Roe平均：

$$
\tilde u=\frac{\sqrt{\rho_L}u_L+\sqrt{\rho_R}u_R}{\sqrt{\rho_L}+\sqrt{\rho_R}},\quad
\tilde v=\frac{\sqrt{\rho_L}v_L+\sqrt{\rho_R}v_R}{\sqrt{\rho_L}+\sqrt{\rho_R}}
$$

$$
\tilde H=\frac{\sqrt{\rho_L}H_L+\sqrt{\rho_R}H_R}{\sqrt{\rho_L}+\sqrt{\rho_R}}
$$

音速：

$$
\tilde a^2=(\gamma-1)\left(\tilde H-\frac{1}{2}(\tilde u^2+\tilde v^2)\right)
$$

---

### 5.2 波の強さ（$x$方向界面）

ジャンプ：

$$
\Delta\rho=\rho_R-\rho_L,\quad
\Delta u=u_R-u_L,\quad
\Delta v=v_R-v_L,\quad
\Delta p=p_R-p_L
$$

代表密度（よく使う形）：

$$
\tilde\rho=\sqrt{\rho_L\rho_R}
$$

係数（$x$方向界面）：

- 接触（エントロピー）成分：
$$
\alpha_2=\Delta\rho-\frac{\Delta p}{\tilde a^2}
$$

- せん断（横速度）成分：
$$
\alpha_3=\tilde\rho\,\Delta v
$$

- 音波（左右）成分：
$$
\alpha_1=\frac{\Delta p-\tilde\rho\,\tilde a\,\Delta u}{2\tilde a^2},\quad
\alpha_4=\frac{\Delta p+\tilde\rho\,\tilde a\,\Delta u}{2\tilde a^2}
$$

---

### 5.3 右固有ベクトル（$x$方向界面）

（保存変数 $U=[\rho,\rho u,\rho v,E]^T$ に対する形）

$$
r_1=
\begin{pmatrix}
1\\
\tilde u-\tilde a\\
\tilde v\\
\tilde H-\tilde u\tilde a
\end{pmatrix},\quad
r_2=
\begin{pmatrix}
1\\
\tilde u\\
\tilde v\\
\frac12(\tilde u^2+\tilde v^2)
\end{pmatrix}
$$

$$
r_3=
\begin{pmatrix}
0\\
0\\
1\\
\tilde v
\end{pmatrix},\quad
r_4=
\begin{pmatrix}
1\\
\tilde u+\tilde a\\
\tilde v\\
\tilde H+\tilde u\tilde a
\end{pmatrix}
$$

---

### 5.4 Roe-FDS 数値流束（$x$方向）

$$
\hat F
=
\frac{F_L+F_R}{2}
-\frac12\sum_{k=1}^{4}|\lambda_k|\alpha_k r_k
$$

entropy fix を入れるなら $|\lambda_k|$ を小さい領域で滑らかにする。

---

### 5.5 $y$方向界面はどうする？

$y$方向界面は「$x$と$y$を入れ替え」ればOK。

- 法線速度：$v$
- 横方向（接線）速度：$u$
- せん断波は「$u$ の段差」が運ばれる

固有値は

$$
\mu_1=\tilde v-\tilde a,\quad
\mu_2=\tilde v,\quad
\mu_3=\tilde v,\quad
\mu_4=\tilde v+\tilde a
$$

係数は

$$
\beta_2=\Delta\rho-\frac{\Delta p}{\tilde a^2}
$$

$$
\beta_3=\tilde\rho\,\Delta u
$$

$$
\beta_1=\frac{\Delta p-\tilde\rho\,\tilde a\,\Delta v}{2\tilde a^2},\quad
\beta_4=\frac{\Delta p+\tilde\rho\,\tilde a\,\Delta v}{2\tilde a^2}
$$

（右固有ベクトルも同様に “$u \leftrightarrow v$、$F\leftrightarrow G$” の対応で作る。）

---

## 6. 2Dの rhs（$L(U)$）の作り方（実装の骨格）

1. 全セルで $U_{i,j}\to W_{i,j}$（primitive）に変換  
2. すべての $x$界面で $\hat F_{i+1/2,j}$ を Roe で計算  
3. すべての $y$界面で $\hat G_{i,j+1/2}$ を Roe で計算  
4. 各セルで

$$
rhs_{i,j}
=
-\frac{\hat F_{i+1/2,j}-\hat F_{i-1/2,j}}{\Delta x}
-\frac{\hat G_{i,j+1/2}-\hat G_{i,j-1/2}}{\Delta y}
$$

---

## 7. CFL条件（2D版の $\Delta t$）

2D の陽解法では「$x$と$y$の両方向に波が進む」ので、典型的に

$$
\Delta t
=
\mathrm{CFL}\;\Bigg/\;
\max_{i,j}\left(
\frac{|u_{i,j}|+a_{i,j}}{\Delta x}
+
\frac{|v_{i,j}|+a_{i,j}}{\Delta y}
\right)
$$

（簡単にしたいなら $\Delta t=\mathrm{CFL}\cdot\min(\Delta x,\Delta y)/\max(|u|+|v|+a)$ でも動くが、上式の方が素直。）

---

## 8. 時間積分（TVD RK2）はそのまま 2D へ

空間離散で $rhs=L(U)$ が作れれば 1D と同じ：

1段目：
$$
U^{(1)}=U^n+\Delta t\,L(U^n)
$$

2段目：
$$
U^{n+1}
=
\frac12U^n+\frac12\left(U^{(1)}+\Delta t\,L(U^{(1)})\right)
$$

---

## 9. 境界条件（ゴーストセル）は2Dに拡張

1Dで「左右2セルコピー」していたのを、2Dでは

- 左端：$U_{0,j}=U_{2,j}$、$U_{1,j}=U_{2,j}$
- 右端：$U_{N_x+2ng-1,j}=U_{N_x+2ng-3,j}$ など
- 下端：$U_{i,0}=U_{i,2}$ …
- 上端：$U_{i,N_y+2ng-1}=U_{i,N_y+2ng-3}$ …

のように **四辺**で行う。

---

## 10. 2D化で「コード上」なにが増えるか

- 変数が4成分になる：$(\rho,\rho u,\rho v,E)$
- 配列が 2D（実際は1Dに潰して保持）になる
- rhs計算が2方向の流束差分になる
- Roe flux が
  - $x$方向用（$\hat F$）
  - $y$方向用（$\hat G$）
  の2種類必要
- CFL計算が2Dになる
- 境界条件が四辺になる
- 出力は $(x,y,\rho,u,v,p)$（必要ならスライスやカラーマップ表示）

---
