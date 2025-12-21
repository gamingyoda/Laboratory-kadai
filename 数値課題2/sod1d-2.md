# Sod問題（1D衝撃波管）に MUSCL法＋Limiter の導入

いま動いている **空間1次（界面の左右状態＝セル中心値そのまま）** のコードに対して、  
**MUSCL（空間2次）＋Limiter（衝撃波での振動抑制）** を入れる流れを、式と実装手順で整理する。

---

## 1. いまの空間1次がやっていること（復習）
有限体積法の更新は

$ \dfrac{dU_i}{dt} = -\dfrac{\hat F_{i+1/2}-\hat F_{i-1/2}}{\Delta x} $

空間1次（MUSCLなし）では界面 $i+1/2$ の左右状態を

- 左状態：$W_L = W_i$
- 右状態：$W_R = W_{i+1}$

と置いて、数値流束を

$ \hat F_{i+1/2} = \mathrm{RoeFlux}(W_i, W_{i+1}) $

で作っていた。  
これは実装は簡単だが、段差（衝撃波・接触面）が **にじみやすい**（解がぼやける）。

---

## 2. MUSCL法とは何か（何が変わる？）
MUSCLは「セル中心値だけ」ではなく、各セルの中で **直線（一次関数）** で状態を近似し、  
界面の左右状態を“少し賢く”作ってからリーマン問題（Roe）に渡す方法。

### 2.1 セル内を直線で近似する
セル $i$ の状態 $W_i$（例：$W=(\rho,u,p)$）に対して、傾き $\sigma_i$ を作り、

- セル右端（界面 $i+1/2$ の左側）  
  $ W_{i+1/2}^{L} = W_i + \dfrac{1}{2}\sigma_i $
- 次セル左端（界面 $i+1/2$ の右側）  
  $ W_{i+1/2}^{R} = W_{i+1} - \dfrac{1}{2}\sigma_{i+1} $

として、界面の左右状態を作る。

つまり **Roeに入れる状態が**
- 旧：$(W_i,\;W_{i+1})$
- 新：$(W_{i+1/2}^L,\;W_{i+1/2}^R)$  
に変わる。

---

## 3. でも“傾き”をそのまま使うと危険（Limiterが必要な理由）
衝撃波のような不連続の近くで、単純に2次精度の傾きを使うと、  
**ギブス現象のような振動（オーバーシュート／アンダーシュート）** が出やすい。

- 密度が一瞬マイナスっぽくなる
- 圧力が負になる
- 計算が破綻する

これを防ぐのが **Limiter（制限関数）**。  
Limiterは「傾きを作るけど、危ないときは傾きを弱める（0に近づける）」役。

---

## 4. Limiterの作り方（超重要：ここが“MUSCL＋Limiter”の本体）
セル $i$ の左右差分を

- 左差分：$ \Delta_- = W_i - W_{i-1} $
- 右差分：$ \Delta_+ = W_{i+1} - W_i $

として、傾き $\sigma_i$ を

$ \sigma_i = \mathrm{Limiter}(\Delta_-,\Delta_+) $

で決める。  
（実装では $\rho,u,p$ の各成分ごとに別々に limiter をかける。）

---

## 5. 代表的なLimiter（式と特徴）
以下は「2つの差分から傾きを返す」形で書く（成分ごとに適用）。

### 5.1 minmod（最も安全・拡散強め）
$ \mathrm{minmod}(a,b)=
\begin{cases}
0 & ab \le 0 \\
\mathrm{sign}(a)\min(|a|,|b|) & ab>0
\end{cases} $

- 符号が違う（極値っぽい）なら傾きを0にして振動を防ぐ
- その分、なめらか部分でも少しぼやけやすい

### 5.2 van Leer（なめらかさ優先・よく使う）
$ \mathrm{vanLeer}(a,b)=
\begin{cases}
0 & ab \le 0 \\
\dfrac{2ab}{a+b} & ab>0
\end{cases} $

- minmodより拡散が弱く、滑らかな部分がきれい

### 5.3 MC（Monotonized Central：バランス型）
$ \mathrm{MC}(a,b)=
\begin{cases}
0 & ab \le 0 \\
\mathrm{sign}(s)\min\left(|s|,2|a|,2|b|\right),\quad s=\dfrac{a+b}{2} & ab>0
\end{cases} $

- “きれいさ”と“安全”のバランスが良いことが多い

---

## 6. MUSCL＋Limiter を入れた「界面流束」の作り方（手順）
界面 $i+1/2$ の数値流束 $\hat F_{i+1/2}$ を作る手順はこうなる：

1. 全セルで primitive を作る：$W_i = (\rho_i,u_i,p_i)$
2. 全セルで limited slope を作る：$\sigma_i$
3. 界面状態を再構成：
   - $W_{i+1/2}^{L} = W_i + \dfrac{1}{2}\sigma_i$
   - $W_{i+1/2}^{R} = W_{i+1} - \dfrac{1}{2}\sigma_{i+1}$
4. Roe流束を計算：
   $ \hat F_{i+1/2} = \mathrm{RoeFlux}(W_{i+1/2}^{L}, W_{i+1/2}^{R}) $
5. 流束差分で rhs（$L(U)$）を作る：
   $ rhs[i] = -\dfrac{\hat F_{i+1/2}-\hat F_{i-1/2}}{\Delta x} $

---

## 7. 実装でどこが変わる？（コード改造ポイント）
あなたのコード構造が

- `zero_slope()`：境界（ゴーストセル）
- `calculate_dt()`：CFL
- `calculete_rhs()`：rhs作成
- `RK2`：時間積分

という形なら、**基本的に変えるのは `calculete_rhs()` の中身だけ**。

### 7.1 いま（空間1次）の rhs 計算
界面流束：
- $F[i] = \mathrm{RoeFlux}(W_i, W_{i+1})$

### 7.2 これから（MUSCL＋Limiter）の rhs 計算
界面流束：
- $F[i] = \mathrm{RoeFlux}(W_{i+1/2}^L, W_{i+1/2}^R)$

そのために `calculete_rhs()` の中で

- `W[i]`（primitive配列）
- `slope[i]`（傾き配列）
- `WL` `WR`（再構成した界面左右状態）

を作る処理が増える。

---

## 8. 実装の疑似コード（この形にすればOK）
（配列の添字範囲はあなたの実装に合わせて調整する。ここではゴースト2枚を想定。）

1) primitiveを作る
- for i=0..n_end-1:  W[i] = cons_to_prim(U[i])

2) slopeを作る（Limiter）
- for i=1..n_end-2:
  - dm = W[i] - W[i-1]
  - dp = W[i+1] - W[i]
  - slope[i] = limiter(dm, dp)   （成分ごと）

3) 界面流束を作る
- for i=1..n_end-3:   （界面 i = i+1/2）
  - WL = W[i]   + 0.5*slope[i]
  - WR = W[i+1] - 0.5*slope[i+1]
  - （もし WL/WR が負圧・負密度になるなら、局所的に 1次へ戻すのが安全）
  - F[i] = RoeFlux(WL, WR)

4) rhs を流束差分で作る
- for i=2..n_end-3:
  - rhs[i] = -(F[i] - F[i-1]) / dx

---

## 9. 「負圧・負密度」対策（MUSCL導入時は特に大事）
MUSCLで再構成した $W_{i+1/2}^L, W_{i+1/2}^R$ が

- $\rho \le 0$
- $p \le 0$

になると Roe の内部で音速が壊れて計算が破綻する。

よくある安全策：
- 再構成が非物理なら、その界面だけ **1次（セル中心値）に戻す**
  - $W_{i+1/2}^L \leftarrow W_i$
  - $W_{i+1/2}^R \leftarrow W_{i+1}$
- あるいは Roeがダメなら **Rusanov にフォールバック**

「まず動かす」なら、局所1次戻しが簡単で効果が高い。

---

## 10. まとめ（導入後に何が良くなる？）
- 空間1次：衝撃波・接触面が太くにじむ（安定だがぼやける）
- MUSCL＋Limiter：  
  - なめらかな部分は2次で精度アップ（勾配がきれい）
  - 不連続付近はLimiterが傾きを抑えて振動を防ぐ
  - 結果として **波面がシャープになりやすい**

次の作業は「`calculete_rhs()` を MUSCL版に置き換える」こと。  
時間積分（RK2）やCFL（dt計算）は基本そのままでOK。

---
