# CHƯƠNG 4 — ĐA THỨC NỘI SUY
## Mã giả các thuật toán

> Lưu ý: để render tốt trên GitHub/MathJax, các công thức dài nên viết ở dạng block `$$ ... $$` hoặc `\[ ... \]`, và các hệ phương trình nên dùng `aligned`/`cases` rõ ràng.

Tài liệu này được xây dựng theo nội dung Chương 4 của giáo trình. Các thuật toán chính gồm:

1. Nội suy Lagrange với mốc bất kỳ.
2. Đánh giá sai số của nội suy Lagrange.
3. Nội suy Lagrange với mốc đều.
4. Lập bảng tỉ sai phân.
5. Nội suy Newton.
6. Spline bậc ba tự nhiên.
7. Tính giá trị spline tại một điểm.

---

# 1. Nội suy Lagrange với mốc bất kỳ

Cho các mốc phân biệt

\[
x_0,x_1,\ldots,x_n
\]

và các giá trị

\[
y_i=f(x_i).
\]

Đa thức nội suy Lagrange là

\[
L_n(x)
=
\sum_{i=0}^n y_i\omega_i(x),
\]

trong đó

\[
\omega_i(x)
=
\prod_{\substack{j=0\\j\ne i}}^n
\frac{x-x_j}{x_i-x_j}.
\]

## Mã giả

```text
THUẬT TOÁN LAGRANGE(x_nodes, y_nodes, x)

Đầu vào:
    x_nodes = [x_0, x_1, ..., x_n]
    y_nodes = [y_0, y_1, ..., y_n]
    x       = điểm cần nội suy

L = 0

Với i = 0,1,...,n:
    omega = 1

    Với j = 0,1,...,n:
        Nếu j khác i:
            omega = omega * (x - x_j)/(x_i - x_j)

    L = L + y_i * omega

Trả về L
```

---

# 2. Xây dựng đa thức Lagrange

Nếu cần biểu thức đa thức thay vì chỉ giá trị tại một điểm:

```text
THUẬT TOÁN ĐA_THỨC_LAGRANGE(x_nodes, y_nodes)

Khởi tạo L(x) = 0

Với i = 0,1,...,n:
    omega_i(x) = 1

    Với j = 0,1,...,n:
        Nếu j khác i:
            omega_i(x)
            = omega_i(x) * (x - x_j)/(x_i - x_j)

    L(x) = L(x) + y_i * omega_i(x)

Khai triển và rút gọn L(x)

Trả về L(x)
```

---

# 3. Đánh giá sai số nội suy Lagrange

Nếu \(f\) khả vi liên tục đến cấp \(n+1\) trên \([a,b]\), thì tại điểm \(\bar x\)

\[
|f(\bar x)-L_n(\bar x)|
\le
\frac{M}{(n+1)!}
\left|
\prod_{i=0}^n(\bar x-x_i)
\right|,
\]

với

\[
M=
\max_{a\le x\le b}|f^{(n+1)}(x)|.
\]

## Mã giả

```text
THUẬT TOÁN SAI_SỐ_LAGRANGE(x_nodes, x_bar, M)

n = số mốc - 1

psi = 1

Với i = 0,1,...,n:
    psi = psi * (x_bar - x_i)

Err = M/factorial(n+1) * |psi|

Trả về Err
```

---

# 4. Lagrange với mốc đều

Giả sử

\[
x_i=a+ih,
\qquad
h=\frac{b-a}{n}.
\]

Đặt

\[
x=a+ht
\quad\Longrightarrow\quad
t=\frac{x-a}{h}.
\]

Ta có

\[
L_n(x)=P_n(t)
=
\sum_{i=0}^n y_i\ell_i(t).
\]

Trong tính toán, có thể dùng trực tiếp công thức Lagrange trên biến \(t\).

## Mã giả

```text
THUẬT TOÁN LAGRANGE_MỐC_ĐỀU(a, h, y_nodes, x)

t = (x-a)/h
n = số phần tử của y_nodes - 1

P = 0

Với i = 0,1,...,n:
    ell = 1

    Với j = 0,1,...,n:
        Nếu j khác i:
            ell = ell * (t-j)/(i-j)

    P = P + y_i*ell

Trả về P
```

---

# 5. Lập bảng tỉ sai phân

Quy ước

\[
D_i^{(0)}=y_i.
\]

Với \(k=1,2,\ldots,n\),

\[
D_i^{(k)}
=
\frac{
D_{i+1}^{(k-1)}-D_i^{(k-1)}
}{
x_{i+k}-x_i
}.
\]

## Mã giả

```text
THUẬT TOÁN BẢNG_TỈ_SAI_PHÂN(x_nodes, y_nodes)

n = số mốc - 1

Tạo bảng D kích thước (n+1) x (n+1)
Gán D[i,0] = y_i

Với k = 1,2,...,n:
    Với i = 0,1,...,n-k:
        D[i,k]
        =
        (D[i+1,k-1] - D[i,k-1])
        /
        (x[i+k] - x[i])

Trả về D
```

Các hệ số của đa thức Newton là

\[
a_k=D_0^{(k)}.
\]

---

# 6. Nội suy Newton

Đa thức Newton:

\[
L_n(x)
=
a_0
+
\sum_{k=1}^{n}
a_k
\prod_{j=0}^{k-1}(x-x_j),
\]

với

\[
a_k=f[x_0;x_1;\ldots;x_k].
\]

## Mã giả

```text
THUẬT TOÁN NEWTON_INTERPOLATION(x_nodes, y_nodes, x)

D = BẢNG_TỈ_SAI_PHÂN(x_nodes, y_nodes)

L = D[0,0]
P = 1

Với k = 1,2,...,n:
    P = P * (x - x[k-1])
    L = L + D[0,k] * P

Trả về L
```

---

# 7. Xây dựng đa thức Newton

```text
THUẬT TOÁN ĐA_THỨC_NEWTON(x_nodes, y_nodes)

D = BẢNG_TỈ_SAI_PHÂN(x_nodes, y_nodes)

L(x) = D[0,0]
P(x) = 1

Với k = 1,2,...,n:
    P(x) = P(x) * (x - x[k-1])
    L(x) = L(x) + D[0,k] * P(x)

Khai triển và rút gọn L(x)

Trả về L(x)
```

---

# 8. Spline bậc ba tự nhiên

Cho

\[
x_0<x_1<\cdots<x_n,
\qquad
y_i=f(x_i).
\]

Đặt

\[
h_i=x_{i+1}-x_i.
\]

Với spline tự nhiên:

\[
m_0=m_n=0.
\]

Các giá trị

\[
m_i=S''(x_i)
\]

được tìm bằng hệ

\[
\frac{h_i}{6}m_i
+
\frac{h_i+h_{i+1}}{3}m_{i+1}
+
\frac{h_{i+1}}{6}m_{i+2}
=
\frac{y_{i+2}-y_{i+1}}{h_{i+1}}
-
\frac{y_{i+1}-y_i}{h_i}.
\]

## Mã giả

```text
THUẬT TOÁN SPLINE_TỰ_NHIÊN(x_nodes, y_nodes)

n = số mốc - 1

Bước 1:
    Với i = 0,...,n-1:
        h_i = x_{i+1} - x_i

Bước 2:
    Đặt m_0 = 0
    Đặt m_n = 0

    Xây dựng hệ tuyến tính cho
    m_1, m_2, ..., m_{n-1}

    Với r = 1,...,n-1:
        hệ số trái   = h_{r-1}/6
        hệ số giữa  = (h_{r-1}+h_r)/3
        hệ số phải  = h_r/6

        vế phải =
            (y_{r+1}-y_r)/h_r
            -
            (y_r-y_{r-1})/h_{r-1}

    Giải hệ để tìm m_1,...,m_{n-1}

Bước 3:
    Với i = 0,...,n-1:
        A_i = y_i - m_i*h_i^2/6
        B_i = y_{i+1} - m_{i+1}*h_i^2/6

Bước 4:
    Trên [x_i, x_{i+1}]:

    S_i(x)
    =
    m_{i+1}(x-x_i)^3/(6h_i)
    +
    m_i(x_{i+1}-x)^3/(6h_i)
    +
    A_i(x_{i+1}-x)/h_i
    +
    B_i(x-x_i)/h_i

Trả về toàn bộ các đa thức S_i(x)
```

---

# 9. Tính giá trị spline tại một điểm

```text
THUẬT TOÁN GIÁ_TRỊ_SPLINE(x_nodes, spline, x)

Tìm i sao cho:
    x_i <= x <= x_{i+1}

Tính:
    S(x) = S_i(x)

Trả về S(x)
```

---

# 10. Cấu trúc hàm Python đề nghị

Notebook của chương sử dụng các hàm:

```text
lagrange_value(...)
lagrange_polynomial(...)
lagrange_error_bound(...)
divided_difference_table(...)
newton_value(...)
newton_polynomial(...)
natural_cubic_spline(...)
spline_value(...)
```

Mỗi hàm được thiết kế để có thể tái sử dụng cho các ví dụ và bài tập của chương.
