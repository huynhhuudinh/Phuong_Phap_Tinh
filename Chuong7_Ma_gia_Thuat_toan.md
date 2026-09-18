# CHƯƠNG 7 — PHƯƠNG TRÌNH VI PHÂN
## Mã giả các thuật toán

> Lưu ý: để render tốt trên GitHub/MathJax, các công thức dài nên viết ở dạng block `$$ ... $$` hoặc `$$ ... $$`, và các hệ phương trình nên dùng `aligned`/`cases` rõ ràng.

Tài liệu này được xây dựng theo nội dung Chương 7 của giáo trình.

Các thuật toán chính gồm:

1. Phương pháp Euler.
2. Phương pháp Euler cải tiến.
3. Phương pháp Runge--Kutta bậc bốn cho phương trình vô hướng.
4. Chuyển phương trình vi phân bậc cao về hệ bậc nhất.
5. Phương pháp Runge--Kutta bậc bốn cho hệ phương trình vi phân.
6. Áp dụng RK4 cho mô hình dịch bệnh SIS.
7. Tính sai số khi biết nghiệm chính xác.

---

# 1. Phương pháp Euler

Xét bài toán Cauchy

$$
\begin{aligned}
y'=f(x;y),\\
y(a)=y_0.
\end{aligned}
$$

Chia đoạn $[a,b]$ thành $n$ đoạn bằng nhau:

$$
h=\frac{b-a}{n},
\qquad
x_i=a+ih.
$$

Công thức Euler:

$$
y_{i+1}
=
y_i+h f(x_i;y_i).
$$

## Mã giả

```text
THUẬT TOÁN EULER(f, a, b, y0, n)

h = (b-a)/n
x = a
y = y0

Ghi lại (x,y)

Với i = 0,...,n-1:
    y = y + h*f(x,y)
    x = x + h

    Ghi lại (x,y)

Trả về toàn bộ các cặp (x_i,y_i)
```

---

# 2. Phương pháp Euler cải tiến

Tại mỗi bước, trước tiên tính giá trị dự báo

$$
\widetilde y_{i+1}
=
y_i+h f(x_i;y_i),
$$

sau đó hiệu chỉnh

$$
y_{i+1}
=
y_i
+
\frac{h}{2}
\left[
f(x_i;y_i)
+
f(x_{i+1};\widetilde y_{i+1})
\right].
$$

## Mã giả

```text
THUẬT TOÁN EULER_CẢI_TIẾN(f, a, b, y0, n)

h = (b-a)/n
x = a
y = y0

Ghi lại (x,y)

Với i = 0,...,n-1:

    y_predict =
        y + h*f(x,y)

    x_next =
        x + h

    y_next =
        y
        + h/2 *
          [f(x,y)
           + f(x_next,y_predict)]

    x = x_next
    y = y_next

    Ghi lại:
        x, y_predict, y

Trả về bảng kết quả
```

---

# 3. Runge--Kutta bậc bốn cho phương trình vô hướng

Với mỗi bước,

$$
k_1
=
h f(x_i;y_i),
$$

$$
k_2
=
h f\left(
x_i+\frac h2;
y_i+\frac{k_1}{2}
\right),
$$

$$
k_3
=
h f\left(
x_i+\frac h2;
y_i+\frac{k_2}{2}
\right),
$$

$$
k_4
=
h f(x_i+h;y_i+k_3).
$$

Cập nhật

$$
y_{i+1}
=
y_i
+
\frac16
(k_1+2k_2+2k_3+k_4).
$$

## Mã giả

```text
THUẬT TOÁN RK4(f, a, b, y0, n)

h = (b-a)/n
x = a
y = y0

Ghi lại (x,y)

Với i = 0,...,n-1:

    k1 = h*f(x,y)

    k2 =
        h*f(
            x+h/2,
            y+k1/2
        )

    k3 =
        h*f(
            x+h/2,
            y+k2/2
        )

    k4 =
        h*f(
            x+h,
            y+k3
        )

    y =
        y
        + (k1+2*k2+2*k3+k4)/6

    x = x+h

    Ghi lại:
        x, y, k1, k2, k3, k4

Trả về bảng kết quả
```

---

# 4. Chuyển phương trình bậc cao về hệ bậc nhất

Xét

$$
y^{(N)}
=
f(x;y,y',\ldots,y^{(N-1)}).
$$

Đặt

$$
y_1=y,\quad
y_2=y',\quad
\ldots,\quad
y_N=y^{(N-1)}.
$$

Khi đó,

$$
\begin{aligned}
y_1' &= y_2, \\
 y_2' &= y_3, \\
 &\vdots \\
 y_{N-1}' &= y_N, \\
 y_N' &= f(x;y_1,\ldots,y_N).
\end{aligned}
$$

## Mã giả

```text
THUẬT TOÁN CHUYỂN_BẬC_CAO_THÀNH_HỆ(N, f)

Định nghĩa vector trạng thái:
    Y = (y1,y2,...,yN)

Định nghĩa hàm vector F:

    F_1 = y2
    F_2 = y3
    ...
    F_(N-1) = yN
    F_N = f(x,y1,y2,...,yN)

Khi đó:
    Y' = F(x,Y)

Trả về hệ bậc nhất tương đương
```

---

# 5. Runge--Kutta bậc bốn cho hệ phương trình vi phân

Xét

$$
\mathbf Y'
=
\mathbf f(x;\mathbf Y).
$$

Với $\mathbf Y_i\approx\mathbf Y(x_i)$,

$$
\mathbf k_1
=
h\mathbf f(x_i;\mathbf Y_i),
$$

$$
\mathbf k_2
=
h\mathbf f
\left(
x_i+\frac h2;
\mathbf Y_i+\frac{\mathbf k_1}{2}
\right),
$$

$$
\mathbf k_3
=
h\mathbf f
\left(
x_i+\frac h2;
\mathbf Y_i+\frac{\mathbf k_2}{2}
\right),
$$

$$
\mathbf k_4
=
h\mathbf f
\left(
x_i+h;
\mathbf Y_i+\mathbf k_3
\right).
$$

Cập nhật

$$
\mathbf Y_{i+1}
=
\mathbf Y_i
+
\frac16
\left(
\mathbf k_1
+
2\mathbf k_2
+
2\mathbf k_3
+
\mathbf k_4
\right).
$$

## Mã giả

```text
THUẬT TOÁN RK4_HỆ(F, a, b, Y0, n)

h = (b-a)/n
x = a
Y = Y0

Ghi lại (x,Y)

Với i = 0,...,n-1:

    K1 = h*F(x,Y)

    K2 =
        h*F(
            x+h/2,
            Y+K1/2
        )

    K3 =
        h*F(
            x+h/2,
            Y+K2/2
        )

    K4 =
        h*F(
            x+h,
            Y+K3
        )

    Y =
        Y
        + (K1+2*K2+2*K3+K4)/6

    x = x+h

    Ghi lại (x,Y)

Trả về bảng kết quả
```

---

# 6. Mô hình SIS

Mô hình:

$$
\begin{aligned}
S' &= -\frac{\beta}{N}SI + \gamma I, \\
I' &= \frac{\beta}{N}SI - \gamma I.
\end{aligned}
$$

trong đó

$$
N=S+I.
$$

Hệ số lây nhiễm cơ bản:

$$
\mathcal R_0
=
\frac{\beta}{\gamma}.
$$

## Mã giả

```text
THUẬT TOÁN SIS_RK4(beta, gamma, S0, I0, a, b, h)

N = S0 + I0

Định nghĩa F(t,Y):

    S = Y[1]
    I = Y[2]

    dS =
        -(beta/N)*S*I
        + gamma*I

    dI =
        (beta/N)*S*I
        - gamma*I

    Trả về (dS,dI)

n = (b-a)/h

Áp dụng:
    RK4_HỆ(F,a,b,(S0,I0),n)

Tính:
    R0 = beta/gamma

Trả về bảng S_i, I_i và R0
```

Nếu

$$
\mathcal R_0\le1,
$$

chương mô tả $I(t)$ tiến dần về $0$.

Nếu

$$
\mathcal R_0>1,\qquad I_0>0,
$$

mức cân bằng được nêu trong chương là

$$
I_\infty
=
\left(
1-\frac{\gamma}{\beta}
\right)N.
$$

---

# 7. Tính sai số khi biết nghiệm chính xác

Nếu biết nghiệm đúng $y(x)$, tại mốc $x_i$,

$$
\Delta_{y_i}
=
|y(x_i)-y_i|.
$$

## Mã giả

```text
THUẬT TOÁN TÍNH_SAI_SỐ(x_values, y_approx, y_exact)

Với i = 0,...,n:
    exact_i = y_exact(x_values[i])
    error_i =
        |exact_i - y_approx[i]|

Trả về:
    exact_i, error_i
```

Để so sánh các phương pháp trên thang logarit:

```text
log_error_i = log10(error_i)
```

chỉ thực hiện với các sai số dương.

---

# 8. Cấu trúc các hàm Python trong notebook

```text
euler(...)
improved_euler(...)
rk4(...)
rk4_system(...)
attach_exact_solution(...)
sis_rhs(...)
run_sis(...)
```

Các hàm trả về bảng lịch sử lặp nhằm dễ dàng đối chiếu với các bảng số trong giáo trình.
