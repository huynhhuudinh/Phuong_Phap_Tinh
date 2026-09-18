# CHƯƠNG 5 — PHƯƠNG PHÁP BÌNH PHƯƠNG BÉ NHẤT
## Mã giả các thuật toán

> Lưu ý: để render tốt trên GitHub/MathJax, các công thức dài nên viết ở dạng block `$$ ... $$` hoặc `$$ ... $$`, và các hệ phương trình nên dùng `aligned`/`cases` rõ ràng.

Tài liệu này được xây dựng theo nội dung Chương 5 của giáo trình. Các thuật toán chính gồm:

1. Xấp xỉ đường thẳng $y=ax+b$.
2. Tính hệ số xác định $R^2$.
3. Xấp xỉ tổ hợp tuyến tính của hai hàm cơ sở $y=ag(x)+bh(x)$.
4. Xấp xỉ đa thức bậc không quá $m$.
5. Xấp xỉ hàm mũ $y=ae^{bx}$ sau phép biến đổi logarit.
6. Xấp xỉ hàm lũy thừa $y=ax^b$ sau phép biến đổi logarit.

---

# 1. Xấp xỉ đường thẳng $y=ax+b$

Với dữ liệu
$$
(x_i,y_i),\qquad i=1,2,\ldots,n,
$$
ta tìm $a,b$ sao cho
$$
S(a,b)=\sum_{i=1}^n(ax_i+b-y_i)^2
$$
nhỏ nhất.

Hệ phương trình chuẩn là
$$
\begin{aligned}
\left(\sum x_i^2\right)a + \left(\sum x_i\right)b &= \sum x_i y_i, \\
\left(\sum x_i\right)a + nb &= \sum y_i.
\end{aligned}
$$

Đặt
$$
S_1=\sum x_i,\quad
S_2=\sum x_i^2,\quad
S_3=\sum y_i,\quad
S_4=\sum x_i y_i.
$$

Khi
$$
D=nS_2-S_1^2\ne0,
$$
ta có
$$
a=\frac{nS_4-S_1S_3}{D},
\qquad
b=\frac{S_2S_3-S_1S_4}{D}.
$$

### Mã giả

```text
THUẬT TOÁN BÌNH_PHƯƠNG_BÉ_NHẤT_ĐƯỜNG_THẲNG(x, y)

Đầu vào:
    x_i, y_i, i = 1,...,n

Khởi tạo:
    S1 = 0
    S2 = 0
    S3 = 0
    S4 = 0

Với i = 1,...,n:
    S1 = S1 + x_i
    S2 = S2 + x_i^2
    S3 = S3 + y_i
    S4 = S4 + x_i*y_i

D = n*S2 - S1^2

Nếu D = 0:
    Trả về "Không xác định được đường thẳng xấp xỉ duy nhất"

a = (n*S4 - S1*S3)/D
b = (S2*S3 - S1*S4)/D

Trả về a, b
```

---

# 2. Hệ số xác định $R^2$

Với các giá trị dự đoán
$$
\widehat y_i=f(x_i)
$$
và
$$
\overline y=\frac1n\sum_{i=1}^n y_i,
$$
hệ số xác định được tính bởi
$$
R^2
=
1-
\frac{\sum_{i=1}^n(y_i-\widehat y_i)^2}
{\sum_{i=1}^n(y_i-\overline y)^2}.
$$

### Mã giả

```text
THUẬT TOÁN HỆ_SỐ_XÁC_ĐỊNH(y, y_hat)

Đầu vào:
    y_i      : giá trị quan sát
    y_hat_i  : giá trị dự đoán

y_bar = trung bình của y_i

SSE = 0
SST = 0

Với i = 1,...,n:
    SSE = SSE + (y_i - y_hat_i)^2
    SST = SST + (y_i - y_bar)^2

Nếu SST = 0:
    Trả về "R^2 không xác định"

R2 = 1 - SSE/SST

Trả về R2
```

---

# 3. Xấp xỉ bằng hai hàm cơ sở

Xét mô hình
$$
y=ag(x)+bh(x).
$$

Tổng bình phương sai số là
$$
S(a,b)
=
\sum_{i=1}^n
\left[
ag(x_i)+bh(x_i)-y_i
\right]^2.
$$

Đặt
$$
A_{11}=\sum g^2(x_i),\quad
A_{12}=\sum g(x_i)h(x_i),\quad
A_{22}=\sum h^2(x_i),
$$
$$
B_1=\sum g(x_i)y_i,\qquad
B_2=\sum h(x_i)y_i.
$$

Hệ phương trình chuẩn là
$$
\begin{aligned}
A_{11}a+A_{12}b=B_1,\\
A_{12}a+A_{22}b=B_2.
\end{aligned}
$$

### Mã giả

```text
THUẬT TOÁN BÌNH_PHƯƠNG_BÉ_NHẤT_HAI_HÀM_CƠ_SỞ(x, y, g, h)

Khởi tạo:
    A11 = 0
    A12 = 0
    A22 = 0
    B1 = 0
    B2 = 0

Với i = 1,...,n:
    gi = g(x_i)
    hi = h(x_i)

    A11 = A11 + gi^2
    A12 = A12 + gi*hi
    A22 = A22 + hi^2
    B1 = B1 + gi*y_i
    B2 = B2 + hi*y_i

D = A11*A22 - A12^2

Nếu D = 0:
    Trả về "Không xác định duy nhất a và b"

a = (B1*A22 - A12*B2)/D
b = (A11*B2 - A12*B1)/D

Trả về a, b
```

---

# 4. Xấp xỉ đa thức bậc không quá $m$

Xét mô hình
$$
y=a_0+a_1x+\cdots+a_mx^m.
$$

Ta đặt
$$
S_k=\sum_{i=1}^n x_i^k,\qquad k=0,1,\ldots,2m,
$$
và
$$
T_k=\sum_{i=1}^n x_i^k y_i,\qquad k=0,1,\ldots,m.
$$

Hệ phương trình chuẩn có dạng
$$
\mathbf A\mathbf a=\mathbf b,
$$
trong đó
$$
A_{jk}=S_{j+k},
\qquad
b_j=T_j,
\qquad
j,k=0,1,\ldots,m.
$$

### Mã giả

```text
THUẬT TOÁN BÌNH_PHƯƠNG_BÉ_NHẤT_ĐA_THỨC(x, y, m)

Với k = 0,...,2m:
    S[k] = 0

Với k = 0,...,m:
    T[k] = 0

Với i = 1,...,n:
    Với k = 0,...,2m:
        S[k] = S[k] + x_i^k

    Với k = 0,...,m:
        T[k] = T[k] + x_i^k*y_i

Với j = 0,...,m:
    Với k = 0,...,m:
        A[j,k] = S[j+k]

    b[j] = T[j]

Nếu det(A) = 0:
    Trả về "Không xác định được đa thức xấp xỉ duy nhất"

Giải hệ:
    A*a = b

Trả về:
    a_0, a_1, ..., a_m
```

---

# 5. Xấp xỉ hàm mũ $y=ae^{bx}$

Giả sử
$$
a>0,\qquad y_i>0.
$$

Lấy logarit tự nhiên:
$$
Y_i=\ln y_i.
$$

Đặt
$$
A=b,\qquad B=\ln a.
$$

Mô hình trở thành
$$
Y=Ax+B.
$$

Do đó, ta chỉ cần áp dụng thuật toán xấp xỉ đường thẳng cho dữ liệu
$$
(x_i,Y_i).
$$

Sau đó,
$$
b=A,\qquad a=e^B.
$$

### Mã giả

```text
THUẬT TOÁN BÌNH_PHƯƠNG_BÉ_NHẤT_HÀM_MŨ(x, y)

Nếu tồn tại y_i <= 0:
    Trả về "Không thể lấy logarit"

Với i = 1,...,n:
    Y_i = ln(y_i)

(A, B) = BÌNH_PHƯƠNG_BÉ_NHẤT_ĐƯỜNG_THẲNG(x, Y)

b = A
a = exp(B)

Trả về a, b
```

**Lưu ý:** Đây là bình phương bé nhất trên dữ liệu đã biến đổi logarit, không phải cực tiểu hóa trực tiếp tổng bình phương sai số trên thang $y$ ban đầu.

---

# 6. Xấp xỉ hàm lũy thừa $y=ax^b$

Giả sử
$$
a>0,\qquad x_i>0,\qquad y_i>0.
$$

Đặt
$$
X_i=\ln x_i,\qquad
Y_i=\ln y_i,
$$
và
$$
A=b,\qquad B=\ln a.
$$

Khi đó,
$$
Y=AX+B.
$$

Sau khi xấp xỉ đường thẳng trên dữ liệu logarit,
$$
b=A,\qquad a=e^B.
$$

### Mã giả

```text
THUẬT TOÁN BÌNH_PHƯƠNG_BÉ_NHẤT_HÀM_LŨY_THỪA(x, y)

Nếu tồn tại x_i <= 0 hoặc y_i <= 0:
    Trả về "Không thể lấy logarit"

Với i = 1,...,n:
    X_i = ln(x_i)
    Y_i = ln(y_i)

(A, B) = BÌNH_PHƯƠNG_BÉ_NHẤT_ĐƯỜNG_THẲNG(X, Y)

b = A
a = exp(B)

Trả về a, b
```

---

# 7. Cấu trúc hàm Python đề nghị

Notebook của chương sử dụng các hàm:

```text
least_squares_line(...)
r_squared(...)
least_squares_two_basis(...)
least_squares_polynomial(...)
least_squares_exponential(...)
least_squares_power(...)
```

Mỗi hàm trả về các tham số của mô hình và các đại lượng trung gian cần thiết để đối chiếu với phép tính trong giáo trình.
