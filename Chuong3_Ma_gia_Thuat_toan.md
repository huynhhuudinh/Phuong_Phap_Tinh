# CHƯƠNG 3 — HỆ PHƯƠNG TRÌNH TUYẾN TÍNH
## Mã giả các thuật toán

> Lưu ý: để render tốt trên GitHub/MathJax, các công thức dài nên viết ở dạng block `$$ ... $$` hoặc `$$ ... $$`, và các hệ phương trình nên dùng `aligned`/`cases` rõ ràng.

Tài liệu này được xây dựng theo nội dung Chương 3 của giáo trình. Các thuật toán chính gồm:

1. Kiểm tra ma trận đường chéo trội nghiêm ngặt theo dòng.
2. Phương pháp khử Gauss.
3. Thế tiến và thế ngược.
4. Phân tích $LU$.
5. Tính các chuẩn vector và chuẩn ma trận thường dùng.
6. Phương pháp lặp đơn (Jacobi).
7. Phương pháp Gauss--Seidel.
8. Kiểm tra hội tụ bằng bán kính phổ.
9. Phương pháp Gradient Descent với bước lặp tối ưu.

---

# 1. Kiểm tra tính đường chéo trội nghiêm ngặt theo dòng

Ma trận $\mathbf A=(a_{ij})$ được gọi là đường chéo trội nghiêm ngặt theo dòng nếu

$$
|a_{ii}| >
\sum_{\substack{j=1\\j\ne i}}^n |a_{ij}|,
\quad i=1,2,\ldots,n.
$$

```text
THUẬT TOÁN KIỂM_TRA_CHÉO_TRỘI(A)

Đầu vào:
    Ma trận vuông A cấp n

Với i = 1,2,...,n:
    S = 0

    Với j = 1,2,...,n:
        Nếu j khác i:
            S = S + |a_ij|

    Nếu |a_ii| <= S:
        Trả về FALSE

Trả về TRUE
```

---

# 2. Phương pháp khử Gauss

## Ý tưởng

Đưa ma trận mở rộng

$$
(\mathbf A|\mathbf b)
$$

về dạng bậc thang bằng các phép biến đổi sơ cấp trên dòng, sau đó giải hệ tam giác trên bằng thế ngược.

## Mã giả

```text
THUẬT TOÁN KHỬ_GAUSS(A, b)

Đầu vào:
    A: ma trận vuông cấp n
    b: vector hệ số tự do

Bước 1:
    Lập ma trận mở rộng M = (A|b)

Bước 2: Khử tiến
    Với k = 1,2,...,n-1:

        Nếu phần tử trụ M[k,k] = 0:
            Tìm dòng p > k sao cho M[p,k] khác 0
            Đổi dòng k và dòng p

        Với i = k+1,...,n:
            m = M[i,k] / M[k,k]

            Với j = k,...,n+1:
                M[i,j] = M[i,j] - m*M[k,j]

Bước 3:
    Giải hệ tam giác trên bằng thế ngược

Bước 4:
    Trả về nghiệm x
```

Trong lập trình thực tế nên sử dụng **pivot từng phần**: tại cột đang xét, đổi dòng hiện tại với dòng có phần tử trị tuyệt đối lớn nhất trong cột đó.

---

# 3. Thế ngược

Giải hệ tam giác trên

$$
\mathbf U\mathbf x=\mathbf y.
$$

```text
THUẬT TOÁN THẾ_NGƯỢC(U, y)

n = số hàng của U

Với i = n,n-1,...,1:
    S = 0

    Với j = i+1,...,n:
        S = S + u_ij*x_j

    x_i = (y_i - S)/u_ii

Trả về x
```

---

# 4. Thế tiến

Giải hệ tam giác dưới

$$
\mathbf L\mathbf y=\mathbf b.
$$

```text
THUẬT TOÁN THẾ_TIẾN(L, b)

n = số hàng của L

Với i = 1,2,...,n:
    S = 0

    Với j = 1,2,...,i-1:
        S = S + l_ij*y_j

    y_i = (b_i - S)/l_ii

Trả về y
```

---

# 5. Phân tích LU

Chương sử dụng biểu diễn

$$
\mathbf A=\mathbf L\mathbf U.
$$

Trong ví dụ, $\mathbf U$ được chọn có các phần tử đường chéo chính bằng $1$, tương ứng với dạng Crout.

## Mã giả phân tích Crout

```text
THUẬT TOÁN CROUT_LU(A)

Đầu vào:
    Ma trận vuông A cấp n

Khởi tạo:
    L = ma trận 0 cấp n
    U = ma trận đơn vị cấp n

Với j = 1,2,...,n:

    # Tính cột j của L
    Với i = j,j+1,...,n:
        S = tổng từ k=1 đến j-1 của L[i,k]*U[k,j]
        L[i,j] = A[i,j] - S

    Nếu L[j,j] = 0:
        Dừng: phân tích theo dạng này không thực hiện được

    # Tính hàng j của U
    Với i = j+1,j+2,...,n:
        S = tổng từ k=1 đến j-1 của L[j,k]*U[k,i]
        U[j,i] = (A[j,i] - S)/L[j,j]

Trả về L, U
```

Sau khi có $\mathbf A=\mathbf L\mathbf U$:

```text
y = THẾ_TIẾN(L, b)
x = THẾ_NGƯỢC(U, y)
```

---

# 6. Chuẩn vector

Với

$$
\mathbf x=(x_1;x_2;\ldots;x_n),
$$

ta có

$$
\|\mathbf x\|_1=\sum_{i=1}^n|x_i|,
$$

$$
\|\mathbf x\|_\infty=\max_i|x_i|,
$$

và

$$
\|\mathbf x\|_2=
\sqrt{\sum_{i=1}^n x_i^2}.
$$

```text
THUẬT TOÁN CHUẨN_VECTOR(x)

norm1 = tổng |x_i|
norm_inf = max |x_i|
norm2 = căn bậc hai của tổng x_i^2

Trả về norm1, norm_inf, norm2
```

---

# 7. Chuẩn ma trận liên kết

Với $\mathbf A=(a_{ij})$,

$$
\|\mathbf A\|_1
=
\max_j\sum_i |a_{ij}|,
$$

$$
\|\mathbf A\|_\infty
=
\max_i\sum_j |a_{ij}|,
$$

và

$$
\|\mathbf A\|_2
=
\sqrt{\lambda_{\max}(\mathbf A^T\mathbf A)}.
$$

```text
THUẬT TOÁN CHUẨN_MA_TRẬN(A)

norm1:
    tính tổng trị tuyệt đối theo từng cột
    lấy giá trị lớn nhất

norm_inf:
    tính tổng trị tuyệt đối theo từng dòng
    lấy giá trị lớn nhất

norm2:
    tính C = A^T A
    tìm trị riêng lớn nhất lambda_max của C
    norm2 = sqrt(lambda_max)

Trả về norm1, norm_inf, norm2
```

---

# 8. Phương pháp lặp đơn (Jacobi)

Từ

$$
\mathbf A\mathbf x=\mathbf b,
$$

với $a_{ii}\ne0$, xây dựng

$$
\mathbf x=\mathbf H\mathbf x+\mathbf g,
$$

trong đó

$$
h_{ii}=0,\qquad
h_{ij}=-\frac{a_{ij}}{a_{ii}},\quad i\ne j,
$$

và

$$
g_i=\frac{b_i}{a_{ii}}.
$$

Nếu

$$
q=\|\mathbf H\|_\infty<1,
$$

thì dãy lặp hội tụ và

$$
\Delta_{\mathbf x_k}
\le
\frac{q}{1-q}
\|\mathbf x_k-\mathbf x_{k-1}\|_\infty.
$$

## Mã giả

```text
THUẬT TOÁN LẶP_ĐƠN_JACOBI(A, b, epsilon)

Bước 1:
    Kiểm tra a_ii khác 0

Bước 2:
    Xây dựng H và g:
        h_ii = 0
        h_ij = -a_ij/a_ii, i khác j
        g_i = b_i/a_ii

Bước 3:
    q = ||H||_inf

    Nếu q >= 1:
        Điều kiện đủ của chương không được thỏa mãn

Bước 4:
    Chọn x_old = g

Bước 5:
    Lặp:
        x_new = H*x_old + g

        Err = q/(1-q) * ||x_new-x_old||_inf

        Nếu Err < epsilon:
            Trả về x_new, Err

        x_old = x_new
```

---

# 9. Phương pháp Gauss--Seidel

Từ dạng

$$
\mathbf x=\mathbf H\mathbf x+\mathbf g,
$$

thành phần thứ $i$ tại bước $k$ được cập nhật bởi

$$
x_i^{(k)}
=
\sum_{j=1}^{i-1}h_{ij}x_j^{(k)}
+
\sum_{j=i}^{n}h_{ij}x_j^{(k-1)}
+
g_i.
$$

Đặt

$$
\alpha_i=\sum_{j=1}^{i-1}|h_{ij}|,
\qquad
\beta_i=\sum_{j=i}^{n}|h_{ij}|,
$$

và

$$
\gamma
=
\max_i
\frac{\beta_i}{1-\alpha_i}.
$$

Sai số được đánh giá bởi

$$
\Delta_{\mathbf x_k}
\le
\frac{\gamma}{1-\gamma}
\|\mathbf x_k-\mathbf x_{k-1}\|_\infty.
$$

## Mã giả

```text
THUẬT TOÁN GAUSS_SEIDEL(A, b, x0, epsilon)

Bước 1:
    Xây dựng H và g từ A, b

Bước 2:
    Tính gamma

Bước 3:
    x_old = x0

Bước 4:
    Lặp:
        x_new = x_old

        Với i = 1,2,...,n:
            S1 = tổng h_ij*x_new[j] với j < i
            S2 = tổng h_ij*x_old[j] với j >= i
            x_new[i] = S1 + S2 + g_i

        Err = gamma/(1-gamma) * ||x_new-x_old||_inf

        Nếu Err < epsilon:
            Trả về x_new, Err

        x_old = x_new
```

---

# 10. Ma trận lặp và tiêu chuẩn bán kính phổ

Viết

$$
\mathbf A=\mathbf D+\mathbf L+\mathbf U.
$$

Ma trận lặp Jacobi:

$$
\mathbf B_J
=
-\mathbf D^{-1}(\mathbf L+\mathbf U).
$$

Ma trận lặp Gauss--Seidel:

$$
\mathbf B_{GS}
=
-(\mathbf D+\mathbf L)^{-1}\mathbf U.
$$

Phương pháp lặp tuyến tính

$$
\mathbf x_k=\mathbf B\mathbf x_{k-1}+\mathbf c
$$

hội tụ với mọi $\mathbf x_0$ khi và chỉ khi

$$
\rho(\mathbf B)<1.
$$

## Mã giả

```text
THUẬT TOÁN BÁN_KÍNH_PHỔ(B)

Tính các trị riêng lambda_1,...,lambda_n của B

rho = max |lambda_i|

Nếu rho < 1:
    Kết luận phương pháp lặp hội tụ với mọi vector ban đầu
Ngược lại:
    Không thỏa tiêu chuẩn hội tụ

Trả về rho
```

---

# 11. Gradient Descent với bước lặp tối ưu

Áp dụng khi $\mathbf A$ đối xứng xác định dương.

Hàm mục tiêu:

$$
F(\mathbf x)
=
\frac12\mathbf x^T\mathbf A\mathbf x
-
\mathbf x^T\mathbf b.
$$

Vector phần dư:

$$
\mathbf r_k
=
\mathbf b-\mathbf A\mathbf x_k.
$$

Bước đi tối ưu:

$$
\gamma_k
=
\frac{\mathbf r_k^T\mathbf r_k}
{\mathbf r_k^T\mathbf A\mathbf r_k}.
$$

Cập nhật:

$$
\mathbf x_{k+1}
=
\mathbf x_k+\gamma_k\mathbf r_k.
$$

Điều kiện dừng trong chương:

$$
\|\mathbf r_k\|_2<\varepsilon.
$$

## Mã giả

```text
THUẬT TOÁN GRADIENT_DESCENT(A, b, epsilon)

Bước 1:
    Kiểm tra A^T = A

Bước 2:
    Kiểm tra A xác định dương

Bước 3:
    Chọn x = vector 0

Bước 4:
    Lặp:
        r = b - A*x
        Err = ||r||_2

        Nếu Err < epsilon:
            Trả về x, Err

        gamma = (r^T r)/(r^T A r)

        x = x + gamma*r
```

---

# 12. Gợi ý cấu trúc các hàm Python

Notebook của chương sử dụng các hàm:

```text
gaussian_elimination(...)
forward_substitution(...)
back_substitution(...)
crout_lu(...)
vector_norms(...)
matrix_norms(...)
jacobi_steps(...)
jacobi_until_tol(...)
gauss_seidel_steps(...)
gauss_seidel_until_tol(...)
iteration_matrices(...)
spectral_radius(...)
gradient_descent_steps(...)
gradient_descent_until_tol(...)
```

Mỗi hàm lặp nên trả về cả:
- nghiệm xấp xỉ;
- số bước lặp;
- sai số hoặc chuẩn phần dư;
- bảng lịch sử các lần lặp.
