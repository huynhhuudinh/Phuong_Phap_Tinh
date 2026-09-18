# CHƯƠNG 2 — PHƯƠNG TRÌNH PHI TUYẾN MỘT BIẾN SỐ
## Mã giả các thuật toán

> Lưu ý: để render tốt trên GitHub/MathJax, các công thức dài nên viết ở dạng block `$$ ... $$` hoặc `\[ ... \]`, và các hệ phương trình nên dùng `aligned`/`cases` rõ ràng.

Tài liệu này được xây dựng bám theo Chương 2 của giáo trình, gồm các nội dung thuật toán:
1. Kiểm tra đoạn ly nghiệm.
2. Phương pháp chia đôi.
3. Phương pháp lặp đơn.
4. Xây dựng hàm lặp theo quy tắc tổng quát.
5. Ước lượng trước số lần lặp của phương pháp lặp đơn.
6. Phương pháp tiếp tuyến (Newton).

---

## 1. Thuật toán kiểm tra đoạn ly nghiệm

### Cơ sở
Với hàm \(f\) có đạo hàm liên tục trên \([a,b]\), nếu

\[
f(a)f(b)<0
\]

và \(f'(x)\) không đổi dấu trên \([a,b]\), thì \([a,b]\) là một đoạn ly nghiệm.

### Mã giả

```text
THUẬT TOÁN KIỂM_TRA_ĐOẠN_LY_NGHIỆM(f, f', a, b)

Đầu vào:
    - Hàm f(x)
    - Đạo hàm f'(x)
    - Hai đầu mút a < b

Bước 1:
    Tính fa = f(a), fb = f(b)

Bước 2:
    Nếu fa = 0:
        Trả về: a là nghiệm đúng
    Nếu fb = 0:
        Trả về: b là nghiệm đúng

Bước 3:
    Nếu fa * fb >= 0:
        Trả về: Không thỏa điều kiện đổi dấu ở hai đầu mút

Bước 4:
    Kiểm tra dấu của f'(x) trên [a,b]

Bước 5:
    Nếu f'(x) > 0 với mọi x thuộc [a,b]
       hoặc f'(x) < 0 với mọi x thuộc [a,b]:
        Trả về: [a,b] là đoạn ly nghiệm
    Ngược lại:
        Trả về: Chưa đủ điều kiện để kết luận bằng tiêu chuẩn này
```

---

## 2. Phương pháp chia đôi

### Cơ sở
Cho \([a,b]\) là đoạn ly nghiệm và \(f(a)f(b)<0\).

Ở bước thứ \(n\),

\[
x_n=\frac{a_n+b_n}{2}.
\]

Sai số được đánh giá bởi

\[
\Delta_{x_n}\le \frac{b-a}{2^{n+1}}.
\]

### Mã giả

```text
THUẬT TOÁN CHIA_ĐÔI(f, a, b, epsilon)

Đầu vào:
    - Hàm f(x)
    - Đoạn ly nghiệm [a,b]
    - Sai số yêu cầu epsilon

Bước 1:
    Kiểm tra f(a)f(b) < 0
    Nếu không thỏa:
        Dừng và báo lỗi

Bước 2:
    n = 0

Bước 3:
    Lặp:
        x = (a + b)/2

        Nếu f(x) = 0:
            Trả về x

        Err = (b - a)/2

        Nếu Err < epsilon:
            Trả về x, Err, n

        Nếu f(a)f(x) < 0:
            b = x
        Ngược lại:
            a = x

        n = n + 1
```

### Trường hợp yêu cầu đúng N lần lặp

```text
THUẬT TOÁN CHIA_ĐÔI_N_BƯỚC(f, a, b, N)

Với n = 0,1,...,N:
    x_n = (a+b)/2
    Ghi lại a_n, b_n, x_n và cận sai số

    Nếu n = N:
        Dừng

    Nếu f(a)f(x_n) < 0:
        b = x_n
    Ngược lại:
        a = x_n
```

---

## 3. Phương pháp lặp đơn

### Cơ sở
Biến đổi phương trình \(f(x)=0\) về dạng

\[
x=\varphi(x).
\]

Nếu \(\varphi([a,b])\subset[a,b]\) và

\[
|\varphi'(x)|\le L<1,
\]

thì dãy

\[
x_n=\varphi(x_{n-1})
\]

hội tụ đến nghiệm duy nhất.

Sai số hậu nghiệm:

\[
\Delta_{x_n}
\le
\frac{L}{1-L}|x_n-x_{n-1}|.
\]

### Mã giả

```text
THUẬT TOÁN LẶP_ĐƠN(phi, x0, L, epsilon)

Đầu vào:
    - Hàm lặp phi(x)
    - Giá trị ban đầu x0
    - Hệ số co L, 0 < L < 1
    - Sai số yêu cầu epsilon

Bước 1:
    x_old = x0
    n = 0

Bước 2:
    Lặp:
        n = n + 1
        x_new = phi(x_old)

        Err = L/(1-L) * |x_new - x_old|

        Ghi lại n, x_new, Err

        Nếu Err < epsilon:
            Trả về x_new, Err, n

        x_old = x_new
```

---

## 4. Ước lượng trước số lần lặp của phương pháp lặp đơn

Sai số tiên nghiệm:

\[
\Delta_{x_n}
\le
\frac{L^n}{1-L}|x_1-x_0|.
\]

Muốn \(\Delta_{x_n}\le\varepsilon\), chọn

\[
n\ge
\left\lceil
\frac{
\ln\left(\dfrac{(1-L)\varepsilon}{|x_1-x_0|}\right)
}{
\ln L
}
\right\rceil.
\]

### Mã giả

```text
THUẬT TOÁN SỐ_BƯỚC_LẶP_TIÊN_NGHIỆM(phi, x0, L, epsilon)

Bước 1:
    x1 = phi(x0)

Bước 2:
    q = ((1-L)*epsilon) / |x1-x0|

Bước 3:
    n = ceil( ln(q) / ln(L) )

Bước 4:
    Trả về n
```

---

## 5. Xây dựng hàm lặp theo quy tắc tổng quát

Giả sử \(f'(x)\) liên tục và không đổi dấu trên đoạn ly nghiệm \([a,b]\). Đặt

\[
m=\min_{a\le x\le b}|f'(x)|,
\qquad
M=\max_{a\le x\le b}|f'(x)|.
\]

### Trường hợp 1: \(f'(x)>0\)

Chọn

\[
\varphi(x)=x-\frac{1}{M}f(x).
\]

Khi đó có thể lấy

\[
L=1-\frac{m}{M}.
\]

### Trường hợp 2: \(f'(x)<0\)

Chọn

\[
\varphi(x)=x+\frac{1}{M}f(x),
\]

và cũng có thể lấy

\[
L=1-\frac{m}{M}.
\]

### Mã giả

```text
THUẬT TOÁN TẠO_HÀM_LẶP_TỔNG_QUÁT(f, f', a, b)

Bước 1:
    Tính
        m = min |f'(x)| trên [a,b]
        M = max |f'(x)| trên [a,b]

Bước 2:
    Xác định dấu của f'(x)

Bước 3:
    Nếu f'(x) > 0 trên [a,b]:
        phi(x) = x - f(x)/M

    Nếu f'(x) < 0 trên [a,b]:
        phi(x) = x + f(x)/M

Bước 4:
    L = 1 - m/M

Bước 5:
    Trả về phi, L, m, M
```

---

## 6. Phương pháp tiếp tuyến (Newton)

### Cơ sở

\[
x_n=x_{n-1}-\frac{f(x_{n-1})}{f'(x_{n-1})}.
\]

Theo điều kiện của chương, chọn \(x_0\in\{a,b\}\) sao cho

\[
f(x_0)f''(x_0)>0.
\]

Đặt

\[
m=\min_{a\le x\le b}|f'(x)|,
\qquad
M=\max_{a\le x\le b}|f''(x)|.
\]

Sai số:

\[
\Delta_{x_n}
\le
\frac{M}{2m}(x_n-x_{n-1})^2.
\]

### Mã giả

```text
THUẬT TOÁN NEWTON(f, f', f'', a, b, epsilon)

Đầu vào:
    - f, f', f''
    - Đoạn ly nghiệm [a,b]
    - Sai số yêu cầu epsilon

Bước 1:
    Kiểm tra f' và f'' liên tục, không đổi dấu trên [a,b]

Bước 2:
    Nếu f(a)f''(a) > 0:
        x_old = a
    Ngược lại:
        x_old = b

Bước 3:
    Tính
        m = min |f'(x)| trên [a,b]
        M = max |f''(x)| trên [a,b]

Bước 4:
    n = 0

Bước 5:
    Lặp:
        n = n + 1

        Nếu f'(x_old) = 0:
            Dừng và báo lỗi

        x_new = x_old - f(x_old)/f'(x_old)

        Err = M/(2m) * (x_new - x_old)^2

        Ghi lại n, x_new, Err

        Nếu Err < epsilon:
            Trả về x_new, Err, n

        x_old = x_new
```

---

## 7. Gợi ý cấu trúc triển khai chương trình

Để thuận tiện khi học và tái sử dụng code, nên viết các hàm:

```text
check_isolating_interval(...)
bisection_steps(...)
bisection_until_tol(...)
fixed_point_steps(...)
fixed_point_until_tol(...)
a_priori_iterations(...)
newton_steps(...)
newton_until_tol(...)
```

Mỗi hàm nên trả về cả:
- nghiệm gần đúng;
- số lần lặp;
- cận sai số;
- bảng lịch sử các lần lặp.

Cách tổ chức này được sử dụng trong notebook `Chuong2_Code_Tat_ca_Vi_du.ipynb`.
