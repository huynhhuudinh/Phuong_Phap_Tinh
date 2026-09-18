# CHƯƠNG 6 — TÍNH GẦN ĐÚNG ĐẠO HÀM VÀ TÍCH PHÂN
## Mã giả các thuật toán

> Lưu ý: để render tốt trên GitHub/MathJax, các công thức dài nên viết ở dạng block `$$ ... $$` hoặc `$$ ... $$`, và các hệ phương trình nên dùng `aligned`/`cases` rõ ràng.

Tài liệu này được xây dựng theo nội dung Chương 6 của giáo trình.

Các thuật toán chính gồm:

1. Tỉ sai phân tiến và lùi.
2. Công thức đạo hàm 3-điểm cuối.
3. Công thức đạo hàm 3-điểm giữa.
4. Công thức đạo hàm 5-điểm cuối.
5. Công thức đạo hàm 5-điểm giữa.
6. Công thức Newton--Cotes đóng.
7. Công thức hình thang.
8. Công thức hình thang mở rộng.
9. Công thức Simpson một phần ba.
10. Công thức Simpson một phần ba mở rộng.
11. Phương pháp tích phân Romberg.
12. Công thức cầu phương Gauss--Legendre trên \([-1,1]\).
13. Công thức cầu phương Gauss--Legendre trên đoạn \([a,b]\).

---

# 1. Tỉ sai phân tiến

Công thức:

$$
f'(x_0)
\approx
\frac{f(x_0+h)-f(x_0)}{h},
\qquad h>0.
$$

Nếu \(f''\) liên tục trên đoạn nối \(x_0\) và \(x_0+h\), thì

$$
\Delta_D
\le
\frac{M|h|}{2},
$$

trong đó

$$
M=
\max |f''(x)|.
$$

## Mã giả

```text
THUẬT TOÁN SAI_PHÂN_TIẾN(f, x0, h)

Nếu h <= 0:
    báo lỗi

D = [f(x0+h) - f(x0)] / h

Trả về D
```

---

# 2. Tỉ sai phân lùi

Công thức có thể viết dưới dạng

$$
f'(x_0)
\approx
\frac{f(x_0)-f(x_0-h)}{h},
\qquad h>0.
$$

## Mã giả

```text
THUẬT TOÁN SAI_PHÂN_LÙI(f, x0, h)

Nếu h <= 0:
    báo lỗi

D = [f(x0) - f(x0-h)] / h

Trả về D
```

---

# 3. Công thức 3-điểm cuối

Theo ký hiệu của chương:

$$
f'(x^*)
\approx
\frac{-3f(x^*)+4f(x^*+h)-f(x^*+2h)}{2h}.
$$

Giá trị \(h\) có thể dương hoặc âm nếu các mốc tương ứng tồn tại.

## Mã giả

```text
THUẬT TOÁN ĐẠO_HÀM_3_ĐIỂM_CUỐI(f, x_star, h)

D =
[-3*f(x_star)
 +4*f(x_star+h)
 -f(x_star+2h)]
/
(2*h)

Trả về D
```

---

# 4. Công thức 3-điểm giữa

$$
f'(x^*)
\approx
\frac{-f(x^*-h)+f(x^*+h)}{2h}.
$$

## Mã giả

```text
THUẬT TOÁN ĐẠO_HÀM_3_ĐIỂM_GIỮA(f, x_star, h)

D =
[-f(x_star-h)+f(x_star+h)]
/
(2*h)

Trả về D
```

---

# 5. Công thức 5-điểm cuối

$$
f'(x^*)
\approx
\frac{
-25f(x^*)
+48f(x^*+h)
-36f(x^*+2h)
+16f(x^*+3h)
-3f(x^*+4h)
}{
12h
}.
$$

## Mã giả

```text
THUẬT TOÁN ĐẠO_HÀM_5_ĐIỂM_CUỐI(f, x_star, h)

D =
[-25*f(x_star)
 +48*f(x_star+h)
 -36*f(x_star+2h)
 +16*f(x_star+3h)
 -3*f(x_star+4h)]
/
(12*h)

Trả về D
```

---

# 6. Công thức 5-điểm giữa

$$
f'(x^*)
\approx
\frac{
f(x^*-2h)-8f(x^*-h)+8f(x^*+h)-f(x^*+2h)
}{
12h
}.
$$

## Mã giả

```text
THUẬT TOÁN ĐẠO_HÀM_5_ĐIỂM_GIỮA(f, x_star, h)

D =
[f(x_star-2h)
 -8*f(x_star-h)
 +8*f(x_star+h)
 -f(x_star+2h)]
/
(12*h)

Trả về D
```

---

# 7. Công thức Newton--Cotes đóng

Chia đoạn \([a,b]\) thành \(n\) đoạn bằng nhau:

$$
x_i=a+ih,
\qquad
h=\frac{b-a}{n}.
$$

Công thức:

$$
I
\approx
I_{\rm NC}
=
h\sum_{i=0}^n H_{n,i}f(x_i).
$$

Các hệ số \(H_{n,i}\) là hệ số Cotes.

## Mã giả

```text
THUẬT TOÁN NEWTON_COTES_ĐÓNG(f, a, b, n, H)

h = (b-a)/n
S = 0

Với i = 0,...,n:
    x_i = a + i*h
    S = S + H[i]*f(x_i)

I = h*S

Trả về I
```

### Cách tính hệ số Cotes bằng hệ moment

Ta có thể xác định \(H_i\) từ các điều kiện

$$
\sum_{i=0}^n H_i i^k
=
\frac{n^{k+1}}{k+1},
\qquad
k=0,1,\ldots,n.
$$

```text
THUẬT TOÁN HỆ_SỐ_COTES(n)

Tạo ma trận A kích thước (n+1)x(n+1)

Với k = 0,...,n:
    Với i = 0,...,n:
        A[k,i] = i^k

    b[k] = n^(k+1)/(k+1)

Giải hệ:
    A*H = b

Trả về H
```

---

# 8. Công thức hình thang

$$
I
=
\int_a^b f(x)\,dx
\approx
\frac{b-a}{2}[f(a)+f(b)].
$$

Sai số:

$$
\Delta_{\rm HT}
\le
\frac{M(b-a)^3}{12},
$$

với

$$
M=\max_{a\le x\le b}|f''(x)|.
$$

## Mã giả

```text
THUẬT TOÁN HÌNH_THANG(f, a, b)

I = (b-a)/2 * [f(a)+f(b)]

Trả về I
```

---

# 9. Công thức hình thang mở rộng

Chia \([a,b]\) thành \(n\) đoạn bằng nhau:

$$
h=\frac{b-a}{n}.
$$

Công thức:

$$
I_{\rm HTMR}
=
\frac{h}{2}
\left[
f(x_0)
+
2\sum_{i=1}^{n-1}f(x_i)
+
f(x_n)
\right].
$$

Sai số:

$$
\Delta_{\rm HTMR}
\le
\frac{(b-a)Mh^2}{12}.
$$

## Mã giả

```text
THUẬT TOÁN HÌNH_THANG_MỞ_RỘNG(f, a, b, n)

h = (b-a)/n

S = f(a) + f(b)

Với i = 1,...,n-1:
    x_i = a + i*h
    S = S + 2*f(x_i)

I = h*S/2

Trả về I
```

---

# 10. Công thức Simpson một phần ba

$$
I
\approx
\frac{b-a}{6}
\left[
f(a)
+
4f\left(\frac{a+b}{2}\right)
+
f(b)
\right].
$$

Sai số thường dùng:

$$
\Delta_{\rm S1/3}
\le
\frac{M(b-a)^5}{2880},
$$

với

$$
M=\max |f^{(4)}(x)|.
$$

## Mã giả

```text
THUẬT TOÁN SIMPSON_1_3(f, a, b)

m = (a+b)/2

I =
(b-a)/6 *
[f(a) + 4*f(m) + f(b)]

Trả về I
```

---

# 11. Công thức Simpson một phần ba mở rộng

Yêu cầu số đoạn chia \(n\) phải chẵn.

$$
h=\frac{b-a}{n}.
$$

$$
I_{\rm S1/3MR}
=
\frac{h}{3}
\left[
f(x_0)+f(x_n)
+
4\sum_{\substack{i=1\\i\ {\rm lẻ}}}^{n-1}f(x_i)
+
2\sum_{\substack{i=2\\i\ {\rm chẵn}}}^{n-2}f(x_i)
\right].
$$

Sai số:

$$
\Delta_{\rm S1/3MR}
\le
\frac{(b-a)Mh^4}{180}.
$$

## Mã giả

```text
THUẬT TOÁN SIMPSON_1_3_MỞ_RỘNG(f, a, b, n)

Nếu n không chẵn:
    báo lỗi

h = (b-a)/n
S = f(a) + f(b)

Với i = 1,...,n-1:
    x_i = a + i*h

    Nếu i lẻ:
        S = S + 4*f(x_i)
    Ngược lại:
        S = S + 2*f(x_i)

I = h*S/3

Trả về I
```

---

# 12. Phương pháp tích phân Romberg

Khởi đầu bằng công thức hình thang mở rộng.

Đặt

$$
R_{i,0}
=
I_{\rm HTMR}(2^i).
$$

Sau đó thực hiện ngoại suy Richardson:

$$
R_{i,k}
=
\frac{
4^kR_{i,k-1}-R_{i-1,k-1}
}{
4^k-1
},
\qquad
k=1,\ldots,i.
$$

## Mã giả

```text
THUẬT TOÁN ROMBERG(f, a, b, m)

Tạo bảng R có kích thước (m+1)x(m+1)

R[0,0] = HÌNH_THANG(f,a,b)

Với i = 1,...,m:
    n = 2^i
    h = (b-a)/n

    S = 0

    Với j = 1,...,2^(i-1):
        x_j = a + (2j-1)*h
        S = S + f(x_j)

    R[i,0] = 0.5*R[i-1,0] + h*S

    Với k = 1,...,i:
        R[i,k]
        =
        [4^k*R[i,k-1] - R[i-1,k-1]]
        /
        (4^k-1)

Trả về bảng R và R[m,m]
```

---

# 13. Cầu phương Gauss--Legendre trên \([-1,1]\)

Công thức:

$$
I
=
\int_{-1}^{1}f(x)\,dx
\approx
\sum_{i=1}^n B_if(x_i),
$$

trong đó \(x_i\) là các nghiệm của đa thức Legendre \(P_n\), còn \(B_i\) là các trọng số Gauss.

## Mã giả

```text
THUẬT TOÁN GAUSS_LEGENDRE_CHUẨN(f, nodes, weights)

I = 0

Với i = 1,...,n:
    I = I + weights[i]*f(nodes[i])

Trả về I
```

---

# 14. Cầu phương Gauss--Legendre trên \([a,b]\)

Đổi biến

$$
x=
\frac{b-a}{2}t
+
\frac{a+b}{2}.
$$

Khi đó

$$
I
=
\frac{b-a}{2}
\int_{-1}^{1}
f\left(
\frac{b-a}{2}t+\frac{a+b}{2}
\right)dt.
$$

## Mã giả

```text
THUẬT TOÁN GAUSS_LEGENDRE(f, a, b, nodes, weights)

c = (a+b)/2
d = (b-a)/2

I = 0

Với i = 1,...,n:
    x_i = c + d*nodes[i]
    I = I + weights[i]*f(x_i)

I = d*I

Trả về I
```

---

# 15. Đánh giá sai số Gauss--Legendre

Với công thức Gauss \(n\)-điểm trên \([-1,1]\),

$$
\Delta_{\rm G}
\le
\frac{
2^{2n+1}(n!)^4
}{
[(2n)!]^3(2n+1)
}
M,
$$

trong đó

$$
M=
\max_{-1\le x\le1}|f^{(2n)}(x)|.
$$

```text
THUẬT TOÁN CẬN_SAI_SỐ_GAUSS(n, M)

C =
2^(2n+1) * (n!)^4
/
[((2n)!)^3 * (2n+1)]

Err = C*M

Trả về Err
```

---

# 16. Cấu trúc các hàm Python trong notebook

```text
forward_difference(...)
backward_difference(...)
three_point_endpoint(...)
three_point_midpoint(...)
five_point_endpoint(...)
five_point_midpoint(...)

cotes_weights(...)
closed_newton_cotes(...)
trapezoid(...)
composite_trapezoid(...)
simpson_one_third(...)
composite_simpson_one_third(...)

romberg(...)
gauss_legendre(...)
gauss_error_bound(...)
```
