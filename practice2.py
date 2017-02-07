k = 15.0
m = 24.0
n = 30.0
N = k+m+n

P1 = ((k/N * (k-1)/(N-1)) + (k/N * m/(N-1)) + (k/N * n/(N-1)))
print P1
P2 = ((m/N * k/(N-1)) + (m/N * (m-1)/(N-1) * (3.0/4.0)) + m/N * n/(N-1) * (2.0/4.0))
print P2
P3 = ((n/N * k/(N-1)) + (n/N * m/(N-1) * (2.0/4.0)))
print P3

P = P1 + P2 + P3
print P
