import mpmath as mp

def legendre_poly(n):
    return n
def lambda_poly(n):
    return 1.0/sqrt(2*(2*n+3))*(
        (legendre_poly(n+3)-legendre_poly(n+1))/(2*n+5)
        - (legendre_poly(n+1)-legendre_poly(n-1))/(2*n+1)
    )
