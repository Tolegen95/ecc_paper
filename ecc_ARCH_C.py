from tinyec import ec
import random
from fastecdsa.curve import Curve
from fastecdsa.point import Point

def inverse(a,p):
    return pow(a,p-2,p)    

def montgomery_ladder_step(XQP, XRP, M, YP, p):
    """
    Performs a single step of the Montgomery ladder based on the provided 
    ladder state and field prime p.
    """
    
    # 1. YR_bar = YP + 2 * M * XRP
    YR_bar = (YP + 2 * M * XRP) % p
    
    # 2. E = XQP - XRP
    E = (XQP - XRP) % p
    
    # 3. F = YR_bar * E
    F = (YR_bar * E) % p
    
    # 4. G = E^2
    G = pow(E, 2, p)
    
    # 5. XRP_prime = XRP * G
    XRP_prime = (XRP * G) % p
    
    # 6. H = YR_bar^2
    H = pow(YR_bar, 2, p)
    
    # 7. M_prime = M * F
    M_prime = (M * F) % p
    
    # 8. YP_prime = YP * F * G
    YP_prime = (YP * F * G) % p
    
    # 9. K = H + M_prime
    K = (H + M_prime) % p
    
    # 10. L = K + M_prime
    L = (K + M_prime) % p
    
    # 11. M_double_prime = XRP_prime - K
    M_double_prime = (XRP_prime - K) % p
    
    # 12. XSP = H * L
    XSP = (H * L) % p
    
    # 13. XTP = XRP_prime^2 + YP_prime
    XTP = (pow(XRP_prime, 2, p) + YP_prime) % p
    
    # 14. YP_double_prime = YP_prime * H
    YP_double_prime = (YP_prime * H) % p
    
    return XSP, XTP, M_double_prime, YP_double_prime
    
def run_testbench(k, priv_key):    
    p = 0x00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD97
    a = 0x00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD94
    b = 0xA6
    m = 0x00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF6C611070995AD10045841B09B761B893
    q = 0x00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF6C611070995AD10045841B09B761B893
    x = 0x01
    y = 0x008D91E471E0989CDA27DF505A453F2B7635294F2DDF23E3B122ACC99C9E9F1E14
    
    gost_256_paramB = Curve(
        'id-tc26-gost-3410-2012-256-paramSetB', 
        p, a, b, q, x, y
    )
    
    G = gost_256_paramB.G
    Q=priv_key*G
    answ = k*Q
        
    print("Nonce for signing",hex(k))
    print("Priv key",hex(priv_key))
    
    k = (k - (1<<256))%q
    
    Z_sq = ((2*Q.y)**2)%p
    mZ = (3*Q.x**2+a)%p
    Xrp = ((mZ**2)-3*Q.x*Z_sq)%p
    Y = (Z_sq**2)%p
    
    XQP = 0
    XRP = (Xrp)%p
    M = (mZ)%p
    YP = (Y)%p
    
    for i in range(256):
        bit = (k >> (255 - i)) & 1
        if bit == 1:
            XQP, XRP, M, YP = montgomery_ladder_step(XQP, XRP, M, YP, p)
        else:
            XRP, XQP, M, YP = montgomery_ladder_step(XRP, XQP, M, YP, p)
            
    numer = (2*Q.y*(M**2 - XQP - XRP))%p
    denom = (3*Q.x * (YP))%p
            
    Z_inv = ((numer)*inverse(denom,p))%p

    x_q = (Q.x + XQP*Z_inv*Z_inv)%p

    assert x_q == answ.x, "Error"
    
    print("got x", hex(x_q))
    

if __name__ == "__main__":
    q = 0x00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF6C611070995AD10045841B09B761B893
    random.seed(42)  # fixed seed for reproducible test-vector generation
    for i in range(100):
        k = random.randint(0,q-1)
        priv_key = random.randint(0,q-1)
        run_testbench(k, priv_key)