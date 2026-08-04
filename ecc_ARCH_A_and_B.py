import random
from fastecdsa.curve import Curve

def inverse(a, p):
    """Modular inversion using Fermat's Little Theorem."""
    return pow(a, p - 2, p)

def modmul(a, b, p, p_not, r_mask, r_shift):
    """
    Montgomery modular multiplication.
    Optimized with bitwise operations since R is a power of 2.
    """
    t = a * b
    m = ((t & r_mask) * p_not) & r_mask
    u = (t + m * p) >> r_shift
    
    return u - p if u >= p else u

def point_double_jacobian_montgomery(x1, y1, z1, a_mont, q):
    """Doubles a point in Jacobian coordinates within the Montgomery domain."""
    p, p_not, r_mask, r_shift = q

    # Pre-calculate squares
    x1_sq = modmul(x1, x1, *q)
    y1_sq = modmul(y1, y1, *q)
    z1_sq = modmul(z1, z1, *q)

    z1_qu = modmul(z1_sq, z1_sq, *q)
    y1_qu = modmul(y1_sq, y1_sq, *q)

    # M = 3*x1^2 + a*z1^4
    three_x1_sq = (x1_sq * 3) % p
    a_z1_qu = modmul(a_mont, z1_qu, *q)
    m = (three_x1_sq + a_z1_qu) % p

    # S = 4*x1*y1^2
    x1_y1_sq = modmul(x1, y1_sq, *q)
    s = (x1_y1_sq * 4) % p

    # T = 8*y1^4
    t = (y1_qu * 8) % p

    # x3 = M^2 - 2S
    m_sq = modmul(m, m, *q)
    x3 = (m_sq - 2 * s) % p

    # y3 = M(S - x3) - 8T
    s_minus_x3 = (s - x3) % p
    m_s_x3 = modmul(m, s_minus_x3, *q)
    y3 = (m_s_x3 - t) % p

    # z3 = 2*y1*z1
    y1_z1 = modmul(y1, z1, *q)
    z3 = (y1_z1 * 2) % p

    return x3, y3, z3

def point_add_jacobian_montgomery(x1, y1, z1, x2, y2, z2, a_mont, q):
    """Adds two points in Jacobian coordinates within the Montgomery domain."""
    # Z == 0 indicates the point at infinity
    if z1 == 0: return x2, y2, z2
    if z2 == 0: return x1, y1, z1

    p, p_not, r_mask, r_shift = q

    z1_sq = modmul(z1, z1, *q)
    z2_sq = modmul(z2, z2, *q)
    z1_cu = modmul(z1_sq, z1, *q)
    z2_cu = modmul(z2_sq, z2, *q)

    u1 = modmul(x1, z2_sq, *q)
    u2 = modmul(x2, z1_sq, *q)
    s1 = modmul(y1, z2_cu, *q)
    s2 = modmul(y2, z1_cu, *q)

    if u1 == u2:
        if s1 == s2:
            return point_double_jacobian_montgomery(x1, y1, z1, a_mont, q)
        else:
            return 0, 0, 0  # Point at infinity

    h = (u1 - u2) % p
    r = (s1 - s2) % p

    h_sq = modmul(h, h, *q)
    h_cu = modmul(h_sq, h, *q)
    v = modmul(u1, h_sq, *q)

    # x3 = r^2 + h^3 - 2v
    r_sq = modmul(r, r, *q)
    x3 = (r_sq + h_cu - 2 * v) % p

    # y3 = r*(v - x3) - s1*h^3
    v_minus_x3 = (v - x3) % p
    r_v_x3 = modmul(r, v_minus_x3, *q)
    s1_h_cu = modmul(s1, h_cu, *q)
    y3 = (r_v_x3 - s1_h_cu) % p

    # z3 = z1 * z2 * h
    z1_z2 = modmul(z1, z2, *q)
    z3 = modmul(z1_z2, h, *q)

    return x3, y3, z3

def point_mult_jacobian(p1_x, p1_y, p, curve_a, scalar):
    """Multiplies a point by a scalar using Jacobian coordinates in Montgomery domain."""
    if scalar == 0:
        return 0, 0

    # Montgomery domain configuration
    r_shift = p.bit_length()  # Will be 512 for GOST-512
    r = 1 << r_shift
    r_mask = r - 1

    r_inv = inverse(r, p)
    p_not = (r * r_inv - 1) // p
    
    q = (p, p_not, r_mask, r_shift)

    # Convert initial point and curve parameters to Montgomery domain
    x_mont = (p1_x * r) % p
    y_mont = (p1_y * r) % p
    z_mont = r % p            # z is initially 1, so 1 * R mod p
    a_mont = (curve_a * r) % p

    # Accumulator representing point at infinity (Z=0)
    ans_x, ans_y, ans_z = 0, 0, 0
    
    # Running base point
    buf_x, buf_y, buf_z = x_mont, y_mont, z_mont

    # Double-and-Add core loop
    while scalar>0:
        if scalar & 1:
            ans_x, ans_y, ans_z = point_add_jacobian_montgomery(ans_x, ans_y, ans_z, buf_x, buf_y, buf_z, a_mont, q)
        
        scalar >>= 1
        buf_x, buf_y, buf_z = point_double_jacobian_montgomery(buf_x, buf_y, buf_z, a_mont, q)

    if ans_z == 0:
        return 0, 0

    # Exit Montgomery Domain for the final answer
    ans_x = modmul(ans_x, 1, *q)
    ans_y = modmul(ans_y, 1, *q)
    ans_z = modmul(ans_z, 1, *q)

    # Convert Jacobian back to Affine coordinates exactly once
    z_inv = inverse(ans_z, p)
    z_inv_sq = (z_inv * z_inv) % p
    z_inv_cu = (z_inv_sq * z_inv) % p
    
    final_x = (ans_x * z_inv_sq) % p
    final_y = (ans_y * z_inv_cu) % p

    return final_x, final_y

def run_testbench():
    p = 0x00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDC7
    curve_a = 0x00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDC4
    b = 0x00E8C2505DEDFC86DDC1BD0B2B6667F1DA34B82574761CB0E879BD081CFD0B6265EE3CB090F30D27614CB4574010DA90DD862EF9D4EBEE4761503190785A71C760
    q = 0x00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF27E69532F48D89116FF22B8D4E0560609B4B38ABFAD2B85DCACDB1411F10B275
    gx = 0x03
    gy = 0x7503CFE87A836AE3A61B8816E25450E6CE5E1C93ACF1ABC1778064FDCBEFA921DF1626BE4FD036E93D75E6A50E3A41E98028FE5FC235F5B889A589CB5215F2A4
    
    gost_512_paramA = Curve(
        'id-tc26-gost-3410-12-512-paramSetA', 
        p, curve_a, b, q, gx, gy
    )
    
    G = gost_512_paramA.G

    random.seed(42)  # fixed seed for reproducible test-vector generation
    for i in range(0, 100):
        scalar = random.randint(0,q-1)
        x, y = point_mult_jacobian(G.x, G.y, p, curve_a, scalar)
        P1 = scalar * G
        
        assert x == P1.x and y == P1.y, f"error P.x, P.y = {P1.x, P1.y}, func x,y = {x,y}\ninitial data {scalar}"
        print(i)
        print("golden model", hex(P1.x), hex(P1.y))
        print("got", hex(x), hex(y))

    print(f"pass")

if __name__ == "__main__":
    run_testbench()