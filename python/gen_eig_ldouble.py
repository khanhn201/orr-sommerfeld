import mpmath as mp
import numpy as np

import numpy.linalg as la
import scipy.linalg as sla

from python.my_timer import tic, toc, print_all_timers

def longdouble_to_mp(x, dps=50):
    """
    Convert np.longdouble / np.clongdouble (scalar or 1D/2D array)
    to mpmath (mpf/mpc or mp.matrix).

    - real longdouble -> mp.mpf
    - complex longdouble -> mp.mpc
    - 1D array -> list of mpf/mpc
    - 2D array -> mp.matrix

    dps: mpmath decimal precision to use (should be >= precision of long double).
    """
    tmr_tag0 = "convert:np_to_mp"
    t0 = tic(tmr_tag0)

    # ensure mpmath has enough precision
    mp.mp.dps = max(mp.mp.dps, dps)
    precision = 36

    def np_format_float_scientific(x):
        return np.format_float_scientific(
                 np.longdouble(x),
                 precision=precision,
                 unique=False,
                 trim='k',
               )

    def _toc(type1=None,type2=None):
        tag = tmr_tag0
        if type1 is not None:
            tag += ":" + type1
        if type2 is not None:
            tag += "(" + type2 + ")"
        toc(tag, t0=t0, out=False)

    # --- scalar case ---
    if np.isscalar(x):
        # complex scalar?
        if np.iscomplexobj(x):
            z = np.clongdouble(x)
            s_re = np_format_float_scientific(z.real)
            s_im = np_format_float_scientific(z.imag)
            s_out = mp.mpc(mp.mpf(s_re), mp.mpf(s_im))
            _toc("scalar","d")
            return s_out
        else:
            s = np_format_float_scientific(x)
            s_out = mp.mpf(s)
            _toc("scalar","c")
            return s_out

    # --- array case ---
    arr = np.asarray(x)

    # decide if array is complex
    is_complex = np.iscomplexobj(arr)
    ctag = "c" if is_complex else "d"

    # handle 1D
    if arr.ndim == 1:
        out = []
        if is_complex:
            arr = arr.astype(np.clongdouble)
            for z in arr:
                s_re = np_format_float_scientific(z.real)
                s_im = np_format_float_scientific(z.imag)
                out.append(mp.mpc(mp.mpf(s_re), mp.mpf(s_im)))
        else:
            arr = arr.astype(np.longdouble)
            for r in arr:
                s = np_format_float_scientific(r)
                out.append(mp.mpf(s))
        _toc("array1d",ctag)
        return out

    # handle 2D -> mp.matrix
    if arr.ndim == 2:
        m, n = arr.shape
        if is_complex:
            arr = arr.astype(np.clongdouble)
            rows = []
            for i in range(m):
                row = []
                for j in range(n):
                    z = arr[i, j]
                    s_re = np_format_float_scientific(z.real)
                    s_im = np_format_float_scientific(z.imag)
                    row.append(mp.mpc(mp.mpf(s_re), mp.mpf(s_im)))
                rows.append(row)
            out = mp.matrix(rows)

        else:
            arr = arr.astype(np.longdouble)
            rows = []
            for i in range(m):
                row = []
                for j in range(n):
                    r = arr[i, j]
                    s = np_format_float_scientific(r)
                    row.append(mp.mpf(s))
                rows.append(row)
            out = mp.matrix(rows)
        _toc("array2d",ctag)
        return out

    raise ValueError("Only 0D, 1D, or 2D inputs are supported.")


def mp_to_longdouble(x, digits=40):
    """
    Convert mpmath objects to NumPy longdouble / clongdouble.

    Supports:
      - mp.mpf  -> np.longdouble
      - mp.mpc  -> np.clongdouble
      - mp.matrix -> 2D np.ndarray (longdouble or clongdouble)
      - 1D iterable of mpf/mpc -> 1D np.ndarray

    Anything else will raise TypeError (on purpose, to avoid
    accidentally feeding containers to np.longdouble).
    """
    tmr_tag0 = "convert:mp_to_np"
    t0 = tic(tmr_tag0)

    is_complex = False

    def _mp_scalar_to_longdouble(v):
        """mpf/mpc -> np.longdouble/np.clongdouble (via high-precision strings)."""
        if isinstance(v, mp.mpf):
            return np.longdouble(mp.nstr(v, digits))
        if isinstance(v, mp.mpc):
            is_complex = True
            re = np.longdouble(mp.nstr(v.real, digits))
            im = np.longdouble(mp.nstr(v.imag, digits))
            return np.clongdouble(re + 1j * im)
        raise TypeError(f"Expected mpf/mpc, got {type(v)}")

    def _toc(type1=None,type2=None):
        tag = tmr_tag0
        if type1 is not None:
            tag += ":" + type1
        if type2 is not None:
            tag += "(" + type2 + ")"
        toc(tag, t0=t0, out=False)

    ctag = "c" if is_complex else tmr_tag0 + "d"
    
    # --- scalar cases ---
    if isinstance(x, (mp.mpf, mp.mpc)):
        out = _mp_scalar_to_longdouble(x)
        _toc("scalar",ctag)
        return out

    # --- mpmath matrix case (2D) ---
    if isinstance(x, mp.matrix):
        m, n = x.rows, x.cols

        # detect complex entries
        has_complex = any(isinstance(x[i, j], mp.mpc)
                          for i in range(m) for j in range(n))
        ctag = "c" if has_complex else tmr_tag0 + "d"
        dtype = np.clongdouble if has_complex else np.longdouble

        out = np.empty((m, n), dtype=dtype)
        for i in range(m):
            for j in range(n):
                out[i, j] = _mp_scalar_to_longdouble(x[i, j])
        _toc("array2d",ctag)
        return out

    # --- 1D iterable of mpf/mpc ---
    if isinstance(x, (list, tuple)):
        if not x:
            return np.array([], dtype=np.longdouble)
        has_complex = any(isinstance(v, mp.mpc) for v in x)
        ctag = "c" if has_complex else tmr_tag0 + "d"
        dtype = np.clongdouble if has_complex else np.longdouble

        out = np.empty(len(x), dtype=dtype)
        for i, v in enumerate(x):
            out[i] = _mp_scalar_to_longdouble(v)
        _toc("array1d",ctag)
        return out


def test_eig(A, lamb, V, dtype=np.longdouble):
    diff = A * V - V * mp.diag(lamb)
    err = mp.norm(diff)
    return mp_to_longdouble(err)

def test_eig2(A, B, lamb, V):
    diff = A * V - B * V * mp.diag(lamb)
    err = mp.norm(diff)
    return mp_to_longdouble(err)
    
def mp_left_solve(B, A):
    """Return C = B^{-1} A without forming B^{-1} explicitly."""
    m, n = A.rows, A.cols
    C = mp.zeros(B.rows, n)
    for j in range(n):
        # solve B x = A[:, j]
        x = mp.lu_solve(B.copy(), A[:, j])
        for i in range(B.rows):
            C[i, j] = x[i]
    return C

def mp_eig(A):
    A = longdouble_to_mp(A)
    lamb, V = mp.eig(A)

    lamb = mp_to_longdouble(lamb)
    V = mp_to_longdouble(V)
    return lamb, V

def mp_eig2(A,B):
    A = longdouble_to_mp(A)
    B = longdouble_to_mp(B)
    
    tic('mp lsolve')
    M = mp_left_solve(B, A)
    tic('mp lsolve')

    print('ttt',mp.norm(B*M-A))
    lamb, V = mp.eig(M)

    err1 = test_eig(M,lamb,V)
    err2 = test_eig2(A,B,lamb,V)
    print('eig solve err',err1,err2)

    lamb = mp_to_longdouble(lamb)
#    print(V)
#    print(mp.nstr(V, 20))
    V = mp_to_longdouble(V)
#    print(V)

    return lamb, V

def eigenpy_eig(A):
    import eigenpy
    from eigenpy import EigenSolver
    A = np.double(A)
    solver = EigenSolver(A)

    w = solver.eigenvalues()   # shape (n,)
    V = solver.eigenvectors() 
    return w, V

def eigenpy_eig2(A,B):
    import eigenpy
    from eigenpy import GeneralizedEigenSolver
    A = np.double(A)
    B = np.double(B)
    solver = GeneralizedEigenSolver(A,B)
    
    w = solver.eigenvalues()   # shape (n,)
    V = solver.eigenvectors()
    return w, V


def my_solve_ext(A,b):
    n = A.shape[0]

    tic("ls_solve:convert:to_mp")
    A = longdouble_to_mp(A.copy())
    b = longdouble_to_mp(b.copy())
    toc("ls_solve:convert:to_mp")

    tic("ls_solve:solve")
    x = mp.lu_solve(A,b)
    toc("ls_solve:solve")

#    x, qr_norm = mp.qr_solve(A,b) # expansive, but more qccurate?
    tic("ls_solve:convert:from_mp")
    x = mp_to_longdouble(x)
    toc("ls_solve:convert:from_mp")
    return x.reshape((-1))

def my_solve_krylov(A,      # callable: w = A(z)
    b,                      # (n,) right-hand side, np array (real/complex)
    tol=None,               # relative tolerance (on weighted 2-norm)
    maxit=200,              # maximum iterations
    dtype=np.clongdouble,
    verbose=False):  # work in complex long double
    """
    Minimal Residual + Krylov subspace (MR) in complex long-double.
    A may be non-symmetric; inner product may be weighted by wt.
    """

    # --- helpers ---
    eps = np.finfo(np.longdouble).eps  # for the small-beta guard
    A = np.asarray(A, dtype=dtype)
    b = np.asarray(b, dtype=dtype)
    n = b.size

    A64 = np.asarray(A, dtype=np.complex128)
    LU, piv = sla.lu_factor(A64, check_finite=True)

    def fun_prec(r):
        return sla.lu_solve((LU, piv), r).astype(np.clongdouble, copy=False)

    def wdot(u_, v_):
        return np.vdot(u_, v_) # Hermitian (conj on first arg)

    # Classical Gram-Schmidt
    def orthogonal_cgs(z,w,P,Q,k):
        if iterk > 1:
            w0 = w.copy()
            betas = np.empty(k - 1, dtype=dtype)
            for j in range(k - 1):
                betas[j] = wdot(Q[:, j], w0)
            z -= P[:, :k-1] @ betas
            w -= Q[:, :k-1] @ betas

        # CGS-2 (second pass)
        if k > 1 and two_pass:
            w1 = w.copy()
            betas2 = np.empty(iterk - 1, dtype=dtype)
            for j in range(iterk - 1):
                betas2[j] = wdot(Q[:, j], w1)
            z -= P[:, :k-1] @ betas2
            w -= Q[:, :k-1] @ betas2
        return z, w

    # Modified Gram-Schmidt
    def orthogonal_mgs(z,w,P,Q,k):
        for j in range(k - 1):
            beta = wdot(Q[:, j], w)
            z -= beta * P[:, j]
            w -= beta * Q[:, j]

        if two_pass: # Pass 2
            for j in range(k - 1):
                beta = wdot(Q[:, j], w)
                z -= beta * P[:, j]
                w -= beta * Q[:, j]
        return z, w

    two_pass = True
    use_cgs = False
    orthogonal = orthogonal_cgs if use_cgs else orthogonal_mgs

    if tol is None:
        tol = eps * 5

    # --- initial residual ---
#    u = np.zeros(n, dtype=dtype)
    u = fun_prec(b)
    r = b.copy()
    if np.max(np.abs(u)) > eps:
        r = b - A @ u

    res0 = np.sqrt(wdot(r, r).real)
    res = res0
    rtol = tol * res0

    if res0 == 0:
        return u

    # storage for Krylov bases
    P = np.zeros((n, maxit), dtype=dtype)
    Q = np.zeros((n, maxit), dtype=dtype)

    iterk = 0
    if verbose:
        print(iterk, res, rtol)

    while res > rtol and iterk < maxit:

        z = fun_prec(r)
        w = A @ z

        # project
        z, w = orthogonal(z,w,P,Q,iterk)
        beta = np.sqrt(wdot(w, w).real)
#        if np.abs(beta) < 5 * eps:
#            break

        # update
        P[:, iterk - 1] = z / beta
        Q[:, iterk - 1] = w / beta

        alpha = wdot(Q[:, iterk - 1], r)
        u = u + alpha * P[:, iterk - 1]
        r = r - alpha * Q[:, iterk - 1]

        res = np.sqrt(wdot(r, r).real)

        iterk += 1
        if verbose:
            print(iterk,res,res/res0)

    return u

def my_solve_ir(A, b, maxit=20, verbose=True):
    A = np.asarray(A, dtype=np.clongdouble)
    b = np.asarray(b, dtype=np.clongdouble)

    A64 = np.asarray(A, dtype=np.complex128)
    LU, piv = sla.lu_factor(A64, check_finite=True)
    
    def wdot(u_, v_):
        return np.vdot(u_, v_)

    def solve64(r):
        return sla.lu_solve((LU, piv), r).astype(np.clongdouble, copy=False)

    x = solve64(b)
    r = b - A @ x
    res = np.sqrt(wdot(r, r).real)
    if verbose:
        print(0,res)

    for i in range(1,maxit+1):

        d = solve64(r)
        x = x + d
        r = b - A @ x
        res = np.sqrt(wdot(r, r).real)

        if verbose:
            print(i,res)
    return x

def my_solve_ir_ext(A, b, maxit=20, verbose=False):
    n = A.shape[0]
    A_mp = longdouble_to_mp(A)
    b_mp = longdouble_to_mp(b.reshape((n,1)))

    A64 = np.asarray(A, dtype=np.complex128)
    LU, piv = sla.lu_factor(A64, check_finite=True)
    def solve64(r):
        r = r.astype(np.complex128, copy=False)
        return sla.lu_solve((LU, piv), r).astype(np.clongdouble, copy=False)

    def solve(r_mp):
        r = mp_to_longdouble(r_mp)
        r = r.reshape((-1))
        x = my_solve_krylov(A, r, verbose=False)
#        x = solve64(r)
        x = x.reshape((n,1))
        return longdouble_to_mp(x)

    def norm(r_mp):
        r = mp_to_longdouble(r_mp)
        return np.sqrt(np.vdot(r, r).real)

    x_mp = solve(b_mp)
    r_mp = b_mp - A_mp * x_mp
    res = norm(r_mp)
    if verbose:
        print(0,res)
    
    for i in range(1,maxit+1):
        d_mp = solve(r_mp)
        x_mp = x_mp + d_mp
        r_mp = b_mp - A_mp * x_mp
        res = norm(r_mp)

        if verbose:
            print(i,res)
    return mp_to_longdouble(x_mp).reshape((-1))

def my_eig(A, w0=None, V0=None): # TODO: general eigs
    A = np.clongdouble(A)

    # initial eigenvalues/vectors from fp64
    # if (w0) or (V0==None): 
    w0, V0 = sla.eig(A)

    w = np.clongdouble(w0)
    V = np.clongdouble(V0)
    
    n = A.shape[0]

    I = np.eye(n, dtype = np.clongdouble)
    def improve(A,v,sigma,maxit=50,upd_shift=False):
        y = v
        for k in range(maxit): # shifted inverse method
            v = y / np.sqrt(np.vdot(y,y))
            # y = my_solve_ext(A - sigma * I, v)    # mp, accurate, slow
#            y = my_solve_ir(A - sigma * I, v)     # long double ir, not converging
#            y = my_solve_ir_ext(A - sigma * I, v) # mp, ir, not converging
            y = my_solve_krylov(A - sigma * I, v) # long double ls solve, fast, lost a bit accuracy
            theta = np.vdot(v, y)

            # update shift results ill-conditioned linear system
            if upd_shift:
                sigma = sigma + 1. / theta
                y = y / theta

        if not upd_shift:
            sigma = sigma + 1. / theta
            y = y / theta
        return y, sigma


    for j in range(n):
        y, sigma = improve(A, V[:,j], w[j])
        V[:,j] = y.reshape((-1))
        w[j] = sigma

    return w, V


def test(A, B, seed=None, symmetric=False):
    def make_random_spd(n, seed=None, dtype=np.double):
        """Random symmetric positive definite matrix B."""
        rng = np.random.default_rng(seed)
        R = rng.standard_normal((n, n),dtype=dtype)
        B = R.T @ R  # symmetric positive definite
        # Add a small multiple of I for conditioning
        B += 1e-1 * np.eye(n,dtype=dtype)
        return B
    
    def make_random_A(n, symmetric=False, seed=None, dtype=np.double):
        """Random A matrix; symmetric or general."""
        rng = np.random.default_rng(seed)
        A = rng.standard_normal((n, n),dtype=dtype)
        if symmetric:
            A = 0.5 * (A + A.T)
        return A

    def test_eig_np(A,lamb,V):
        D = np.diag(lamb.reshape((-1)))
        R = A @ V - V @ D
        res_global = sla.norm(R, ord='fro')
        return res_global

    def test_eig2_np(A,B,lamb,V):
        D = np.diag(lamb.reshape((-1)))
        R = A @ V - B @ V @ D
        res_global = sla.norm(R, ord='fro')
        return res_global

    dt_in = np.cdouble       # input data type
    dt_test = np.clongdouble # desired data type

    # A = make_random_A(n, symmetric=symmetric, seed=seed, dtype=dt_in)
    # B = make_random_spd(n, seed=seed, dtype=dt_in)

    A = dt_test(A)
    B = dt_test(B)
    M = sla.solve(B.copy(), A.copy())
    print('M dtype', M.dtype) # sla produces only fp64
    M = dt_test(M)

    tic('eig:numpy')
    #w_ref, V_ref = sla.eig(A, B)  # generalized EVP
    w_np, V_np = sla.eig(M.copy())
    toc('eig:numpy',out=True)

    tic('eig:mpmath')
#    #w_mp, V_mp = mp_eig2(A,B)
    w_mp, V_mp = mp_eig(M)
    toc('eig:mpmath')

#     tic('eig:eigenpy')
#     w_ep, V_ep = eigenpy_eig(M.copy())
# #    w_ep, V_ep = eigenpy_eig2(A.copy(),B.copy())
#     toc('eig:eigenpy')

    tic('eig:ours')
    w_my, V_my = my_eig(M.copy())
    toc('eig:ours',out=True)

    res_np = test_eig_np(M, w_np, V_np)
    res_mp = test_eig_np(M, w_mp, V_mp)
    # res_ep = test_eig_np(M, w_ep, V_ep)
    res_my = test_eig_np(M, w_my, V_my)

    res2_np = test_eig2_np(A, B, w_np, V_np)
    res2_mp = test_eig2_np(A, B, w_mp, V_mp)
    # res2_ep = test_eig2_np(A, B, w_ep, V_ep)
    res2_my = test_eig2_np(A, B, w_my, V_my)
    print("-"*70)
    print("||res||_fro   %-6s %-12s %-12s %-12s" % ("n","numpy","mpmath","ours"))
    print("eig:          %-6d %-12.4e %-12.4e %-12.4e" % (A.shape[0], res_np, res_mp, res_my))
    print("general eig:  %-6d %-12.4e %-12.4e %-12.4e" % (A.shape[0], res2_np, res2_mp, res2_my))
    print("-"*70)


#test(10,seed=10) # < 1 sec
# test(50,seed=10) # 12 sec


