import ctypes
import numpy as np
import timeit
import time

G = 6.67430e-11
total_t = ctypes.c_long()

class Vector(ctypes.Structure):
    _fields_ = [("x", ctypes.c_double), ("y", ctypes.c_double)]

_gravity = ctypes.CDLL("./gravitylib.so", winmode=2)
_gravity.acc.argtypes = [
    np.ctypeslib.ndpointer([('x', '<f8'), ('y', '<f8')], 1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(np.double, 1, flags='C_CONTIGUOUS'),
    ctypes.c_int,
    np.ctypeslib.ndpointer([('x', '<f8'), ('y', '<f8')], 1, flags='C_CONTIGUOUS'),
    #ctypes.POINTER(ctypes.c_long)
]

_gravity.acc_x_n.argtypes = [
    np.ctypeslib.ndpointer([('x', '<f8'), ('y', '<f8')], 1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(np.double, 1, flags='C_CONTIGUOUS'),
    ctypes.c_int,
    np.ctypeslib.ndpointer([('x', '<f8'), ('y', '<f8')], 1, flags='C_CONTIGUOUS'),
    #ctypes.POINTER(ctypes.c_long),
    ctypes.c_int
]
_gravity.acc_x_n.restype = ctypes.c_long

def c_acc(x, m, N):
    a = np.empty_like(x)
    _gravity.acc(x, m, N, a)
    return a

def py_acc(x, m):
    r = x[np.newaxis, :, :] - x[:, np.newaxis, :]
    abs_r = np.linalg.norm(r, axis=2)
    np.fill_diagonal(abs_r, 1)
    a_mat = m[np.newaxis, :, np.newaxis] * r / abs_r[:, :, np.newaxis] ** 3
    return G * np.sum(a_mat, axis=1)


def main():
    N = 1000
    # c_x = np.array([(0, 0), (0, 3.85e+8)], dtype=[('x', '<f8'), ('y', '<f8')])
    c_xx = np.random.uniform(-1e+9, 1e+9, N)
    c_xy = np.random.uniform(-1e+9, 1e+9, N)
    c_x = np.array(list(zip(c_xx, c_xy)), dtype=[('x', '<f8'), ('y', '<f8')])
    #print(c_x)
    # v = np.array([(0, 0), (1022, 0)], dtype=[('x', '<f8'), ('y', '<f8')])
    #m = np.array([5.9724e+24, 0.07346e+24])
    m = np.random.uniform(1e+22, 1e+26, N)

    c_t = timeit.timeit("c_acc(c_x, m, 2)", globals={"c_acc": c_acc, "c_x": c_x, "m": m}, number=100_000)
    print("C with python calls:", c_t / 100_000)

    py_x = np.array(list(zip(c_xx, c_xy)))
    py_t = timeit.timeit("py_acc(py_x, m)", globals={"py_acc": py_acc, "py_x": py_x, "m": m}, number=1_00)
    print("Pure python:", py_t / 1_00)

    start = time.time()
    a = np.empty_like(c_x)
    it = _gravity.acc_x_n(c_x, m, 2, a, 100_000)
    end = time.time() - start
    print("Pure C:", end / 100_000, it / 10**11)

    print(total_t.value)

main()