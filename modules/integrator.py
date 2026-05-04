import numpy as np
from typing import Callable
from debug import pick
import time as pytime

def Euler(y_i: np.ndarray, dt: float, tf: float, acc_function: Callable[[np.ndarray], np.ndarray], time=False):
    """
    Computes position and velocity of the satellite using Euler integrator.
    acc_function can be a simple gravity model f(r) or a full RHS f(t, y).
    """
    start = pytime.process_time()
    N = int(round(tf / dt))
    dim = N + 1
    y = np.zeros((dim, 6))
    y[0, :] = y_i

    # Check if acc_function is a simple gravity model f(r) or full RHS f(t, y)
    import inspect
    sig = inspect.signature(acc_function)
    is_rhs = len(sig.parameters) == 2

    for t_idx in range(N):
        if is_rhs:
            y_dot = acc_function(t_idx * dt, y[t_idx])
        else:
            acc = acc_function(y[t_idx, 0:3])
            y_dot = np.concatenate((y[t_idx, 3:6], acc))
        
        y[t_idx + 1, :] = y[t_idx, :] + y_dot * dt

    elapsed = pytime.process_time() - start
    if time:
        return y, elapsed
    return y

def Heun(y_i: np.ndarray, dt: float, tf: float, acc_function: Callable[[np.ndarray], np.ndarray]):
    N = int(round(tf / dt))
    dim = N + 1
    y = np.zeros((dim, 6))
    y[0, :] = y_i

    import inspect
    sig = inspect.signature(acc_function)
    is_rhs = len(sig.parameters) == 2

    for t_idx in range(N):
        t = t_idx * dt
        if is_rhs:
            y_dot_k1 = acc_function(t, y[t_idx])
            y_predict = y[t_idx] + y_dot_k1 * dt
            y_dot_k2 = acc_function(t + dt, y_predict)
            y[t_idx + 1] = y[t_idx] + 0.5 * dt * (y_dot_k1 + y_dot_k2)
        else:
            acc_k1 = acc_function(y[t_idx, 0:3])
            y_dot_k1 = np.concatenate((y[t_idx, 3:6], acc_k1))
            y_predict = y[t_idx] + y_dot_k1 * dt
            acc_k2 = acc_function(y_predict[0:3])
            y_dot_k2 = np.concatenate((y_predict[3:6], acc_k2))
            y[t_idx + 1] = y[t_idx] + 0.5 * dt * (y_dot_k1 + y_dot_k2)

    return y

def RK4(y_i: np.ndarray, dt: float, tf: float, acc_function: Callable[[np.ndarray], np.ndarray]):
    N = int(round(tf / dt))
    dim = N + 1
    y = np.zeros((dim, 6))
    y[0, :] = y_i

    import inspect
    sig = inspect.signature(acc_function)
    is_rhs = len(sig.parameters) == 2

    for i in range(N):
        t = i * dt
        if is_rhs:
            k1 = acc_function(t, y[i])
            k2 = acc_function(t + 0.5*dt, y[i] + 0.5*dt*k1)
            k3 = acc_function(t + 0.5*dt, y[i] + 0.5*dt*k2)
            k4 = acc_function(t + dt, y[i] + dt*k3)
        else:
            r, v = y[i, 0:3], y[i, 3:6]
            k1 = np.concatenate((v, acc_function(r)))
            
            r2, v2 = r + 0.5*dt*k1[:3], v + 0.5*dt*k1[3:]
            k2 = np.concatenate((v2, acc_function(r2)))
            
            r3, v3 = r + 0.5*dt*k2[:3], v + 0.5*dt*k2[3:]
            k3 = np.concatenate((v3, acc_function(r3)))
            
            r4, v4 = r + dt*k3[:3], v + dt*k3[3:]
            k4 = np.concatenate((v4, acc_function(r4)))

        y[i + 1] = y[i] + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    return y
