import numpy as np

# 参数
RT = 0.5  # hopping
mu = 0.3  # 化学势
beta = 28  # 逆温度
Nlx = Nly = 10

# 二维正方格点的能带（简单紧束缚）
# ε(k) = -2t(cos(kx) + cos(ky))
kx = 2*np.pi*np.arange(Nlx)/Nlx
ky = 2*np.pi*np.arange(Nly)/Nly
KX, KY = np.meshgrid(kx, ky)
epsilon = -2*RT*(np.cos(KX) + np.cos(KY))

# 费米函数
def fermi(E, mu, beta):
    x = beta*(E - mu)
    # 避免溢出
    return np.where(x > 100, 0.0, np.where(x < -100, 1.0, 1/(1 + np.exp(x))))

# 计算密度（每个自旋）
n_per_spin = np.mean(fermi(epsilon, mu, beta))
n_total = 2 * n_per_spin  # 上下自旋

print(f"参数: RT={RT}, mu={mu}, beta={beta}")
print(f"能带范围: [{epsilon.min():.2f}, {epsilon.max():.2f}]")
print(f"在T=1/β时的密度（每自旋）: {n_per_spin:.6f}")
print(f"总密度（两自旋）: {n_total:.6f}")
print()

# 测试不同的mu值
print("不同μ值的密度:")
for mu_test in [0.0, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0]:
    n = 2 * np.mean(fermi(epsilon, mu_test, beta))
    print(f"  μ={mu_test:.1f}: 密度={n:.6f}")
