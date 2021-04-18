import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import tkinter.messagebox  # 弹窗库
from scipy import optimize as op  # 函数的曲线拟合
from math import sqrt
import timeit

DNA_SIZE = 24  # 基因长度
POP_SIZE = 200  # 种群个数
CROSSOVER_RATE = 0.8  # 交叉概率
MUTATION_RATE = 0.005  # 变异概率
N_GENERATIONS = 50  # 迭代次数
X_BOUND = [-3, 3]  # x取值范围
Y_BOUND = [-3, 3]  # y取值范围
Z_BOUND = [-10, 10]  # z取值范围
ACCURACY = 100  # 取点精度
FIT_CHOICE = True  # 适应方向
EXPRESSION_3D = ""  # 函数表达式
CURVE_POINT_3D = "(2,2,1),(2,0.3,0.15),(2,1,0.5),\
    (1,1,1),(1,2,2),(1,0.3,0.3),\
    (3,2,0.67),(3,1,0.33),(3,0.3,0.1)"  # 3D曲线拟合点


# 编码
def set_3d_data(x, y):
    return eval(EXPRESSION_3D.replace("\n", "").replace("\\", "").replace("\"", "").replace("\'", "").lower())


# 曲线拟合函数
def set_curve_func(X, a, b):
    x = X[0]
    y = X[1]
    z = a * y / x + b
    return z.ravel()  # 多维数组转为一维数组


def get_3d_popt():
    strxyz = eval(CURVE_POINT_3D.replace("\n", ",").replace(",,", ",").replace("\"", "").replace("\'", ""))
    x_group = [strxyz[i][0] for i in range(len(strxyz))]
    y_group = [strxyz[j][1] for j in range(len(strxyz))]
    z_group = [strxyz[k][2] for k in range(len(strxyz))]
    # 格式化数据
    for i in range(len(x_group) - 1):
        for j in range(i + 1, len(x_group)):
            if y_group[i] > y_group[j]:
                x_group[i], x_group[j] = x_group[j], x_group[i]
                y_group[i], y_group[j] = y_group[j], y_group[i]
                z_group[i], z_group[j] = z_group[j], z_group[i]
    for i in range(round(sqrt(len(x_group))) - 1):
        for j in range(i + 1, round(sqrt(len(x_group)))):
            if x_group[i] > x_group[j]:
                x_group[i], x_group[j] = x_group[j], x_group[i]
                y_group[i], y_group[j] = y_group[j], y_group[i]
                z_group[i], z_group[j] = z_group[j], z_group[i]

    # 需要拟合的数据组，转为二维数组
    x_group = np.array(list(set(x_group))).reshape(1, -1)
    y_group = np.array(list(set(y_group))).reshape(1, -1)
    z_group = np.array(z_group)

    # 根据坐标矩阵生成网格点
    X, Y = np.meshgrid(x_group, y_group)
    # 提升到三维数组
    XX = np.expand_dims(X, 0)
    YY = np.expand_dims(Y, 0)
    # 两个三维数组加在一起,Z由X,Y两个点决定，所以要在三维数组里选择二维数组
    XY = np.append(XX, YY, axis=0)

    popt, pcov = op.curve_fit(set_curve_func, XY, z_group)
    return x_group, y_group, z_group, popt, pcov


# 解码
def translateDNA(pop):  # pop表示种群矩阵，一行表示一个二进制编码表示的DNA，矩阵的行数为种群长度
    x_pop = pop[:, 1::2]  # 奇数列表示X
    y_pop = pop[:, ::2]  # 偶数列表示Y
    x_sum = np.zeros(POP_SIZE)
    y_sum = np.zeros(POP_SIZE)

    # 为将二进制串映射到指定范围，首先先将二进制串按权展开，将二进制数转化为十进制数，
    # 然后将转换后的实数压缩,通过以上这些步骤所有二进制串表示都可以转换为[0,1]之间的小数，现在只需要将[0,1] 区间内的数映射到我们要的区间即可。
    for i in range(DNA_SIZE):
        for j in range(POP_SIZE):
            x_sum[j] += x_pop[j][i] * (2 ** (DNA_SIZE - 1 - i))
    x = x_sum / float(2 ** DNA_SIZE - 1) * (X_BOUND[1] - X_BOUND[0]) + X_BOUND[0]

    for i in range(DNA_SIZE):
        for j in range(POP_SIZE):
            y_sum[j] += y_pop[j][i] * (2 ** (DNA_SIZE - 1 - i))
    y = y_sum / float(2 ** DNA_SIZE - 1) * (Y_BOUND[1] - Y_BOUND[0]) + Y_BOUND[0]

    # dot矩阵乘法
    # x = x_pop.dot(2 ** np.arange(DNA_SIZE)[::-1]) / float(2 ** DNA_SIZE - 1) * (X_BOUND[1] - X_BOUND[0]) + X_BOUND[0]
    # y = y_pop.dot(2 ** np.arange(DNA_SIZE)[::-1]) / float(2 ** DNA_SIZE - 1) * (Y_BOUND[1] - Y_BOUND[0]) + Y_BOUND[0]
    return x, y


# 适应度函数 确定每个个体被保留下来的概率，而概率不能是负值
def get_fitness(pop):
    if EXPRESSION_3D:
        x, y = translateDNA(pop)
        z = set_3d_data(x, y)
    else:
        x_group, y_group, z_group, popt, pcov = get_3d_popt()
        X_BOUND[0] = x_group[0][0]
        X_BOUND[1] = x_group[0][len(x_group[0]) - 1]
        Y_BOUND[0] = y_group[0][0]
        Y_BOUND[1] = y_group[0][len(y_group[0]) - 1]
        x, y = translateDNA(pop)
        x = x.reshape(1, -1)
        y = y.reshape(1, -1)
        xy = np.append(x, y, 1).reshape(2, -1)
        z = set_curve_func(xy, *popt)

    if FIT_CHOICE:
        # 向最大值适应
        # 减去最小的适应度是为了防止适应度出现负数，通过这一步fitness的范围为[0, np.max(z)-np.min(z)],最后在加上一个很小的数防止出现为0的适应度
        return (z - np.min(z)) + 1e-3
    else:
        # 向最小值适应
        return (np.max(z) - z) + 1e-3


# 交叉
def crossover(pop, CROSSOVER_RATE=0.8):
    new_pop = []
    for father in pop:  # 遍历种群中的每一个个体，将该个体作为父亲
        child = father  # 孩子先得到父亲的全部基因（这里我把一串二进制串的那些0，1称为基因）
        if np.random.rand() < CROSSOVER_RATE:  # 产生子代时不是必然发生交叉，而是以一定的概率发生交叉
            mother = pop[np.random.randint(POP_SIZE)]  # 再种群中选择另一个个体，并将该个体作为母亲
            cross_points = np.random.randint(low=0, high=DNA_SIZE * 2)  # 随机产生交叉的点
            child[cross_points:] = mother[cross_points:]  # 孩子得到位于交叉点后的母亲的基因
        mutation(child, MUTATION_RATE)  # 每个后代有一定的机率发生变异
        new_pop.append(child)
    return new_pop


# 变异
def mutation(child, MUTATION_RATE=0.005):
    while np.random.rand() < MUTATION_RATE:  # 以MUTATION_RATE的概率进行变异
        mutate_point = np.random.randint(0, DNA_SIZE)  # 随机产生一个实数，代表要变异基因的位置
        child[mutate_point] = child[mutate_point] ^ 1  # 将变异点的二进制为反转


# 选择
def select(pop, fitness):
    # fitness为一个列表，用来规定选取np.arange(POP_SIZE)中每个元素的概率，默认为概率相同，未被选择的会被同化
    idx = np.random.choice(np.arange(POP_SIZE), size=POP_SIZE, replace=True, p=fitness / (fitness.sum()))

    # 选择最大值
    # idx = np.full(POP_SIZE, np.argmax(set_2d_data(translateDNA(pop))))

    # 选择最小值
    # idx = np.full(POP_SIZE, np.argmin(set_2d_data(translateDNA(pop))))
    return pop[idx]


# 打印结果
def print_info(pop, t2):
    fitness = get_fitness(pop)
    # 返回最大适应值的下标
    max_fitness_index = np.argmax(fitness)
    print("最后一代的最大适应值:", fitness[max_fitness_index])
    x, y = translateDNA(pop)
    print("最优的基因型：", pop[max_fitness_index])
    print("最优解(x, y):", (x[max_fitness_index], y[max_fitness_index]))

    if EXPRESSION_3D:
        print("z:", set_3d_data(x[max_fitness_index], y[max_fitness_index]))
        mes = "最优的基因型: \n" + str(pop[max_fitness_index]) \
              + "\n\n最优解: \n(x, y): " + str((x[max_fitness_index], y[max_fitness_index])) \
              + "\nz: " + str(set_3d_data(x[max_fitness_index], y[max_fitness_index]))\
              + "\n\n 结果：在成本为 " + str(x[max_fitness_index]) \
              + " ,效率为 " + str(y[max_fitness_index]) \
              + " 的情况下\n价值最优为 " + str(set_3d_data(x[max_fitness_index], y[max_fitness_index]))\
              + "\n\n耗时:" + str(t2) + " 秒"
    else:
        x_group, y_group, z_group, popt, pcov = get_3d_popt()
        X = x[max_fitness_index].reshape(1, -1)
        Y = y[max_fitness_index].reshape(1, -1)
        # 提升到三维数组
        XX = np.expand_dims(X, 0)
        YY = np.expand_dims(Y, 0)
        # 两个三维数组加在一起,Z由X,Y两个点决定，所以要在三维数组里选择二维数组
        XY = np.append(XX, YY, axis=0)
        print("z:", set_curve_func(XY, *popt))
        mes = "最优的基因型: \n" + str(pop[max_fitness_index]) \
              + "\n\n最优解: \n(x, y): " + str((x[max_fitness_index], y[max_fitness_index])) \
              + "\nz: " + str(set_curve_func(XY, *popt)) \
              + "\n\n 结果：在成本为 " + str(x[max_fitness_index]) \
              + " ,效率为 " + str(y[max_fitness_index]) \
              + " 的情况下\n价值最优为 " + str(set_curve_func(XY, *popt)[0])\
              + "\n\n耗时:" + str(t2) + " 秒"

    window = tkinter.Tk()
    window.title('结果')
    # window.geometry("500x300") # 弹窗大小
    tkinter.Label(window,
                  text=mes,  # 标签的文字
                  bg='white',  # 背景颜色
                  font=('黑体', 20),  # 字体和字体大小
                  width=100, height=20  # 标签长宽
                  ).pack()  # 固定窗口位置
    window.mainloop()


# 绘制3d图像
def plot_3d(ax):
    if EXPRESSION_3D:
        X = np.linspace(*X_BOUND, ACCURACY)  # *列表解引用
        Y = np.linspace(*Y_BOUND, ACCURACY)
        # 两个矩阵分别成为相应的矩阵填充的n阶矩阵的网格点
        X, Y = np.meshgrid(X, Y)
        Z = set_3d_data(X, Y)
        # 绘制一个三维曲面
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm, label=u'自定义函数: y=%s' % EXPRESSION_3D)
        ax.set_xlim(X_BOUND[0], X_BOUND[1])
        ax.set_ylim(Y_BOUND[0], Y_BOUND[1])
        ax.set_zlim(Z_BOUND[0], Z_BOUND[1])
        plt.title(u'自定义函数: z= %s' % EXPRESSION_3D)

    else:
        x_group, y_group, z_group, popt, pcov = get_3d_popt()
        X_BOUND[0] = x_group[0][0]
        X_BOUND[1] = x_group[0][len(x_group[0]) - 1]
        Y_BOUND[0] = y_group[0][0]
        Y_BOUND[1] = y_group[0][len(y_group[0]) - 1]

        X = np.linspace(x_group[0][0], x_group[0][len(x_group[0]) - 1], ACCURACY)
        Y = np.linspace(y_group[0][0], y_group[0][len(y_group[0]) - 1], ACCURACY)
        X, Y = np.meshgrid(X, Y)
        # 提升到三维数组
        XX = np.expand_dims(X, 0)
        YY = np.expand_dims(Y, 0)
        # 两个三维数组加在一起,Z由X,Y两个点决定，所以要在三维数组里选择二维数组
        XY = np.append(XX, YY, axis=0)
        Z = set_curve_func(XY, *popt).reshape(ACCURACY, -1)
        # 生成拟合函数
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm)
        plt.title(u'拟合函数: z=%5.3f * y/x %+5.3f' % tuple(popt))

    ax.set_xlabel(u' 成 本')
    ax.set_ylabel(u' 效 率')
    ax.set_zlabel(u' 价 值')

    # 解决中文显示问题
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    plt.pause(3)
    plt.show()


def print_3d():
    t1 = timeit.default_timer()

    fig = plt.figure("三维APS排程")
    ax = fig.add_subplot(projection='3d')
    plt.ion()  # 将画图模式改为交互模式，程序遇到plt.show不会暂停，而是继续执行
    plot_3d(ax)

    # pop表示种群矩阵，一行表示一个二进制编码表示的DNA，矩阵的行数为种群数目,DNA_SIZE为编码长度,模拟由父母生成个体，一半的DNA分别来自父母
    pop = np.random.randint(2, size=(POP_SIZE, DNA_SIZE * 2))
    for _ in range(N_GENERATIONS):  # 迭代N代
        # 去除之前生成的黑点，生成新的个体
        if 'dna' in locals():
            dna.remove()
        if EXPRESSION_3D:
            x, y = translateDNA(pop)
            # 数据点与原先的进行画图比较，在曲面上生成个体黑点
            dna = ax.scatter(x, y, set_3d_data(x, y), c='black', marker='o')
        else:
            x_group, y_group, z_group, popt, pcov = get_3d_popt()
            X_BOUND[0] = x_group[0][0]
            X_BOUND[1] = x_group[0][len(x_group[0]) - 1]
            Y_BOUND[0] = y_group[0][0]
            Y_BOUND[1] = y_group[0][len(y_group[0]) - 1]
            x, y = translateDNA(pop)
            X = x.reshape(1, -1)
            Y = y.reshape(1, -1)
            # 提升到三维数组
            XX = np.expand_dims(X, 0)
            YY = np.expand_dims(Y, 0)
            # 两个三维数组加在一起,Z由X,Y两个点决定，所以要在三维数组里选择二维数组
            XY = np.append(XX, YY, axis=0)
            dna = ax.scatter(x, y, set_curve_func(XY, *popt), c='black', marker='o', label=u'数据点')

        # 解决中文显示问题
        plt.rcParams['font.sans-serif'] = ['SimHei']
        plt.rcParams['axes.unicode_minus'] = False

        plt.show()
        plt.pause(0.01)
        pop = np.array(crossover(pop, CROSSOVER_RATE))
        fitness = get_fitness(pop)
        # 选择生成新的种群
        pop = select(pop, fitness)

    t2 = timeit.default_timer() - t1
    print(u'耗时: {} s'.format(t2))

    print_info(pop, t2)
    plt.ioff()
    # plot_3d(ax)


if __name__ == "__main__":
    print_3d()
