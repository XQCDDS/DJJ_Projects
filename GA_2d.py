import numpy as np
import matplotlib.pyplot as plt
import tkinter.messagebox  # 弹窗库
import timeit

DNA_SIZE = 24  # 编码长度
POP_SIZE = 200  # 种群长度
CROSSOVER_RATE = 0.8  # 交叉概率
MUTATION_RATE = 0.005  # 变异概率
N_GENERATIONS = 50  # 迭代次数
X_BOUND = [-3, 3]  # x取值范围
Y_BOUND = [-3, 3]  # y取值范围
ACCURACY = 100  # 取点精度
FIT_CHOICE = True       # 适应方向
EXPRESSION_2D = "-x ** 2"  # 函数表达式


# 编码
def set_2d_data(x):
    return eval(EXPRESSION_2D.replace("\n", "").replace("\\", "").replace("\"", "").replace("\'", "").lower())


# 解码
def translateDNA(pop):  # pop表示种群矩阵，一行表示一个二进制编码表示的DNA，矩阵的行数为种群长度
    x_pop = pop[::]  # 奇数列表示X
    x_sum = np.zeros(POP_SIZE)
    # 为将二进制串映射到指定范围，首先先将二进制串按权展开，将二进制数转化为十进制数，
    # 然后将转换后的实数压缩,通过以上这些步骤所有二进制串表示都可以转换为[0,1]之间的小数，现在只需要将[0,1] 区间内的数映射到我们要的区间即可。
    for i in range(DNA_SIZE):
        for j in range(POP_SIZE):
            x_sum[j] += x_pop[j][i] * (2 ** (DNA_SIZE - 1 - i))
    x = x_sum / float(2 ** DNA_SIZE - 1) * (X_BOUND[1] - X_BOUND[0]) + X_BOUND[0]

    # dot矩阵乘法
    # x = x_pop.dot(2 ** np.arange(DNA_SIZE)[::-1]) / float(2 ** DNA_SIZE - 1) * (X_BOUND[1] - X_BOUND[0]) + X_BOUND[0]
    return x


# 适应度函数 确定每个个体被保留下来的概率，而概率不能是负值
def get_fitness(pop):
    x = translateDNA(pop)
    y = set_2d_data(x)

    if FIT_CHOICE:
        # 向最大值适应
        # 减去最小的适应度是为了防止适应度出现负数，通过这一步fitness的范围为[0, np.max(z)-np.min(z)],最后在加上一个很小的数防止出现为0的适应度
        return (y - np.min(y)) + 1e-3
    else:
        # 向最小值适应
        return (np.max(y) - y) + 1e-3


# 交叉
def crossover(pop, CROSSOVER_RATE=0.8):
    new_pop = []
    for father in pop:  # 遍历种群中的每一个个体，将该个体作为父亲
        child = father  # 孩子先得到父亲的全部基因（这里我把一串二进制串的那些0，1称为基因）
        if np.random.rand() < CROSSOVER_RATE:  # 产生子代时不是必然发生交叉，而是以一定的概率发生交叉
            mother = pop[np.random.randint(POP_SIZE)]  # 再种群中选择另一个个体，并将该个体作为母亲
            cross_points = np.random.randint(low=0, high=DNA_SIZE)  # 随机产生交叉的点
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
    x = translateDNA(pop)
    print("最优的基因型：", pop[max_fitness_index])
    print("最优解x:", x[max_fitness_index])
    print("y:", set_2d_data(x[max_fitness_index]))

    mes = "最优的基因型: \n" + str(pop[max_fitness_index]) \
          + "\n\n最优解:\nx: " + str(x[max_fitness_index]) \
          + "\ny: " + str(set_2d_data(x[max_fitness_index]))\
          + "\n\n耗时:" + str(t2) + " 秒"
    window = tkinter.Tk()
    window.title('结果')
    # window.geometry("500x300") # 弹窗大小
    tkinter.Label(window,
                  text=mes,  # 标签的文字
                  bg='white',  # 背景颜色
                  font=('黑体', 20),  # 字体和字体大小
                  width=80, height=20  # 标签长宽
                  ).pack()  # 固定窗口位置
    window.mainloop()


# 绘制2d图像
def plot_2d(ax):
    X = np.linspace(*X_BOUND, ACCURACY)  # *列表解引用
    Y = set_2d_data(X)

    ax.plot(X, Y, color="red")

    ax.set_xlabel(u' X 轴')
    ax.set_ylabel(u' y 轴')
    ax.set_xlim(X_BOUND[0], X_BOUND[1])
    ax.set_ylim(Y_BOUND[0], Y_BOUND[1])
    # 用来正常显示中文标签
    plt.rcParams['font.sans-serif'] = ['SimHei']
    # 用来正常显示负号
    plt.rcParams['axes.unicode_minus'] = False
    plt.pause(3)
    plt.show()


def print_2d():
    t1 = timeit.default_timer()

    fig = plt.figure("二维遗传算法(单倍体)")
    ax = fig.add_subplot()
    plt.ion()  # 将画图模式改为交互模式，程序遇到plt.show不会暂停，而是继续执行
    plot_2d(ax)

    # pop表示种群矩阵，一行表示一个二进制编码表示的DNA，矩阵的行数为种群数目,DNA_SIZE为编码长度,模拟由父母生成个体，一半的DNA分别来自父母
    pop = np.random.randint(2, size=(POP_SIZE, DNA_SIZE))

    for _ in range(N_GENERATIONS):  # 迭代N代
        x = translateDNA(pop)
        # 去除之前生成的黑点，生成新的个体
        if 'dna' in locals():
            dna.remove()
        # 生成散点图，数据点与原先的进行画图比较，在曲面上生成个体黑点
        dna = ax.scatter(x, set_2d_data(x), c='black', marker='o')
        plt.show()
        plt.pause(0.1)
        pop = np.array(crossover(pop, CROSSOVER_RATE))

        fitness = get_fitness(pop)
        # 选择生成新的种群
        pop = select(pop, fitness)

    t2 = timeit.default_timer() - t1
    print(u'耗时: {} s'.format(t2))

    print_info(pop, t2)
    plt.ioff()
    plot_2d(ax)


if __name__ == "__main__":
    print_2d()
